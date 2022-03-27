#!/usr/bin/env python3

import argparse
import sys

from itertools import permutations

import logging

from .graph import Graph, DiGraph, load_graph
from .pattern import Pattern
from .helpers import vertex_str, decode_vertices, suborder

log = logging.getLogger("mandoline")


class TD:
    @staticmethod
    def decompose(G, order, extra_edges=None):
        order = list(order)
        assert set(G) == set(order), "Order {} is not incompatible with vertices {}".format(order, list(G))

        GG = G.copy()
        if extra_edges != None:
            GG.add_edges(extra_edges)
        return TD._decompose_rec(G, GG, order, [])

    @staticmethod
    def _decompose_rec(G, subG, order, sep):
        sep = list(sep)
        depth = len(sep)
        R = Graph.from_graph(subG)

        in_neighbours = []
        for v in order:
            if v not in R:
                continue

            Nv = sorted([sep.index(u) for u in G.neighbours(v) if u in sep])
            in_neighbours.append(tuple(Nv))
            sep.append(v)
            R.remove_node(v)
            if not R.is_connected():
                break

        children = []
        for CC in R.connected_components():
            children.append(TD._decompose_rec(G, CC, order, tuple(sep)))

        res = TD(sep, in_neighbours, children, depth)
        return res

    @staticmethod
    def from_string(s):
        res, _ = TD._from_string_rec(s, tuple(), 0)
        return res

    @staticmethod
    def _split_child_string(s):
        if len(s) == 0:
            return []

        # Splits a string by the '|' symbole but ignore portions
        # that are nested inside curly braces.
        splits = [-1]
        depth = 0
        for i,c in enumerate(s):
            if c == '{':
                depth += 1
            elif c == '}':
                depth -= 1
            elif c == '|' and depth == 0:
                splits.append(i)
        splits.append(len(s))
        res = [s[i+1:j] for i,j in zip(splits[:-1], splits[1:])]
        return res

    @staticmethod
    def _from_string_rec(s, sep, start_index, rec_depth=0):
        # Split string into separator / children
        i = s.find('{')
        if i == -1:
            in_neighbours, child_str = s, ""
        else:
            in_neighbours, child_str = s[:i], s[i+1:-1]

        in_neighbours = list(map(decode_vertices, in_neighbours.split(',')))
        in_neighbours = [tuple(N) for N in in_neighbours]
        bag = tuple(range(start_index, start_index+len(in_neighbours)))
        child_index = start_index + len(bag) # Next free index for children
        depth = len(sep)
        sep = sep + bag

        assert bag == sep[depth:]

        children = []
        for c in TD._split_child_string(child_str):
            child, max_index = TD._from_string_rec(c, sep, child_index, rec_depth+1)
            children.append(child)
            child_index = max_index

        return TD(sep, in_neighbours, children, depth), child_index


    def __init__(self, sep, in_neighbours, children, depth):
        self.parent = None
        self._sep = tuple(sep)
        self._bag = tuple(sep[depth:])
        self._in = tuple(in_neighbours)
        self.depth = depth
        self.children = tuple(sorted(children, key=lambda c: c.max_postfix()))
        for c in self.children:
            c.parent = self

    def nodes(self):
        res = set(self._sep)
        for c in self.children:
            res |= c.nodes()
        return res

    def height(self):
        if len(self.children) == 0:
            return len(self._sep)

        return max([c.height() for c in self.children])

    def depth_of(self, x):
        for c in self.children:
            try:
                return c.depth_of(x)
            except ValueError:
                pass

        # This raises a ValueError if x is not in _bag
        return self.depth + self._bag.index(x)

    def upwards_closure(self, nodes):
        """
            Returns all nodes of the decomposition that lie above
            'nodes', e.g. the union of all root-paths leading to each
            member of 'nodes'. The result will contain all members of 'nodes'
            that are contained in this decomposition.
        """
        res = set()
        for i,x in enumerate(self._bag):
            if x in nodes:
                ix = i + self.depth
                res = set(self._sep[:ix+1]) # We want to return thus, but there might be
                                            # yet another member of 'nodes' below us.

        for c in self.children:
            res |= c.upwards_closure(nodes)
        return res

    def __iter__(self):
        return iter(self.nodes())

    def in_neighbours(self, u):
        if u in self._bag:
            i = self._bag.index(u)
            return set([self._sep[j] for j in self._in[i]])

        for c in self.children:
            res = c.in_neighbours(u)
            if res != None:
                return res
        return None

    def descendants(self):
        res = set(self._bag)
        for c in self.children:
            res |= c.descendants()
        return res

    def max_postfix(self):
        if len(self.children) == 0:
            return self._in
        else:
            return self._in + self.children[-1].max_postfix()

    def prefix(self):
        raise RuntimeError("Deprecated")
        return tuple(self._sep[:self.depth])

    def adhesion_size(self):
        """
            Returns the length of the root-path up to the
            first branching vertex.
        """
        return len(self._sep)

    def compatible_with(self, order, level=0):
        prefix = "  "*level
        # print("{}Testing {} in {}".format(prefix, order, self))
        order = suborder(order, self.descendants())
        # print("{}Restricted to {}".format(prefix, order))

        if len(order) == 0:
            return True

        # Match nodes in current bag
        for x in self._bag:
            if len(order) == 0:
                return True
            if x == order[0]:
                order = order[1:]

        # print("{}Matched all but {}".format(prefix, order))

        if len(set(order) & set(self._bag)) != 0:
            # Could not match a vertex from the current bag,
            # so order is incompatible.
            return False

        if len(order) == 0:
            return True

        # Order must be compatbile with _all_ children
        for c in self.children:
            if not c.compatible_with(order, level=level+1):
                return False
        return True

    def to_graph(self):
        """
            Returns the graph described by this decompositon.
        """
        res = Graph()
        res.add_nodes(self._sep)

        # Add upward edges from bag vertices to root-path
        for u, iN in zip(self._bag,self._in):
            N = [self._sep[i] for i in iN] # Translate from indices to vertices
            for v in N:
                res.add_edge(u,v)

        # Recurse
        for c in self.children:
            H = c.to_graph()
            res.add_nodes(H)
            res.add_edges(H.edges())
        return res

    def to_ditree(self):
        """
            Returns the tree underlying this decompositon
            as a directed graph. Arcs are oriented _away_ from
            the root, towards the leaves.
        """
        res = DiGraph()
        for u,v in zip(self._bag[:-1], self._bag[1:]):
            res.add_arc(u,v)

        for c in self.children:
            H = c.to_ditree()
            res.add_arcs(H.arcs())
            res.add_arc(self._bag[-1], c._bag[0])
        return res

    def to_piece(self, pattern_size):
        """
            Turns a _linear_ TD decomposition into a piece;
            this enables us to use the counting infrastructure of
            the enumeration part.
        """
        assert self.is_linear()
        order = list(self.orders())[0]
        LG, _ = self.to_graph().to_lgraph(order)

        # This TD decomposition might describe a disconnected graph.
        # We do know, however, that the final pattern is connected; hence
        # every vertex that is not wreachable inside this decomposition
        # _must_ be wreachable in order to be part of the pattern.
        # The 'pattern_size' here informs the Pattern construction what
        # the maximum possible wr-distance is for such disconnected parts.
        pattern = Pattern.from_disconnected(LG, pattern_size)

        piece = pattern.to_single_piece()
        return piece

    def orders(self):
        for o in self.to_ditree().embeddings():
            yield o

    def suborders(self, vertices):
        for o in self.to_ditree().subembeddings(vertices):
            yield o

    def copy(self):
        copy_children = []
        for c in self.children:
            copy_children.append(c.copy())

        res = TD(self._sep, self._in, copy_children, self.depth)
        for c in res.children:
            assert(c.parent == res)
        return res

    def __hash__(self):
        res = 219787954134**(self.depth+1)
        for N in self._in:
            res += hash(N)
        for c in self.children:
            res += hash(c)
        return res

    def __eq__(self, other):
        if other == None:
            return False

        if self.depth != other.depth or len(self.children) != len(other.children):
            return False

        if self._in != other._in:
            return False

        for child, other_child in zip(self.children, other.children):
            if child != other_child:
                return False

        return True

    def is_linear(self):
        assert(len(self.children) != 1)
        return len(self.children) == 0

    def split(self):
        """
            Splits the current decomposition along the common prefix
            of its children: we return every child-decomposition with
            the common root prefix attached to it.
        """
        if len(self.children) == 0:
            yield self.copy()
            return

        assert(len(self.children) != 1)

        for c in self.children:
            res = c.copy()
            res._in = self._in + res._in
            res.depth = self.depth
            res._bag = res._sep[res.depth:]
            res.parent = self.parent

            yield res

    def chop(self, size):
        assert size >= 1 and size <= len(self._in)

        child = self.copy()

        if size == len(self._in):
            return child

        in_neighbours = self._in[:size]
        sep = self._sep[:size]

        child.depth += size
        child._bag = child._sep[child.depth:]
        child._in = child._in[size:]

        return TD(sep, in_neighbours, [child], self.depth)

    def merge(self, other, merge_depth):
        assert(self.depth == other.depth)
        assert(self._in[:merge_depth] == other._in[:merge_depth])

        left = self.chop(merge_depth)
        right = other.chop(merge_depth)

        res = TD(left._sep, left._in, left.children + right.children, left.depth)
        assert(len(res.children) != 1)

        return res

    def overlap_with(self, other):
        """
            Computes three sets of nodes (L,R,U): 
            L contains nodes exclusive to `self`,
            R nodes exclusive to `other` and U contains all
            joint nodes.
        """
        rootPath = set(self._sep) & set(other._sep)
        nodesA, nodesB = self.nodes(), other.nodes()
        common_nodes = nodesA & nodesB
        nodesA -= common_nodes
        nodesB -= common_nodes
        nodesAll = nodesA | nodesB | rootPath
        return (nodesA, nodesB, nodesAll)        

    def __repr__(self):
        return self.order_string()

    def td_string(self):
        if len(self.children) == 0:
            return ','.join(map(lambda N: vertex_str(N), self._in))
        else:
            return ','.join(map(lambda N: vertex_str(N), self._in)) + '{' + '|'.join(map(lambda c: c.td_string(), self.children)) + '}'

    def order_string(self):
        if len(self.children) == 0:
            return ''.join(map(str, self._bag))
        else:
            return ''.join(map(str, self._bag)) + '{' + ','.join(map(lambda c: c.order_string(), self.children)) + '}'

def main():
    parser = argparse.ArgumentParser(description='Exhaustively decomposes H into tree decompositions')

    parser.add_argument('H', help='Pattern graph H')
    # parser.add_argument('G', help='Host graph G')
    parser.add_argument('--debug', action='store_true')
    parser.add_argument('--no-reduction', action='store_true' )
    parser.add_argument('--quiet', action='store_true' )

    args = parser.parse_args()

    # Set up logging
    ch = logging.StreamHandler(sys.stdout)
    if args.quiet:
        # Mute on top level, don't add handler
        log.setLevel(logging.CRITICAL)
    elif args.debug:
        ch.setLevel(logging.DEBUG)
        log.setLevel(logging.DEBUG)
        log.addHandler(ch)
    else:
        ch.setLevel(logging.INFO)
        log.setLevel(logging.INFO)
        log.addHandler(ch)

    # Load pattern and graph
    H = load_graph(args.H)
    log.info("Loaded pattern graph with {} vertices and {} edges".format(len(H), H.num_edges()))
    log.info(H)

    seen = set()
    for order in permutations(H):
        tdH = TD.decompose(H, order)
        if tdH in seen:
            continue
        seen.add(tdH)

        print("\nOrder:", ''.join(map(str,order)))
        print(H.to_lgraph(order)[0])
        print(tdH.order_string())

        print("Decomposition", tdH)
        print("Graph", tdH.td_string())

        print("Represented orders:")
        for o in tdH.orders():
            print(" ", ''.join(map(str,o)))

        splits = list(tdH.split())
        if len(splits) == 0:
            continue

        print("Splits:")
        for td in splits:
            print("  ", td, td.nodes())

        merge_depth = len(tdH._bag)
        print("Merging back at depth {}:".format(merge_depth))
        merged = splits[0]
        print("  ", merged)
        for td in splits[1:]:
            merged = merged.merge(td, merge_depth)
            print("  ", merged)
        assert(merged == tdH)

    print("\n")
    print("Computed {} tree decompositions".format(len(seen)))


if __name__ == "__main__":
    main()
