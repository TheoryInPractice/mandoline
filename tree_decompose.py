#!/usr/bin/env python3

from graph import Graph, load_graph, short_str
from pattern import PatternBuilder, Pattern

import argparse
import itertools
import sys

from collections import defaultdict
from sortedcontainers import SortedSet
from itertools import permutations, product, chain
import bisect
import math, random
import cairo

import logging

log = logging.getLogger("mandoline")

class TD:
    @staticmethod
    def decompose(G, order):
        assert(G.is_connected())
        assert set(G) == set(order), "Order {} is not incompatible with vertices {}".format(order, list(G)) 
        return TD._decompose_rec(G, G, order, [])

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

    def __init__(self, sep, in_neighbours, children, depth):
        self.parent = None
        self._sep = tuple(sep)
        self._bag = tuple(sep[depth:])
        self.in_neighbours = tuple(in_neighbours)
        self.depth = depth
        self.children = tuple(sorted(children, key=lambda c: c.in_neighbours))
        for c in self.children:
            c.parent = self

    def nodes(self):
        res = set(self._sep) 
        for c in self.children:
            res |= c.nodes()
        return res

    def prefix(self):
        return tuple(self._sep[:self.depth])

    def orders(self):
        res = tuple(self._bag)
        if len(self.children) == 0:
            yield res
            return

        child_orders = [list(c.orders()) for c in self.children]

        for perm in permutations(child_orders):
            for prod in product(*perm):
                yield res + tuple(chain(*prod))

    def copy(self):
        copy_children = []
        for c in self.children:
            copy_children.append(c.copy())

        res = TD(self._sep, self.in_neighbours, copy_children, self.depth)
        for c in res.children:
            assert(c.parent == res)
        return res

    def __hash__(self):
        res = 219787954134**(self.depth+1)
        for N in self.in_neighbours:
            res += hash(N)
        for c in self.children:
            res += hash(c)
        return res

    def __eq__(self, other):
        if other == None:
            return False

        if self.depth != other.depth or len(self.children) != len(other.children):
            return False

        if self.in_neighbours != other.in_neighbours:
            return False

        for child, other_child in zip(self.children, other.children):
            if child != other_child:
                return False

        return True

    def adhesion_size(self):
        return len(self.in_neighbours)

    def split(self):
        if len(self.children) == 0:
            yield self.copy() 
            return

        for c in self.children:
            res = c.copy()
            res.in_neighbours = self.in_neighbours + res.in_neighbours
            res.depth = self.depth
            res._bag = res._sep[res.depth:]
            res.parent = self.parent

            yield res

    def chop(self, size):
        assert size >= 1 and size <= len(self.in_neighbours)

        child = self.copy()

        if size == len(self.in_neighbours):
            return child

        in_neighbours = self.in_neighbours[:size]
        sep = self._sep[:size]

        child.depth += size
        child._bag = child._sep[child.depth:]
        child.in_neighbours = child.in_neighbours[size:]

        return TD(sep, in_neighbours, [child], self.depth)

    def merge(self, other, merge_depth):
        assert(self.depth == other.depth)
        assert(self.in_neighbours[:merge_depth] == other.in_neighbours[:merge_depth])

        left = self.chop(merge_depth)
        right = other.chop(merge_depth)

        res = TD(left._sep, left.in_neighbours, left.children + right.children, left.depth)

        return res

    def __repr__(self):
        return self.order_string()

    def td_string(self):
        if len(self.children) == 0:
            return ','.join(map(lambda N: short_str(N), self.in_neighbours))
        else:
            return ','.join(map(lambda N: short_str(N), self.in_neighbours)) + '{' + '|'.join(map(lambda c: c.td_string(), self.children)) + '}'       

    def order_string(self):
        if len(self.children) == 0:
            return ''.join(map(str, self._bag))
        else:
            return ''.join(map(str, self._bag)) + '{' + ','.join(map(lambda c: c.order_string(), self.children)) + '}'

if __name__ == "__main__":
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