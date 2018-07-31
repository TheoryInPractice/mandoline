#!/usr/bin/env python3
import itertools
from collections import defaultdict as defaultdict
import gzip
from operator import itemgetter
from itertools import chain, combinations, permutations, product

from sortedcontainers import SortedSet
from datastructures import Indexmap

from helpers import singlehash, pairhash, short_str

import bisect
import math, random
import unittest
import logging

log = logging.getLogger("mandoline")

class DAGError(Exception):
    """
        Thrown by top. embedding algorithms if underlying
        digraph is not a DAG.
    """
    pass

def load_graph(filename):
    import os.path
    res = Graph()

    _, ext = os.path.splitext(filename)

    if ext == '.gz':
        with gzip.open(filename, 'r') as filebuf:
            for l in filebuf:
                u, v = l.decode().split()
                u, v = int(u), int(v)
                if u == v:
                    continue
                res.add_edge(u,v)
    else:
        with open(filename, 'r') as filebuf:
            for l in filebuf:
                u, v = l.split()
                u, v = int(u), int(v)
                if u == v:
                    continue
                res.add_edge(u,v)
    res.remove_loops()
    return res


class Graph:
    def __init__(self):
        self.adj = defaultdict(set)
        self.nodes = set()
        self.hash = 14695981039346656037

    @staticmethod
    def from_graph(other):
        res = Graph()
        for v in other:
            res.add_node(v)
        for u,v in other.edges():
            res.add_edge(u,v)
        return res

    def copy(self):
        return Graph.from_graph(self)

    def __contains__(self, u):
        return u in self.nodes

    def __iter__(self):
        return iter(self.nodes)

    def __len__(self):
        return len(self.nodes)

    def __repr__(self):
        return short_str(sorted(self.nodes)) + "{" + ','.join(map(lambda e: '{}{}'.format(*e), self.edges())) + "}"
    
    def __hash__(self):
        return self.hash

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False

        if self.hash != other.hash:
            return False
        if self.nodes != other.nodes:
            return False
        if set(self.edges()) != set(other.edges()):
            return False

        return True

    def get_max_id(self):
        return max(self.nodes)

    def edges(self):
        for u in self:
            for v in self.adj[u]:
                if u <= v:
                    yield (u,v)

    def num_edges(self):
        return sum(1 for _ in self.edges())

    def add_node(self,u):
        if u in self.nodes:
            return # Important for keeping hash consistent!
        self.nodes.add(u)
        self.hash ^= singlehash(u)

    def add_nodes(self, nodes):
        for u in nodes:
            self.add_node(u) # Safest way to keep the hash intact

    def remove_node(self,u):
        if u not in self.nodes:
            return # Important for keeping hash consistent!
        for v in self.neighbours(u):
            self.adj[v].remove(u)
        del self.adj[u]
        self.nodes.remove(u)
        self.hash ^= singlehash(u)

    def add_edge(self,u,v):
        if self.adjacent(u,v):
            return # Important for keeping hash consistent! 
        self.add_node(u)
        self.add_node(v)        
        self.adj[u].add(v)
        self.adj[v].add(u)
        self.hash ^= pairhash(u,v)
        self.hash ^= pairhash(v,u)

    def add_edges(self, edges):
        for u,v in edges:
            self.add_edge(u,v)

    def remove_edge(self,u,v):
        if not self.adjacent(u,v):
            return # Important for keeping hash consistent!
        self.adj[u].discard(v)
        self.adj[v].discard(u)
        self.hash ^= pairhash(u,v)
        self.hash ^= pairhash(v,u)        

    def merge(self, u, v):
        """ 
            Merges node u _onto_ v, meaning we identify the two
            nodes and label the resulting node 'v'.
        """
        assert u in self.nodes and v in self.nodes
        neighs = self.neighbours(u)
        self.remove_node(u)

        for w in neighs:
            self.add_edge(v,w)

    def merge_pairs(self, pairs):
        for u,v in pairs:
            self.merge(u,v)

    def remove_loops(self):
        for v in self:
            self.remove_edge(v,v)

    def adjacent(self,u,v):
        return v in self.adj[u]

    def neighbours(self,u):
        return frozenset(self.adj[u])

    def neighbours_set(self, centers):
        res = set()
        for v in centers:
            res = res | self.neighbours(v)
        return (res - centers)

    def neighbours_set_closed(self, centers):
        res = set()
        for v in centers:
            res = res | self.neighbours(v)
        return res        

    def rneighbours(self,u,r):
        res = set([u])
        for _ in range(r):
            res |= self.neighbours_set_closed(res)
        return res

    def degree(self,u):
        return len(self.adj[u])

    def degree_sequence(self):
        return [ self.degree(u) for u in self.nodes] 

    def degree_dist(self):
        res = defaultdict(int)
        for u in self.nodes:
            res[self.degree(u)] += 1
        return res

    def calc_average_degree(self):
        num_edges = len( [e for e in self.edges()] )
        return float(2*num_edges) / len(self.nodes)

    def subgraph(self, vertices):
        res = Graph()
        selected = set(vertices)
        for v in self:
            if v in selected:
                res.add_node(v)

        for u,v in self.edges():
            if u in selected and v in selected:
                res.add_edge(u,v)
        return res

    def relabel(self, mapping):
        res = Graph()
        for v in self:
            res.add_node(mapping[v])

        for u,v in self.edges():
            res.add_edge(mapping[u], mapping[v])
        return res

    def to_lgraph(self, order=None):
        """
            Returns an ordered graph with node indices from 0 to n-1 and 
            an Indexmap to recover the original ids.
            The nodes are sorted by (degree, label) in descending order
            in order to keep the wcol number small.
        """
        lgraph = LGraph()

        if order == None:
            order = [(self.degree(u), u) for u in self.nodes] 
            order.sort(reverse=True)
            order = [u for _,u in order]

        mapping = Indexmap(len(order))
        for iu,u  in enumerate(order):
            iN = mapping.indices_of(self.neighbours(u))
            iN = [x for x in iN if x != None] # Remove 'None' values
            lgraph._add_node(iu, iN)
            mapping.put(iu, u)

        return lgraph, mapping

    def compute_core(self, minDeg):
        smallDegree = set()
        largeDegree = set(self.nodes)
        changed = True
        while changed:
            changed = False
            for v in largeDegree:
                remaining_N = set(self.neighbours(v)) & largeDegree
                if len(remaining_N) < minDeg:
                    smallDegree.add(v)
                    changed = True
            largeDegree -= smallDegree
        largeDegree = self.nodes - smallDegree
        return self.subgraph(largeDegree)

    def enum_patterns(self):
        found = set()
        for order in permutations(sorted(self)):
            pat, Indexmap = self.to_pattern(order)
            if pat in found:
                continue
            found.add(pat)
            yield pat, Indexmap

    def to_pattern(self, order):
        from pattern import Pattern
        """
            Returns a `pattern' (and ordered graph with node indices from 0 to n-1)
            and an Indexmap to recover the original ids. 
        """
        pat = LGraph()

        mapping = Indexmap(len(order))
        for iu,u  in enumerate(order):
            iN = mapping.indices_of(self.neighbours(u))
            iN = [x for x in iN if x != None] # Remove 'None' values
            pat._add_node(iu, iN)
            mapping.put(iu, u)
        pat.compute_wr(len(pat))
        return Pattern(pat), mapping

    def is_connected(self):
        if len(self.nodes) <= 1:
            return True
        root = next(iter(self.nodes))
        return len(self.bfs(root)) == len(self.nodes)

    def connected_components(self): 
        n = len(self)
        res = []
        if n == 0:
            return res

        remainder = set(self.nodes)
        while len(remainder) > 0:
            root = next(iter(remainder))
            comp_nodes = self.bfs(root)
            component = Graph()
            res.append(self.subgraph(comp_nodes))
            remainder -= comp_nodes
        return res

    def subgraph(self, nodes):
        nodes = set(nodes)
        res = Graph()
        for v in nodes:
            res.add_node(v)
            for u in self.neighbours(v) & nodes:
                res.add_edge(u,v)
        return res

    def bfs(self, root):
        res = set()
        for layer in self.bfs_layers(root):
            res |= layer
        return res

    def bfs_layers(self, root):
        curr = set([root])
        yield frozenset([root])

        seen = set()
        while True:
            seen |= curr
            _next = set()
            for u in curr:
                _next |= set(self.neighbours(u))
            _next -= seen
            if len(_next) == 0:
                break
            curr = _next
            yield frozenset(curr)

    @staticmethod
    def path(n):
        res = Graph()
        for u in range(n):
            res.add_node(u)
        for u,v in zip(range(n-1), range(1,n)):
            res.add_edge(u,v)
        return res

class DiGraph:
    def __init__(self):
        self.out = defaultdict(set)
        self._in = defaultdict(set)
        self.nodes = set()

    @staticmethod
    def from_digraph(other):
        res = DiGraph()
        for v in other:
            res.add_node(v)
        for u,v in other.arcs():
            res.add_arc(u,v)
        return res

    def copy(self):
        return DiGraph.from_digraph(self)

    def __contains__(self, u):
        return u in self.nodes

    def __iter__(self):
        return iter(self.nodes)

    def __len__(self):
        return len(self.nodes)

    def __repr__(self):
        return short_str(sorted(self.nodes)) + "{" + ','.join(map(lambda e: '{}{}'.format(*e), self.arcs())) + "}"
    
    def arcs(self):
        for u in self:
            for v in self.out[u]:
                yield (u,v)

    def num_arcs(self):
        return sum(1 for _ in self.arcs())

    def add_node(self,u):
        self.nodes.add(u)

    def remove_node(self,u):
        if u not in self.nodes:
            return
        for v in self.out[u]:
            self._in[v].remove(u)
        for v in self._in[u]:
            self.out[v].remove(u)
        del self.out[u]
        del self._in[u]
        self.nodes.remove(u)

    def dissolve_node(self,u):
        """
            Removes a node an connects all its
            former in-neighbours to all its former
            out-neighbours. 
        """
        if u not in self.nodes:
            return
        in_neighs = self.in_neighbours(u)
        out_neighs = self.out_neighbours(u)
        self.remove_node(u)

        for v,w in product(in_neighs, out_neighs):
            self.add_arc(v,w)

    def merge(self, u, v):
        """ 
            Merges node u _onto_ v, meaning we identify the two
            nodes and label the resulting node 'v'.
        """
        assert u in self.nodes and v in self.nodes
        in_neighs = self.in_neighbours(u)
        out_neighs = self.out_neighbours(u)
        self.remove_node(u)

        for w in in_neighs:
            self.add_arc(w,v)
        for w in out_neighs:
            self.add_arc(v,w)

    def merge_pairs(self, pairs):
        for u,v in pairs:
            self.merge(u,v)

    def add_arc(self,u,v):
        self.nodes.add(u)
        self.nodes.add(v)        
        self.out[u].add(v)
        self._in[v].add(u)

    def add_arcs(self, arcs):
        for u,v in arcs:
            self.add_arc(u,v)

    def remove_arc(self,u,v):
        self.out[u].discard(v)
        self._in[v].discard(u)

    def remove_loops(self):
        for v in self:
            self.remove_arc(v,v)

    def adjacent(self,u,v):
        return v in self.out[u]

    def out_neighbours(self,u):
        return frozenset(self.out[u])

    def out_neighbours_set(self, centers):
        res = set()
        for u in centers:
            res = res | self.out_neighbours(u)
        return (res - centers)

    def in_neighbours(self,u):
        return frozenset(self._in[u])

    def out_degree(self, u):
        return len(self.out[u])

    def in_degree(self, u):
        return len(self._in[u])       

    def degree(self, u):
        return self.in_degree(u) + self.out_degree(u)   

    def normalize(self):
        """
            Returns a digraph on vertices [0,n-1] and an Indexmap. Moreover,
            if this digraph is a DAG, the vertices will be in a topological
            order (sources first).
        """
        n = len(self)
        imap = Indexmap(n)
        in_degrees = dict((u,self.in_degree(u)) for u in self)
        buckets = defaultdict(set)
        for u, indeg in in_degrees.items():
            buckets[indeg].add(u)

        index = 0
        while index < n:
            # Find a lowest-degree vertex. 
            d = 0
            while len(buckets[d]) == 0:
                d += 1

            # 'Delete' v
            v = buckets[d].pop()
            del in_degrees[v]
            imap.put(index, v)

            # Update in-degrees of remaining vertices
            for u in self.out_neighbours(v):
                if u not in in_degrees:
                    continue # Already processed
                buckets[in_degrees[u]].remove(u)
                in_degrees[u] -= 1
                buckets[in_degrees[u]].add(u)

            index += 1

        # Construct resulting graph
        H = DiGraph()
        for iv in range(n):
            H.add_node(iv)

        for u,v in self.arcs():
            iu, iv = imap.index_of(u), imap.index_of(v)
            H.add_arc(iu, iv)

        return H, imap

    def subembeddings(self, selection):
        """
            Returns all orders (embeddings) of the vertex subset 'selection'
            which appear as suborders in topological orderings of the whole
            graph.
        """
        selection = set(selection)
        H = self.copy()
        remove = self.nodes - selection
        for u in remove:
            H.dissolve_node(u)

        for emb in H.embeddings():
            yield emb

    def embeddings(self):
        """
            Varol-Rotem algorithm to enumerate all topological embeddings of a poset.
            It takes as input the oriented arcs of a digraph.
            Yaakov L. Varol and Doron Rotem, An Algorithm to Generate All Topological Sorting Arrangements.
                Computer J., 24 (1981) pp. 83-84. row
            Adapted from https://github.com/dbasden/python-digraphtools/blob/master/digraphtools/topsort.py.
        """
        H, imap = self.normalize()
        n = len(H)
        arcs = list(H.arcs())

        # Compute (undirected) adjacency matrix. 
        adj = [[False for j in range(n)]+[True] for i in range(n+1)]
        for i,j in arcs:
            if i > j:
                raise DAGError()
            adj[i][j] = True
            adj[j][i] = True

        # Enumerate topological embeddings
        loc = list(range(n))
        p = list(range(n+1))
        yield imap.vertices_at(p[0:n])
        i = 0
        k = 0
        while i < n-1:
            k = loc[i]
            kk = k + 1
            if adj[i][p[kk]]:
                p[i:k+1] = [p[k]]+p[i:k] # roll-right a[i:k+1]
                loc[i] = i
                i += 1
            else:
                p[k],p[kk] = p[kk],p[k]
                loc[i] = kk
                i = 0
                yield imap.vertices_at(p[0:n])

    def is_acyclic(self):
        """
            Tests whether the graph is acyclic by attempting
            to compute a topological embedding.
        """
        H, imap = self.normalize()
        for i,j in H.arcs():
            if i > j:
                return False
        return True


class LGraph:
    """
        A graph imbued with a linear order. Optimized
        fro wcol computation.
    """
    def __init__(self):
        self.wr = [[]]
        self.outN = [] 

    def __len__(self):
        return len(self.wr[0])

    def __iter__(self):
        return iter(range(len(self)))

    def depth(self):
        return len(self.wr)

    def in_neighbours(self, iu):
        return self.wr[0][iu]

    def common_in_neighbours(self, vertices):
        if len(vertices) == 0:
            return []
        if len(vertices) == 1:
            return self.wr[0][vertices[0]]

        # Find vertex with smallest in-neighbourhood  
        temp = [(len(self.wr[0][iv]),iv) for iv in vertices]
        smallest = min(temp)[1]
        smallN = self.wr[0][smallest]

        # Check for each vertex in smallN whether it is connected
        # to all supplied vertices.
        res = []
        for iu in smallN:
            for iv in vertices:
                if not self.adjacent(iu, iv):
                    break
            else: # Did not break
                res.append(iu)
        return res

    def adjacent(self, iu, iv):
        if iv < iu:
            return iv in self.wr[0][iu]
        else:
            return iu in self.wr[0][iv]

    def adjacent_ordered(self, iu, iv):
        # Assume that iu < iv!
        return iu in self.wr[0][iv]

    def edges(self):
        for u in self:
            for v in self.wr[0][u]:
                yield (v,u)

    def wreach_iter(self):
        for iu in self:
            yield iu, sorted(self.wreach_all(iu))
    
    def wreach_union(self, iu, min_dist, max_dist):
        res = []
        for d in range(min_dist, max_dist+1):
            res.extend(self.wr[d][iu])
        return res

    def wreach_all(self, iu):
        res = []
        for wr in self.wr:
            res.extend(wr[iu])
        return res

    def _add_node(self, iu, iN):
        iN = SortedSet(iN)
        self.wr[0].append(iN)
        self.outN.append(SortedSet())

        for iv in iN:
            self.outN[iv].add(iu)

    def forward_bfs(self, iroot, r):
        curr = set([iroot])
        seen = set()
        for d in range(1,r+1):
            seen |= curr
            _next = set()
            for iu in curr:
                for iv in chain(self.wr[0][iu], self.outN[iu]):
                    if iv > iroot:
                        _next.add(iv)
            _next -= seen
            curr = _next
            yield d, frozenset(curr)

    def compute_wr(self, depth):
        # Initialize wr-arrays
        old_depth = len(self.wr)
        while len(self.wr) < depth:
            self.wr.append([[] for _ in self])

        # We compute each WR as follows: for every vertex u
        # we compute a depth-bounded bfs in the graph induced
        # by all vertices _larger_ than u. For each vertex v we find
        # in this manner, we add u to v's WR-set at the respective depth.
        for iu in self:
            for d, layer in self.forward_bfs(iu, depth):
                if d <= old_depth:
                    continue
                for iv in layer:
                    self.wr[d-1][iv].append(iu)

        return self

    def brute_force_count(self, pattern):
        """
            Counts the number of times pattern appears as an ordered
            subgraph by brute force.
        """
        count = sum(1 for _ in self.brute_force_enumerate(pattern))
        return count

    def brute_force_enumerate(self, _pattern):
        from pattern import PatternMatch
        """
            Enumerate ordered subgraphs that are isomorphic to the
            given pattern.
        """
        for cand in itertools.combinations(self, len(_pattern)):
            mapping = list(zip(cand, _pattern))

            match = PatternMatch(self, _pattern)
            for iu, i in zip(cand, _pattern):
                match = match.extend(i, iu)
                if not match:
                    break
            else:
                assert match
                yield match            

    def match(self, iu, piece, partial_match=None, filtered_leaves=None, allowed_matches=None):
        """
            Returns all ordered sets X \\subseteq WR(iu)
            such that ordered graph induced by X \\cup \\{iu\\} 
            matches the provided piece.
        """
        from pattern import PatternMatch
        assert self.depth() >= len(piece.pattern)-1

        if filtered_leaves == None:
            filtered_leaves = set() 

        # Prepare basic match by adding the root
        missing_leaves = list(piece.leaves)

        log.debug("Wreach-sets of vertex %d:", iu)
        for d,wr in enumerate(self.wr):
            log.debug("(%d) %s", d, wr[iu])

        if partial_match:
            base_match = partial_match.extend(piece.root, iu)
            if base_match == None:
                return

            missing_leaves = [i for i in missing_leaves if not base_match.is_fixed(i)]
        else:
            base_match = PatternMatch(self, piece.pattern).extend(piece.root, iu)
            if base_match == None:
                return

        for match in self._match_rec(piece, base_match, missing_leaves, filtered_leaves, allowed_matches, 0):
            yield match

    def complete(self, piece, partial_match):
        """
            Completes the provided partial match to contain
            all vertices of the given piece. All roots
            must habe been set for this to work.
        """
        from pattern import PatternMatch        
        assert self.depth() >= len(piece.pattern)-1

        missing_leaves = [i for i in piece.leaves if not partial_match.is_fixed(i)]

        for match in self._match_rec(piece, partial_match, missing_leaves, set(), None, 0):
            yield match

    def _match_rec(self, piece, match, missing_leaves, filtered_leaves, allowed_matches, depth):
        assert(match[piece.root] != None)

        prefix = " "*(depth+2)
        log.debug("%sCALL %s %s %s", prefix, piece, match, missing_leaves)
        if len(missing_leaves) == 0:
            log.debug(prefix+"RETURN match %s", match)
            yield match
            return

        i = missing_leaves[-1]
        log.debug("%sMatching index %d for match %s", prefix, i, match)
        candidates = None
        if i in piece.nroots:
            # We need to pick this vertex from all wreach vertices.
            assert(piece.root_dist[i] != None)
            wr_dist = piece.root_dist[i]
            wr_alt = sorted(self.wreach_union(match[piece.root], 1, wr_dist))
            candidates = wr_alt

            log.debug("%sIndex is an nroot at wr-dist %s", prefix, wr_dist)
            log.debug("%s Restricted wreach set %s", prefix, wr_alt)
        else:
            # This vertex lies within the back-neighbourhood of
            # at least one already matched vertex. 
            anchors = match.find_anchors(i)
            assert(len(anchors) > 0)
            candidates = self.common_in_neighbours(anchors)

            log.debug("%sIndex has anchors %s", prefix, anchors)
            log.debug("%s Restricted wreach set %s", prefix, candidates)

        # Restrict to range as defined by already matched vertices.
        lower, upper = match.get_range(i)
        ilower = bisect.bisect_left(candidates, lower)
        iupper = bisect.bisect_right(candidates, upper)

        log.debug("%sRestricting candidates to range (%d,%d)", prefix, lower,upper)
        log.debug("%s Candidates %s", prefix, candidates)
        log.debug("%s Restricted to %s (indices (%d:%d))", prefix, candidates[ilower:iupper], ilower, iupper)

        # TODO: Optimize by testing whether allowe_leaves != None outside the loop
        for iv in candidates[ilower:iupper]:
            if i in filtered_leaves and i not in allowed_matches[iv]:
                log.debug("%sIndex %d is not allowed for vertex %d (allowed are %s)", prefix, i, iv, allowed_matches[iv])
                continue
            log.debug("%sPutting %d on index %d", prefix, iv, i)
            next_match = match.extend(i, iv)
            if next_match == None:
                continue
            for m in self._match_rec(piece, next_match, missing_leaves[:-1], filtered_leaves, allowed_matches, depth+1):
                yield m

    def _compare(self, mleaves, mroot, piece):
        """
            Tests whether the ordered graph induces by
            mleaves \\cup \\{mroot\\} matches the given piece.
        """
        mapping = list(zip(mleaves, piece.leaves))

        # Test whether edges between the root/rest are correct
        for mleaf, leaf in mapping:
            if piece.adjacent(leaf, piece.root) != self.adjacent(mleaf, mroot):
                return False

        # TODO: check remaining edges/non-edges
        for mappedA, mappedB in itertools.combinations(mapping, 2):
            mleafA, leafA = mappedA
            mleafB, leafB = mappedB
            if piece.adjacent(leafA, leafB) != self.adjacent(mleafA, mleafB):
                return False

        return True 

    def __repr__(self):
        return ','.join(map(lambda s: str(list(s)),self.wr[0]))
     

class TestGraphMethods(unittest.TestCase):

    def test_eq(self):
        from datastructures import Bimap
        G = Graph()
        G.add_edge('a','b')
        G.add_edge('b','c')
        G.add_edge('c','d')
        G.add_edge('a','d')

        self.assertEqual(G, G.copy())
        H = Graph()
        H.add_edges(G.edges())

        self.assertEqual(G, H)

        H.add_node('e')
        self.assertNotEqual(G, H)
        
        H.remove_node('e')
        self.assertEqual(G, H)
        
        H.add_edge('a', 'a')
        self.assertNotEqual(G, H)

        H.remove_loops()
        self.assertEqual(G, H)

        # Test __eq__ under relabeling
        mapping = Bimap()
        mapping.put('a', 1)
        mapping.put('b', 2)
        mapping.put('c', 3)
        mapping.put('d', 4)

        H = H.relabel(mapping)
        self.assertNotEqual(G, H)
        H = H.relabel(mapping)
        self.assertEqual(G, H)

    def test_hashing(self):
        edges = [(0,1),(1,2),(2,0),(3,1)]
        G = Graph()
        H = Graph()
        self.assertEqual(hash(G), hash(H))

        for u,v in edges:
            G.add_edge(u,v)

        for u,v in reversed(edges):
            H.add_edge(v,u)

        self.assertEqual(hash(G), hash(H))

        for u,v in edges:
            H.add_edge(u,v)

        self.assertEqual(hash(G), hash(H))

    def test_hashing_rnd(self):
        edges = []
        n = 100
        for _ in range(1000):
            edges.append((random.randint(0,n-1),random.randint(0,n-1)))

        G = Graph()
        H = Graph()
        for u,v in edges:
            G.add_edge(u,v)
            H.add_edge(v,u)
            self.assertEqual(hash(G), hash(H))

        for u,v in edges[:len(edges)//2]:
            G.remove_edge(u,v)
            H.remove_edge(v,u)
            self.assertEqual(hash(G), hash(H))

        sup_nodes = list(range(n+1, 2*n))
        H.add_nodes(sup_nodes)

        for _ in range(1000):
            H.add_edge(random.randint(n+1, 2*n-1), random.randint(n+1, 2*n-1))
        H = H.subgraph(range(n))

        self.assertEqual(hash(G), hash(H))        

        for i in range(n):
            G.remove_node(i)
            H.remove_node(i)
            self.assertEqual(hash(G), hash(H))


class TestLGraphMethods(unittest.TestCase):

    def test_ordering(self):
        G = Graph()
        G.add_edge('a','b')
        G.add_edge('b','c')
        G.add_edge('c','d')

        LG, mapping = G.to_lgraph()
        self.assertEqual(len(G), len(LG))

        for u,v in G.edges():
            iu, iv = mapping.index_of(u), mapping.index_of(v)
            self.assertTrue(LG.adjacent(iu,iv))


        for iu,iv in LG.edges():
            u, v = mapping.vertex_at(iu), mapping.vertex_at(iv)
            self.assertTrue(G.adjacent(u,v))

    def test_pattern_enum(self):
        G = Graph()
        G.add_edge('a', 'b')
        G.add_edge('b', 'c')

        patterns = [p for p,m in G.enum_patterns()]
        self.assertEqual(len(patterns), 3)

        num_pieces = [1,1,2]

        for num, pattern in zip(num_pieces, patterns):
            pieces = pattern.decompose()
            self.assertEqual(len(pieces), num)

class TestDiGraphMethods(unittest.TestCase):            
    
    def test_basics(self):
        G = DiGraph()
        G.add_arc(0,1)
        G.add_arc(0,2)
        G.add_arc(0,3)

        self.assertEqual(set([1,2,3]), G.out_neighbours(0))
        self.assertEqual(G.out_degree(0), 3)
        self.assertEqual(G.in_degree(0), 0)
        self.assertEqual(G.in_degree(1), 1)
        self.assertEqual(G.out_degree(1), 0)

if __name__ == '__main__':
    unittest.main()