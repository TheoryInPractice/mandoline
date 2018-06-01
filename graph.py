import itertools
from collections import defaultdict as defaultdict
import gzip
from operator import itemgetter
from itertools import chain, combinations

import unittest

def load_graph(filename):
    res = Graph()
    with gzip.open(filename, 'r') as filebuf:
        for l in filebuf:
            u, v = l.decode().split()
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

    def __contains__(self,u):
        return u in self.nodes

    def __iter__(self):
        return iter(self.nodes)

    def __len__(self):
        return len(self.nodes)
    
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
        self.nodes.add(u)

    def add_edge(self,u,v):
        self.nodes.add(u)
        self.nodes.add(v)        
        self.adj[u].add(v)
        self.adj[v].add(u)

    def remove_edge(self,u,v):
        self.adj[u].discard(v)
        self.adj[v].discard(u)

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

    def to_lgraph(self):
        """
            Returns graph with node indices from 0 to n-1 and 
            a mapping dictionary to recover the original ids.
            The nodes are sorted by (degree, label) in descending
            order.
        """
        lgraph = LGraph()

        order = [(self.degree(u), u) for u in self.nodes] 
        order.sort(reverse=True)
        for _,u in order:
            lgraph._add_node(u, self.neighbours(u))

        return lgraph

class LGraph:
    def __init__(self):
        self.positions = {}
        self.labels = {}
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

    def adjacent(self, iu, iv):
        # TODO: Use that wr-sets are ordered 
        # for quicker access!
        if iv < iu:
            return iv in self.wr[0][iu]
        else:
            return iu in self.wr[0][iv]

    def edges(self):
        for u in self:
            for v in self.wr[0][u]:
                yield (v,u)

    def wreach_all(self, iu):
        res = []
        for wr in self.wr:
            res.extend(wr[iu])
        return res

    def label(self, iu):
        return self.labels[iu]

    def position(self, u):
        return self.positions[u]

    def _add_node(self, u, N):
        iu = self.positions[u] = len(self) 
        self.labels[iu] = u

        inN = []
        for v in N:
            if v not in self.positions:
                continue
            inN.append(self.positions[v])

        inN.sort()
        self.wr[0].append(inN)
        self.outN.append([])

        for iv in inN:
            self.outN[iv].append(iu)

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
        while len(self.wr) <= depth:
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
                    self.wr[d][iv].append(iu)

        return self

    def match(self, iu, piece):
        """
            Returns all ordered sets X \subseteq WR(iu)
            such that X \cup \{iu\} matches (in edges/non-edges
            in the provided order) the provided piece.
        """
        assert self.depth() >= piece.depth()
        wreach = sorted(self.wreach_all(iu))

        for cand in itertools.combinations(wreach, len(piece.leaves)):
            matched = self._compare(cand, iu, piece)
            if matched:
                yield cand

    def _compare(self, mleaves, mroot, piece):
        """
            Tests whether the ordered graph induces by
            mleaves \cup \{mroot\} matches the given piece.
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
        return ','.join(map(str,zip(self.wr[0], self.outN)))


class TestLGraphMethods(unittest.TestCase):

    def test_ordering(self):
        G = Graph()
        G.add_edge('a','b')
        G.add_edge('b','c')
        G.add_edge('c','d')

        LG = G.to_lgraph()
        self.assertEqual(len(G), len(LG))

        print(LG.labels)
        print(LG)

        for u,v in G.edges():
            iu, iv = LG.position(u), LG.position(v)
            print((u,v), 'mapped to',(iu,iv))
            self.assertTrue(LG.adjacent(iu,iv))


        for iu,iv in LG.edges():
            u, v = LG.label(iu), LG.label(iv)
            print((iu,iv), 'mapped to',(u,v))
            self.assertTrue(G.adjacent(u,v))

if __name__ == '__main__':
    unittest.main()