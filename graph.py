import itertools
from collections import defaultdict as defaultdict
import gzip
from operator import itemgetter
from itertools import chain, combinations, permutations

import bisect
import math

import unittest

class indexmap:
    def __init__(self, size):
        self.vertex_to_index = {}
        self.index_to_vertex = [None] * size

    def __len__(self):
        return len(self.index_to_vertex)

    def put(self, index, vertex):
        if index >= len(self) or index < 0:
            raise IndexError()
        self.vertex_to_index[vertex] = index
        self.index_to_vertex[index] = vertex

    def vertex_at(self, index):
        if index >= len(self) or index < 0 or self.index_to_vertex[index] == None:
            raise IndexError()
        return self.index_to_vertex[index]

    def index_of(self, vertex):
        return self.vertex_to_index[vertex]

    def indices_of(self, vertices):
        return [self.vertex_to_index[v] if v in self.vertex_to_index else None for v in vertices]

    def __repr__(self):
        return 'IM['+','.join(self.index_to_vertex)+']'

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

    def to_lgraph(self, order=None):
        """
            Returns an ordered graph with node indices from 0 to n-1 and 
            an indexmap to recover the original ids.
            The nodes are sorted by (degree, label) in descending order
            in order to keep the wcol number small.
        """
        lgraph = LGraph()

        if order == None:
            order = [(self.degree(u), u) for u in self.nodes] 
            order.sort(reverse=True)
            order = [u for _,u in order]

        mapping = indexmap(len(order))
        for iu,u  in enumerate(order):
            iN = mapping.indices_of(self.neighbours(u))
            iN = [x for x in iN if x != None] # Remove 'None' values
            lgraph._add_node(iu, iN)
            mapping.put(iu, u)

        return lgraph, mapping

    def enum_patterns(self):
        found = set()
        for order in permutations(sorted(self)):
            pat, indexmap = self.to_pattern(order)
            if pat in found:
                continue
            found.add(pat)
            yield pat, indexmap

    def to_pattern(self, order):
        from pattern import Pattern
        """
            Returns a `pattern' (and ordered graph with node indices from 0 to n-1)
            and an indexmap to recover the original ids. 
        """
        pat = LGraph()

        mapping = indexmap(len(order))
        for iu,u  in enumerate(order):
            iN = mapping.indices_of(self.neighbours(u))
            iN = [x for x in iN if x != None] # Remove 'None' values
            pat._add_node(iu, iN)
            mapping.put(iu, u)
        pat.compute_wr(len(pat))
        return Pattern(pat), mapping

class LGraph:
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
            i = bisect.bisect_left(self.wr[0][iu], iv)
            return i != len(self.wr[0][iu]) and self.wr[0][iu][i] == iv
            # return iv in self.wr[0][iu]
        else:
            i = bisect.bisect_left(self.wr[0][iv], iu)
            return i != len(self.wr[0][iv]) and self.wr[0][iv][i] == iu
            # return iu in self.wr[0][iv]

    def edges(self):
        for u in self:
            for v in self.wr[0][u]:
                yield (v,u)

    def wreach_iter(self):
        for iu in self:
            yield iu, sorted(self.wreach_all(iu))

    def wreach_all(self, iu):
        res = []
        for wr in self.wr:
            res.extend(wr[iu])
        return res

    def _add_node(self, iu, iN):
        iN = list(sorted(iN))
        self.wr[0].append(iN)
        self.outN.append([])

        for iv in iN:
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
                match = match.extend(iu,i)
                if not match:
                    break
            else:
                assert match
                yield match            

            # for (mu,u), (mv,v) in itertools.combinations(mapping, 2):
            #     if _pattern.adjacent(u,v) != self.adjacent(mu, mv):
            #         break
            # else:
            #     yield cand

    def match(self, iu, piece, partial_match=None):
        from pattern import PatternMatch
        """
            Returns all ordered sets X \\subseteq WR(iu)
            such that ordered graph induced by X \\cup \\{iu\\} 
            matches the provided piece.
        """
        assert self.depth() >= piece.depth()
        # TODO: do we really have to go to maximal depth wreach?
        # It seems like we could be more discerning here.
        wreach = sorted(self.wreach_all(iu))

        # Prepare basic match by adding the root
        missing_leaves = list(piece.leaves)

        if partial_match:
            base_match = partial_match.extend(iu, piece.root)
            if base_match == None:
                return

            missing_leaves = [i for i in missing_leaves if not base_match.is_fixed(i)]

        else:
            base_match = PatternMatch(self, piece.pattern).extend(iu, piece.root)
            if base_match == None:
                return

        # Early out: match already complete
        if len(missing_leaves) == 0:
            yield base_match
            return

        for match in self._match_rec(piece, wreach, base_match, missing_leaves):
            yield match

    def _match_rec(self, piece, wreach, match, missing_leaves):
        assert(match.is_fixed(piece.root))

        if len(missing_leaves) == 0:
            yield match
            return

        i = missing_leaves[-1]
        candidates = None
        if i in piece.nroots:
            # We need to pick this vertex from all wreach vertices.
            candidates = wreach
        else:
            # This vertex lies within the back-neighbourhood of
            # at least one already matched vertex. 
            anchors = match.find_anchors(i)
            assert(len(anchors) > 0)
            candidates = self.common_in_neighbours(anchors)

        # Restrict to range as defined by already matched vertices.
        lower, upper = match.get_range(i)
        ilower = bisect.bisect_left(candidates, lower)
        iupper = bisect.bisect_right(candidates, upper)
        for iv in candidates[ilower:iupper]:
            next_match = match.extend(iv, i)
            if next_match == None:
                continue
            for m in self._match_rec(piece, wreach, next_match, missing_leaves[:-1]):
                yield m


    def _match_rec_old(self, i, leaf_index, match, wreach, wreach_indexed, missing_leaves, debug):
        if debug:
            print("    Wreach: ", wreach)
            print("    Current match is", match)
            print("    Extension candidates for index {}:".format(missing_leaves[leaf_index]), wreach)
            print("      -- trimmed by index:", wreach[i:])

        # Narrow down search range using the bounds dictated 
        # by the current match.    
        lower, upper = match.get_range(missing_leaves[leaf_index])

        ilower = bisect.bisect_left(wreach, lower)
        ilower = max(ilower, i)
        iupper = bisect.bisect_right(wreach, upper)
        candidates = wreach_indexed[ilower:iupper]
       
        if debug:
            print("    Trimmed to range {},{}:".format(lower,upper), wreach[ilower:iupper], "(", candidates ,")")

        if leaf_index == len(missing_leaves)-1:
            # Every match here is a complete match and we return it
            for j,iv in candidates:
                if debug:
                    print("      Trying", iv)
                next_match = match.extend(iv, missing_leaves[leaf_index])
                if next_match:
                    yield next_match
        else:
            # Recurse on partial matches
            for j,iv in candidates:
                next_match = match.extend(iv, missing_leaves[leaf_index])
                if next_match:
                    for m in self._match_rec(j+1, leaf_index+1, next_match, wreach, wreach_indexed, missing_leaves, debug):
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
        return ','.join(map(str,self.wr[0]))
     


class TestLGraphMethods(unittest.TestCase):

    def test_ordering(self):
        G = Graph()
        G.add_edge('a','b')
        G.add_edge('b','c')
        G.add_edge('c','d')

        LG, mapping = G.to_lgraph()
        self.assertEqual(len(G), len(LG))

        print(LG)
        print(mapping)

        for u,v in G.edges():
            iu, iv = mapping.index_of(u), mapping.index_of(v)
            print((u,v), 'mapped to',(iu,iv))
            self.assertTrue(LG.adjacent(iu,iv))


        for iu,iv in LG.edges():
            u, v = mapping.vertex_at(iu), mapping.vertex_at(iv)
            print((iu,iv), 'mapped to',(u,v))
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

if __name__ == '__main__':
    unittest.main()