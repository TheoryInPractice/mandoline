import itertools
from collections import defaultdict as defaultdict
import gzip
from operator import itemgetter
from itertools import chain, combinations

import bisect

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

    def to_lgraph(self, order=None):
        """
            Returns graph with node indices from 0 to n-1 and 
            a mapping dictionary to recover the original ids.
            The nodes are sorted by (degree, label) in descending
            order.
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
            Returns all ordered sets X \subseteq WR(iu)
            such that X \cup \{iu\} matches (in edges/non-edges
            in the provided order) the provided piece.
        """
        assert self.depth() >= piece.depth()
        wreach = sorted(self.wreach_all(iu))
        wreach_indexed = list(zip(range(len(wreach)), wreach))

        # Initialize stack with basic match (mapping root of piece to iu),
        # or extend the provided partial match. In the latter case, we also
        # need to determine which leaves still have to be fixed (missing_leaves)
        debug = False
        leaf_index = 0 
        missing_leaves = list(piece.leaves)

        if partial_match:
            base_match = partial_match.extend(iu, piece.root)
            if base_match == None:
                return

            missing_leaves = [i for i in missing_leaves if not base_match.is_fixed(i)]

            debug = True

            if debug:
                print("    base match", base_match)
                print("    missing leaves", missing_leaves)
        else:
            base_match = PatternMatch(self, piece.pattern).extend(iu, piece.root)
            if base_match == None:
                return

        if len(missing_leaves) == 0:
            yield base_match
            return

        stack = [(0, leaf_index, base_match)]

        if debug:
            print("    Wreach: ", wreach)

        while len(stack) > 0:
            if debug:
                print("   ", stack)
            i, leaf_index, match = stack.pop()

            if debug:
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
                # A match here means that we put the current match on
                # the stack and work on the new match instead
                for j,iv in candidates:
                    next_match = match.extend(iv, missing_leaves[leaf_index])
                    if next_match:
                        # Save were to continue with previous match
                        stack.append((j+1,leaf_index,match))
                        # Save current match on stack
                        stack.append((j+1,leaf_index+1,next_match))
                        break

        # for cand in itertools.combinations(wreach, len(piece.leaves)):
        #     matched = self._compare(cand, iu, piece)
        #     if matched:
        #         yield cand

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

if __name__ == '__main__':
    unittest.main()