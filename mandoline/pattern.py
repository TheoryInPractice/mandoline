from .graph import Graph
import math
import unittest

from sortedcontainers import SortedSet

import logging
from .helpers import short_str

log = logging.getLogger(__name__)

class BoundingBox:
    """
        Helper class for drawing only.
    """
    def __init__(self):
        self.left = None
        self.right = None
        self.top = None
        self.bottom = None

    @property
    def width(self):
        return self.right-self.left

    @property
    def height(self):
        return self.bottom-self.top

    def include(self, x, y):
        self.left = x if self.left == None else min(self.left, x)
        self.right = x if self.right == None else max(self.right, x)
        self.top = y if self.top == None else min(self.top, y)
        self.bottom = y if self.bottom == None else max(self.bottom, y)

    def move(self, dx, dy):
        self.left += dx
        self.right += dx
        self.top += dy
        self.bottom += dy

    def __repr__(self):
        return "BBox[{},{},{},{}]".format(self.top, self.right, self.bottom, self.left)


class PatternMatch:
    def __init__(self, LG, pattern):
        self.LG = LG
        self.pattern = pattern
        self.vertex_to_index = {}
        self.index_to_vertex = tuple([None] * len(pattern))

        self._range_cache = {}

    def __len__(self):
        return len(self.index_to_vertex)

    def __getitem__(self, i):
        return self.index_to_vertex[i]

    def is_complete(self):
        """ Returns true if all indices have been filled with vertices """
        return len(self.vertex_to_index) == len(self.index_to_vertex)

    def covers_piece(self, piece):
        """ Returns whether all indices of the provided piece are filled in this match """
        if self.index_to_vertex[piece.root] == None:
            return False

        for i in piece.leaves:
            if self.index_to_vertex[i] == None:
                return False

        return True

    def matched_vertices(self):
        indices = range(len(self.index_to_vertex))
        return [(i,self.index_to_vertex[i]) for i in indices if self.index_to_vertex[i] != None ]

    def rightmost_unmatched(self):
        # TODO unused, deprecate?
        for i,iv in reversed(enumerate(self.index_to_vertex)):
            if iv == None:
                return i
        return None

    def find_anchors(self, i):
        """
            For a given index i, returns vertices of _matched_ indices to
            the right of i that are neighbours in the underlying pattern.
        """
        res = [self.index_to_vertex[j] for j in self.pattern.out_neighbours[i]]
        res = [iv for iv in res if iv != None ]
        return res

    def extend_multiple(self, extensions):
        """
            Performs multiple extensions, provided as an array of the form
            [(index, vertex) ...]. Returns None if the matches are incompatible
            with the pattern or order.
        """
        res = self

        for i, u in extensions:
            res = res.extend(i, u)
            if res == None:
                return None

        return res

    def extend(self, i, u):
        """
            Extends the match by putting u on the index i. If the extension
            is invalid (u has the wrong edges towards already matched vertices
            or u break the order of the match) this function returns None.
        """
        if self.index_to_vertex[i] != None or u in self.vertex_to_index:
            return None # Invalid: index or vertex already assigned

        left, right = self.get_range(i)
        if u < left or u > right:
            log.debug("Warning: tried extending %d with %d at %d, invalid range.", self, u, i)
            log.debug("   (Range: %d, %d, %d)",  left, u, right )
            return None

        # Test whether extension is valid
        for v,j in self.vertex_to_index.items():
            if self.pattern.adjacent(i,j) != self.LG.adjacent(u,v):
                return None

        # Perform extension
        res = PatternMatch(self.LG, self.pattern)
        res.vertex_to_index = dict(self.vertex_to_index)
        res.vertex_to_index[u] = i

        temp = list(self.index_to_vertex)
        temp[i] = u
        res.index_to_vertex = tuple(temp)

        return res

    def restrict_to(self, indices):
        """
            Restricts the current match to the provided
            indices, all others are 'forgotten'
        """
        res = PatternMatch(self.LG, self.pattern)
        for i in indices:
            iu = self.index_to_vertex[i]
            assert iu != None, "Cannot restrict {} to {}, index {} is not set.".format(self, indices, i)
            res.vertex_to_index[iu] = i

        temp = list(self.index_to_vertex)
        indices_set = set(indices)
        for i in range(len(temp)):
            temp[i] = temp[i] if i in indices_set else None
        res.index_to_vertex = tuple(temp)

        return res

    def contains(self, u):
        return u in self.vertex_to_index

    def is_fixed(self, i):
        return self.index_to_vertex[i] != None

    def get_range(self, i):
        """
            Returns the possible range of values for the index i,
            as dictated by the already set indices to the left and
            right of i. The returned range [a,b] is inclusive on both
            sides, so values a,b are allowed. If b < a, no possible values
            exist.
        """
        if i in self._range_cache:
            return self._range_cache[i]

        assert(self.index_to_vertex[i] == None)
        # v = self.index_to_vertex[i]
        # if v != None:
        #     return (v,v)

        ileft = iright = i
        while ileft >= 0 and self.index_to_vertex[ileft] == None:
            ileft -= 1

        k = len(self)
        while iright < k and self.index_to_vertex[iright] == None:
            iright += 1

        lbound = i
        if ileft >= 0:
            lbound = self.index_to_vertex[ileft]+(i-ileft)

        rbound = len(self.LG)-(k-i)
        if iright < k:
            rbound = self.index_to_vertex[iright]-(iright-i)

        # Cache result
        self._range_cache[i] = (lbound, rbound)
        return (lbound, rbound)

    def indices_to_vertices(self, indices):
        return [self.index_to_vertex[i] for i in sorted(indices)]

    def get_adhesion(self, size):
        """
            Returns the 'size' first vertices matched here as a tuple.
        """
        return tuple([self.index_to_vertex[i] for i in range(size)])

    def __hash__(self):
        return self.index_to_vertex.__hash__()

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.LG == other.LG and self.pattern == other.pattern and self.index_to_vertex == other.index_to_vertex
        # return False

    def __repr__(self):
        return "PM "+ ",".join(map(lambda x: "_" if x == None else str(x), self.index_to_vertex))

class PatternBuilder:
    def __init__(self, size):
        self.graph = Graph()
        for i in range(size):
            self.graph.add_node(i)
        self.size = size

    def add_edge(self, u, v):
        self.graph.add_edge(u,v)
        return self

    def build(self):
        # LG,_ = self.graph.to_lgraph(range(self.size))
        # res = Pattern.from_lgraph(LG)
        res,_ = self.graph.to_pattern(range(self.size))

        return res

class Pattern:
    def __init__(self, LG):
        n = len(LG)
        LG.compute_wr(n)

        self.wreach = [None] * n
        self.in_neighbours = [None] * n
        self.wr_dist = [None] * n
        self.out_neighbours = [SortedSet() for _ in range(n)]
        for iu in LG:
            self.wreach[iu] = LG.wreach_all(iu)
            self.in_neighbours[iu] = SortedSet(LG.in_neighbours(iu))
            for iv in self.in_neighbours[iu]:
                self.out_neighbours[iv].add(iu)

        # Compute pairwise 'wr-distance'. Each piece will
        # inherit a row of this matrix (technically, we only need
        # one row for every root but I don't want to compute them here).
        # Note that a a distance 0 here indicates that a vertex is
        # in the in-neigbourhood of another vertex, a distance d indicates
        # that we need to compute wr[d] in order to see the connection.
        for root in range(n):
            self.wr_dist[root] = [None] * n
            for d in range(n):
                for iu in LG.wr[d][root]:
                    self.wr_dist[root][iu] = d

        self.cached_hash = self._compute_hash()

    @staticmethod
    def from_disconnected(LG, max_wr_dist):
        res = Pattern(LG)
        for root,row in enumerate(res.wr_dist):
            # Every node past 'root' (including 'root') will have
            # 'None' by definition. If we find 'None' _before_ 'root',
            # we change it to max_wr_dist instead since we are promised
            # that the underlying pattern is connected.
            for i in range(root):
                if row[i] == None:
                    assert i < root
                    row[i] = max_wr_dist

        # In case the hash uses the wr_dist at some point
        res.cached_hash = res._compute_hash()
        return res

    def __len__(self):
        return len(self.wreach)

    def __iter__(self):
        return iter(range(len(self.wreach)))

    def edges(self):
        for u in self:
            for v in self.in_neighbours[u]:
                yield (v,u)

    def __hash__(self):
        return self.cached_hash

    def monotone_bfs(self, iu):
        seen = set([iu])
        layers = [SortedSet([iu])]
        while True:
            next_layer = SortedSet()
            for iv in layers[-1]:
                next_layer.update(set(self.in_neighbours[iv]) - seen)
            if len(next_layer) == 0:
                break
            seen.update(next_layer)
            layers.append(next_layer)
        return layers

    def monotone_reachable(self, iu):
        res = set()
        frontier = [iu]
        while len(frontier) > 0:
            iv = frontier.pop()
            if iv in res:
                continue

            res.add(iv)
            frontier += self.in_neighbours[iv]
        return res

    def _compute_hash(self):
        # fnv-style 64 bit hash
        fnv_prime = 1099511628211
        fnv_offset = 14695981039346656037
        modulo = 1 << 64 # 2**64

        graph_hash = fnv_offset
        layer_hash = fnv_offset
        for iu, iN in enumerate(self.in_neighbours):
            vertex_hash = fnv_offset
            for iv in iN:
                vertex_hash = vertex_hash ^ iu
                vertex_hash = (vertex_hash * fnv_prime) % modulo
                vertex_hash = vertex_hash ^ iv
                vertex_hash = (vertex_hash * fnv_prime) % modulo
            layer_hash = layer_hash ^ vertex_hash
            layer_hash = (layer_hash * fnv_prime) % modulo
        graph_hash = graph_hash ^ layer_hash
        graph_hash = (graph_hash * fnv_prime) % modulo
        return graph_hash

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False

        if len(self.in_neighbours) != len(other.in_neighbours):
            return False
        if len(self.wreach) != len(other.wreach):
            return False
        if self.in_neighbours != other.in_neighbours:
            return False
        if self.wreach != other.wreach:
            return False

        return True

    def adjacent(self, iu, iv):
        if iv < iu:
            return iv in self.in_neighbours[iu]
        else:
            return iu in self.in_neighbours[iv]

    def adjacent_ordered(self, iu, iv):
        # Assumes iu < iv
        return iu in self.in_neighbours[iv]

    def decompose(self):
        """
            Decompose ordered pattern into pieces, e.g.
            a set of (independent) vertices whose joint
            WReach-sets cover all of the graph.
        """

        # First compute roots
        seen = set()
        roots = []

        for iu in reversed(range(len(self))):
            if iu in seen:
                continue # Already covered
            seen |= set(self.wreach[iu])
            roots.insert(0, iu)

        pieces = []

        for i,iu in enumerate(roots):
            wreach = self.wreach[iu]
            pieces.append(Piece(self, iu, roots[:i], wreach))

        return pieces

    def to_single_piece(self):
        """
            Converts this pattern into a single piece. We need this
            for the counting part, where we encounter disconnected Patterns
            for which we know that they are nonetheless inside a single wreach
            set.
        """
        root = len(self)-1
        leaves = list(range(root))
        return Piece(self, root, [], leaves)

    def __repr__(self):
        return 'PT'+','.join(map(lambda s: short_str(list(s)),self.in_neighbours))

    def draw_subgraph(self, ctx, nodes, colors):
        node_col = (0,0,0)
        node_col_inactive = (0.75,0.75,0.75)
        edge_col = (0,0,0)
        edge_col_inactive = (0.75,0.75,0.75)

        rad = 5 # Node radius
        offy = 0
        offx = 0
        diffx = 30

        def coords(i):
            return offx + i*diffx, offy

        # Compute bounding box
        bbox = BoundingBox()
        bbox.include(0,0)
        for iu in self:
            ux, uy = coords(iu)
            bbox.include(ux, uy)
            for iv in self.in_neighbours[iu]:
                vx, vy = coords(iv)
                if iv % 2:
                    bbox.include((ux+vx)/2,uy+(ux-vx)/2)
                else:
                    bbox.include((ux+vx)/2,uy-(ux-vx)/2)

        offy = -bbox.top
        offx = -bbox.left

        # Draw edges
        ctx.set_line_width(1.5)
        for iu in self:
            ux, uy = coords(iu)
            for iv in self.in_neighbours[iu]:
                vx, vy = coords(iv)

                if iu not in nodes or iv not in nodes:
                    ctx.set_source_rgb(*edge_col_inactive)
                else:
                    ctx.set_source_rgb(*edge_col)

                if iv == iu-1:
                    ctx.move_to(ux, uy)
                    ctx.line_to(vx, vy)
                elif iv % 2:
                    ctx.arc((ux+vx)/2,uy,(ux-vx)/2,0,math.pi)
                else:
                    ctx.arc((ux+vx)/2,uy,(ux-vx)/2,math.pi,0)
                ctx.stroke()

        # Draw nodes
        for iu in self:
            ux, uy = coords(iu)
            ctx.arc(ux, uy, rad, 0, 2*math.pi)
            if iu in colors:
                ctx.set_source_rgb(*colors[iu])
            elif iu not in nodes:
                ctx.set_source_rgb(*node_col_inactive)
            else:
                ctx.set_source_rgb(*node_col)
            ctx.fill()

        # Return bounding box
        bbox.move(offx, offy)
        return bbox

    def draw(self, ctx):
        return self.draw_subgraph(ctx, set(self), set())

class Piece:
    def __init__(self, pattern, root, previous_roots, leaves):
        super().__init__()
        self.root = root
        assert len(previous_roots) == 0 or previous_roots[-1] != self.root

        self.pattern = pattern
        self.leaves = list(leaves)

        # wr-distance from every leaf to the root
        self.root_dist = pattern.wr_dist[root]

        # At this point this is pure paranoia
        self.leaves.sort()

        # Compute 'covering' set of vertices from which _all_
        # leaves can be reached monotonically. The idea is that
        # these vertices are the only ones that we have to find
        # in WReach, all other vertices can then be found via
        # bfs through back_neighbour sets (which are much smaller)
        nroots = []
        remainder = SortedSet(leaves)
        remainder.add(root)

        while len(remainder) > 0:
            nr = remainder[-1]
            del remainder[-1]
            nroots.append(nr)

            covered = pattern.monotone_reachable(nr)
            remainder -= covered
        self.nroots = SortedSet(nroots)

        # Pre-compute interval in which the very previous root
        # needs to be placed.
        self.interval = None
        if len(previous_roots) > 0:
            left = right = previous_roots[-1]
            while left not in self.leaves and left > 0:
                left -= 1
            while right not in self.leaves and right != self.root:
                right += 1
            self.interval = (left, right)

        # This objects are (conceptually) immutable, so we can afford
        # to compute a slightly expensive hash.
        self.cached_hash = self._compute_hash()

    def __hash__(self):
        return self.cached_hash

    def _compute_hash(self):
        # fnv-style 64 bit hash
        fnv_prime = 1099511628211
        fnv_offset = 14695981039346656037
        modulo = 2 << 64

        h = fnv_offset
        h = ((h ^ self.root) * fnv_prime ) % modulo
        for u in self.leaves:
            h = ((h ^ u) * fnv_prime ) % modulo
        for u,v in sorted(self.edges()):
            h = ((h ^ u) * fnv_prime ) % modulo
            h = ((h ^ v) * fnv_prime ) % modulo

        return h

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False
        if self.root != other.root or self.leaves != other.leaves:
            return False

        if self.pattern == other.pattern:
            return True

        return set(self.edges()) == set(other.edges())

    def __iter__(self):
        return iter(self.leaves)

    def __repr__(self):
        return "Piece {} ({}) {{{}}}".format(
                " ".join(map(str,self.leaves)),
                 self.root,
                ", ".join(map(str,self.root_dist)))

    def adjacent(self, u, v):
        assert u == self.root or u in self.leaves
        assert v == self.root or v in self.leaves
        return self.pattern.adjacent(u,v)

    def root_equivalent(self, other):
        """
            Tests whether two piece are isomorphic up to the
            index of the root. Two conditions have to be met:
            the leave sets are the same and the roots have the
            same neighbourhood.
        """
        if self.pattern != other.pattern:
            return False
        pattern = self.pattern
        if self.leaves != other.leaves:
            return False

        assert(len(pattern.out_neighbours[self.root]) == 0)
        assert(len(pattern.out_neighbours[other.root]) == 0)
        return pattern.in_neighbours[self.root] == pattern.in_neighbours[other.root]

    def insertion_interval(self):
        return self.interval

    def edges(self):
        for u,v in self.pattern.edges():
            if u not in self.leaves and u != self.root:
                continue
            if v not in self.leaves and v != self.root:
                continue
            yield (u,v)

    def num_leaves(self):
        return len()

    def depth(self):
        return len(self.pattern)

    def draw(self, ctx):
        rootset = set([self.root])
        active = set(self.leaves)
        root_col = (216/255.0,20/255.0,149/255.0)
        nroot_col = (77/255.0, 166/255.0, 255/255.0)

        colors = {}
        for ir in self.nroots:
            colors[ir] = nroot_col
        colors[self.root] = root_col

        return self.pattern.draw_subgraph(ctx, active | rootset, colors)

class TestPatternMethods(unittest.TestCase):

    def test_decompose(self):
        H = PatternBuilder(4) \
                .add_edge(0,1).add_edge(0,2).add_edge(1,3) \
                .build()
        self.assertEqual(len(list(H.decompose())), 2)

    def test_equality(self):
        H1 = PatternBuilder(4) \
                .add_edge(0,1).add_edge(0,2).add_edge(1,3) \
                .build()
        H2 = PatternBuilder(4) \
                .add_edge(1,3).add_edge(2,0).add_edge(1,0) \
                .build()
        H3 = PatternBuilder(4) \
                .add_edge(0,1).add_edge(1,2).add_edge(1,3) \
                .build()

        self.assertEqual(H1, H1)
        self.assertEqual(H2, H2)
        self.assertEqual(H3, H3)

        self.assertEqual(H1.__hash__(), H2.__hash__())
        self.assertEqual(H1, H2)

        self.assertNotEqual(H1.__hash__(), H3.__hash__())
        self.assertNotEqual(H1, H3)

if __name__ == '__main__':
    unittest.main()
