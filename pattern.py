from graph import Graph, LGraph 
import math
import unittest

class PatternMatch:
    def __init__(self, LG, pattern):
        self.LG = LG
        self.pattern = pattern
        self.vertex_to_index = dict()
        self.index_to_vertex = tuple([None] * len(pattern))

    def __len__(self):
        return len(self.index_to_vertex)

    def extend(self, u, i):
        if self.index_to_vertex[i] != None or u in self.vertex_to_index:
            return None # Invalid: index or vertex already assigned

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
            Restricts the current pattern to the provided
            indices, all others are 'forgotten'
        """
        res = PatternMatch(self.LG, self.pattern)
        for i in indices:
            iu = self.index_to_vertex[i]
            assert iu != None
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

    def get_range(self,i):
        """
            Returns the possible range of values for the index i,
            as dictated by the already set indices to the left and
            right of i. The returned range [a,b] is inclusive on both
            sides, so values a,b are allowed. If b < a, no possible values
            exist.
        """
        ileft, iright = i, i
        while ileft >= 0 and self.index_to_vertex[ileft] == None:
            ileft -= 1

        k = len(self)
        while iright < k and self.index_to_vertex[iright] == None:
            iright += 1
        
        # TODO: this can be slightly improved, since every slot between
        # i and ileft (iright) is empty, this narrows the range further.
        lbound = 0 if ileft < 0 else self.index_to_vertex[ileft]+1
        rbound = len(self.LG) if iright >= k else self.index_to_vertex[iright]-1

        return (lbound, rbound)


    def indices_to_vertices(self, indices):
        return [self.index_to_vertex[i] for i in sorted(indices)]

    def __hash__(self):
        return self.index_to_vertex.__hash__()

    def __eq__(self, other):    
        if isinstance(other, self.__class__):
            return self.LG == other.LG and self.pattern == other.pattern and self.index_to_vertex == other.index_to_vertex
        return False

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
        res.compute_wr(self.size)

        return res

class BoundingBox:
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


class Pattern(LGraph):
    def __init__(self):
        super().__init__()

    def __hash__(self):
        # fnv-style 64 bit hash
        fnv_prime = 1099511628211
        fnv_offset = 14695981039346656037
        modulo = 2 << 64 

        graph_hash = fnv_offset
        for r,wr in enumerate(self.wr):
            layer_hash = fnv_offset
            for iu, iN in enumerate(self.wr[0]):
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

        if len(self.wr) != len(other.wr):
            return False
        for layerSelf, layerOther in zip(self.wr, other.wr):
            if layerSelf != layerOther:
                return False
        return True

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
            wreach = self.wreach_all(iu)
            seen |= set(wreach)
            roots.insert(0, iu)

        pieces = []
        prev_leaves = []
        for i,iu in enumerate(roots):
            wreach = self.wreach_all(iu)
            pieces.append(Piece(self, iu, roots[:i], wreach))
            prev_leaves = wreach
        return pieces

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
            for iv in self.in_neighbours(iu):
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
            for iv in self.in_neighbours(iu):
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

        # At this point this is pure paranoia
        self.leaves.sort() 

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
                ", ".join(map(str,self.edges())))

    def adjacent(self, u, v):
        assert u == self.root or u in self.leaves
        assert v == self.root or v in self.leaves
        return self.pattern.adjacent(u,v)

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

        colors = {}
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