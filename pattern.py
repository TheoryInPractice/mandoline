from graph import Graph, LGraph 
import math
import unittest

class PatternBuilder:
    def __init__(self, size):
        self.graph = Graph()
        for i in range(size):
            self.graph.add_node(i)            

    def add_edge(self, u, v):
        self.graph.add_edge(u,v)
        return self

    def build(self):
        res = Pattern()

        for u in range(self.graph.get_max_id()+1):
            res._add_node(u, self.graph.neighbours(u))

        res.compute_wr(len(self.graph))

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
            adhesion = []
            if i > 0:
                adhesion = set(wreach) & set(prev_leaves)
                adhesion = sorted(adhesion)
            pieces.append(Piece(self, iu, roots[:i], wreach, adhesion))
            prev_leaves = wreach
        return pieces

class Piece:
    def __init__(self, pattern, root, previous_roots, leaves, adhesion):
        super().__init__()
        self.root = root
        assert len(previous_roots) == 0 or previous_roots[-1] != self.root

        self.pattern = pattern
        self.leaves = list(leaves)
        self.adhesion = list(adhesion)
        self.previous_roots = list(previous_roots)

        # At this point this is pure paranoia
        self.leaves.sort() 
        self.adhesion.sort() 
        self.previous_roots.sort()

        # Pre-compute interval in which the very previous root
        # needs to be placed.
        self.interval = None
        if len(self.previous_roots) > 0:
            left = right = self.previous_roots[-1]
            while left not in self.leaves and left > 0:
                left -= 1
            while right not in self.leaves and right != self.root:
                right += 1
            self.interval = (left, right)


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
        prev_root_col = (100/255.0, 100/255.0, 100/255.0 )
        adhesion_col = (9/255.0, 167/255.0, 216/255.0)

        colors = {}
        colors[self.root] = root_col
        for iu in self.adhesion:
            colors[iu] = adhesion_col

        for ir in self.previous_roots:
            colors[ir] = prev_root_col

        return self.pattern.draw_subgraph(ctx, active | rootset, colors)

class TestPatternMethods(unittest.TestCase):

    def test_decompose(self):
        H = PatternBuilder(4) \
                .add_edge(0,1).add_edge(0,2).add_edge(1,3) \
                .build() 
        self.assertEqual(len(list(H.decompose())), 2)


if __name__ == '__main__':
    unittest.main()