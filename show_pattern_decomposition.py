#!/usr/bin/env python3

from graph import Graph, load_graph
from pattern import PatternBuilder, Pattern
import math, random
import cairo


# G = load_graph('lesmiserables.txt.gz')
# G = load_graph('karate.txt.gz')

# LG = G.to_lgraph()
# LG.compute_wr(3)

# for i in LG:
# 	print(	','.join([str(LG.wr[d][i]) for d in range(LG.depth())]))


if __name__ == '__main__':
    width_in_inches, height_in_inches = 10, 10
    width_in_points, height_in_points = \
        width_in_inches * 72, height_in_inches * 72
    width, height = width_in_points, height_in_points

    offx, offy = 50, 50

    surface = cairo.SVGSurface("test.svg", width_in_points, height_in_points)
    ctx = cairo.Context(surface)

    ctx.save()
    
    ctx.set_source_rgb(1,1,1)
    ctx.rectangle(0, 0, width, height)
    ctx.fill()

    # Draw stuff
    H = None

    # Test fixed graphs
    # This one has three leaves
    # H = PatternBuilder(4) \
    #         .add_edge(0,1).add_edge(2,3).add_edge(0,3) \
    #         .add_edge(0,4).add_edge(2,4) \
    #         .build() 

    # Tree where adhesion sets decrease in size
    # H = PatternBuilder(6) \
    #         .add_edge(0,1).add_edge(1,2) \
    #         .add_edge(2,3).add_edge(1,4).add_edge(0,5) \
    #         .build()

    # Triangle with tail
    # H = PatternBuilder(4) \
    #     .add_edge(0,1).add_edge(1,2).add_edge(0,2).add_edge(0, 3) \
    #     .build()

    # Triangle with two tails
    # H = PatternBuilder(5) \
    #     .add_edge(0,1).add_edge(1,2).add_edge(0,2).add_edge(0, 3).add_edge(1,4) \
    #     .build()


    if not H: # Random graph
        n = 6
        builder = PatternBuilder(n) 
        for _ in range(2*n):
            u = random.randint(0, n)
            v = random.randint(0, n)
            builder.add_edge(u,v)
        H = builder.build()

    ctx.save()
    ctx.translate(offx, offy)
    bbox = H.draw(ctx) 
    offy += bbox.height + 20
    ctx.restore()

    for piece in H.decompose():
        ctx.save()
        ctx.translate(offx, offy)
        bbox = piece.draw(ctx)
        offy += bbox.height + 20
        ctx.restore()

    ctx.restore()
    ctx.show_page()
    surface.finish()

