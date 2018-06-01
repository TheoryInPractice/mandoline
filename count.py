#!/usr/bin/env python3

from graph import Graph, load_graph
from pattern import PatternBuilder, Pattern

from collections import defaultdict
import bisect
import math, random
import cairo

class Increment:
    def __init__(self):
        self.events = []

    def record(self, time):
        assert(len(self.events) == 0 or self.events[-1] < time)
        self.events.append(time)

    def query(self, time):
        # Find the number of events that happened
        # before 'time'.
        index = bisect.bisect_right(time)
        return index 

class DPTable:
    def __init__(self):
        self.events = defaultdict(lambda: Increment())

    def record(self, roots, adhesion, pos):
        adhesion = tuple(adhesion)
        self.events[adhesion].record(pos)

    def query(self, roots, adhesion, start, end):
        pass




G = load_graph('example-graphs/path.txt.gz')

H = PatternBuilder(4) \
        .add_edge(0,1).add_edge(0,2).add_edge(1,3) \
        .build() 

LG = G.to_lgraph()
LG.compute_wr(len(H)-1)

pieces = H.decompose()

DP = DPTable()

for iu in LG:
    wreach = sorted(LG.wreach_all(iu))

    print(iu, wreach)
    prev_adhesion = []
    for i, piece in enumerate(pieces):
        print(i, piece)
        for ileaves in LG.match(iu, piece):
            # The set ileaves together with iu induces
            # a graph isomorphic to the current piece
            mapping = dict(zip(piece.leaves, ileaves))

            # Determine adhesion sets under current mapping:
            # the previous set in order to query the number of partial
            # subgraphs we can assemble, the current ahdesion set to 
            # store the result.
            prev_adhesion_mapped = [mapping[iv] for iv in prev_adhesion]
            curr_adhesion_mapped = [mapping[iv] for iv in piece.adhesion]

            count = 0
            if len(prev_adhesion_mapped) == 0:
                assert i == 0, "{}".format(i) # For now!
                count = 1 
                DP.record([piece.root], curr_adhesion_mapped, iu)
            else:
                start, end = piece.insertion_interval()
                stard, end = mapping[start], mapping[end]

                print("Search match in [{},{}] with adhesion {}".format(start, end, prev_adhesion_mapped))

                # count = DP.query(piece.previous_roots, prev_adhesion_mapped, start, end )

        prev_adhesion = list(piece.adhesion)
        print()