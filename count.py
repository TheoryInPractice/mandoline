#!/usr/bin/env python3

from graph import Graph, load_graph
from pattern import PatternBuilder, Pattern

import sys

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




# Small path in longer path
G = load_graph('example-graphs/path.txt.gz')
H = PatternBuilder(4) \
        .add_edge(0,1).add_edge(0,2).add_edge(1,3) \
        .build() 

# Triangle with tail in karate
G = load_graph('example-graphs/karate.txt.gz')
H = PatternBuilder(4) \
        .add_edge(0,1).add_edge(1,2).add_edge(0,2).add_edge(0, 3) \
        .build()

LG,mapping = G.to_lgraph()
LG.compute_wr(len(H)-1)

pieces = list(H.decompose())
for i,piece in enumerate(pieces):
    print(i, piece)
    print("  Leaves:", piece.leaves)
    print("  Adhesion:", piece.adhesion)

DP = DPTable()
candidates = [defaultdict(list) for _ in pieces]

for iu in LG:
    wreach = sorted(LG.wreach_all(iu))

    prev_adhesion = []
    for i, piece in enumerate(pieces):
        for ileaves in LG.match(iu, piece):
            # The set ileaves together with iu induces
            # a graph isomorphic to the current piece
            mapping = dict(zip(piece.leaves, ileaves))

            # Determine adhesion sets under current mapping
            curr_adhesion_mapped = tuple([mapping[iv] for iv in piece.adhesion])

            candidates[i][curr_adhesion_mapped].append(iu)

print()
for i, piece in enumerate(pieces):
    print("Matches for", piece, ":")
    print(" ",candidates[i])

