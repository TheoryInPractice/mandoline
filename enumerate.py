#!/usr/bin/env python3

from graph import Graph, load_graph
from pattern import PatternBuilder, Pattern

import itertools
import sys

from collections import defaultdict
from sortedcontainers import SortedSet
import bisect
import math, random
import cairo

def find_matches(LG, piece, adhesion):
    matches = defaultdict(SortedSet)
    for iu, wreach in LG.wreach_iter():
        for match in LG.match(iu, piece):
            mapped_adhesion = match.restrict_to(adhesion)
            matches[mapped_adhesion].add(iu)
    return matches

# Triangle with tail in karate
G = load_graph('example-graphs/karate.txt.gz')
H = PatternBuilder(4) \
        .add_edge(0,1).add_edge(1,2).add_edge(0,2).add_edge(0, 3) \
        .build()

print(H)

LG = G.to_lgraph()
LG.compute_wr(len(H)-1)

cand = [0, 2, 5, 8]
mapping = list(zip(cand, H))

truth = list(LG.brute_force_enumerate(H))
print("Found pattern {} times as ordered subgraph by brute force, e.g.".format(len(truth)))
print(truth[:5], "\n")

truth = set(truth)

print("Decomposing pattern:")
pieces = list(H.decompose())
assert len(pieces) == 2 # Only this special case right now

for i,piece in enumerate(pieces):
    print(i, piece)
    print("  Leaves:", piece.leaves)

print("\nCounting triangles and leaves separate:")
print("Triangles:")
triangle_matches = find_matches(LG, pieces[0], [0])
print(triangle_matches)

print("\nLeaves:")

count = 0
errors = 0
found = set()
for iu, wreach in LG.wreach_iter():
    print("\n")
    print(iu,":")

    uIN = set(LG.in_neighbours(iu))
    for iumatch in LG.match(iu, pieces[1]):
        # Restrict search to candidates that appear
        # before iu in the global order

        adhesion = iumatch.restrict_to([0])
        cands = triangle_matches[adhesion]
        cands = cands[:cands.bisect_right(iu-1)]

        print("Attempting to extend", iumatch)
        print("  adhesion:", adhesion)
        print("  candidate leaves:", cands)

        for iv in cands:
            if iv in uIN:
                continue # Abort: iu, iv are neighbours
            for ivmatch in LG.match(iv, pieces[0], partial_match=iumatch):
                count += 1
                found.add(ivmatch)
                if ivmatch in truth:
                    print(ivmatch)
                else:
                    errors += 1
                    print(">>>", ivmatch, "<<<")

            # for ivleaves in LG.match(iv, pieces[0]):
            #     ivleaves_set = set(ivleaves)
            #     if ivleaves[0] != iumatch[0]:
            #         continue

            #     overlap = (uIN & ivleaves_set) - iuleaves_set
            #     if len(overlap) == 0:
            #         count += 1
            #         match = tuple(sorted(iuleaves_set | ivleaves_set | set([iu,iv])))
            #         if match in truth:
            #             print(match)
            #         else:
            #             errors += 1
            #             print(">>>", match, "<<<")
                    

print("\n\n")
print("Total count:", count)
print("Not in truth:", errors)

missing = list(truth - found)
print("Not found:", len(missing))
print("Examples:")
print(missing[:min(len(missing), 20)])

