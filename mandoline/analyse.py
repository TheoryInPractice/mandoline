#!/usr/bin/env python3

from .graph import Graph, load_graph
from .pattern import PatternBuilder, Pattern

import argparse
import itertools
import sys

from collections import defaultdict, Counter
from sortedcontainers import SortedSet
import bisect
import math, random


from .enumerate import find_matches

import logging

log = logging.getLogger("mandoline")

def find_matches_adh(LG, piece, adhesion):
    matches = defaultdict(SortedSet)
    for iu in LG:
        for match in LG.match(iu, piece):
            yield match

def main():
    parser = argparse.ArgumentParser(description='Enumerates H in G')

    parser.add_argument('H', help='Pattern graph H')
    parser.add_argument('G', help='Host graph G')
    parser.add_argument('--validate', action='store_true')
    parser.add_argument('--debug', action='store_true')

    args = parser.parse_args()

    # Set up logging
    ch = logging.StreamHandler(sys.stdout)
    if args.debug:
        ch.setLevel(logging.DEBUG)
    else:
        ch.setLevel(logging.INFO)
    log.addHandler(ch)
    log.setLevel(logging.DEBUG)

    # Load pattern and graph
    H = load_graph(args.H)
    log.info("Loaded pattern graph with {} vertices and {} edges".format(len(H), H.num_edges()))
    log.info(H)

    G = load_graph(args.G)
    log.info("Loaded host graph with {} vertices and {} edges".format(len(G), G.num_edges()))

    # Reduce to \delta(H)-core
    minDeg = min(H.degree_sequence())
    G = G.compute_core(minDeg)

    log.info("Reduced host graph to {} vertices and {} edges".format(len(G), G.num_edges()))

    LG, mapping = G.to_lgraph()
    LG.compute_wr(len(H)-1)

    marked = set() # Mark vertices that are useful for at least one patter
    
    pieces = set()
    for P,indexmap in H.enum_patterns():
        pieces.update(P.decompose())

    log.info("Computed all {} pieces".format(len(pieces)))

    matches = [defaultdict(SortedSet) for _ in range(len(pieces))]
    total = 0
    for i, piece in enumerate(pieces):
        log.info("Computed matches for piece {}".format(piece))
        count = 0
        for iu in LG:
            for match in LG.match(iu, piece):
                boundary = match.restrict_to(piece.leaves)
                matches[i][boundary].add(iu)
                count += 1
        total += count
        log.info("  Found {} matches".format(count))
    log.info("Collected {} total matches".format(total))


if __name__ == "__main__":
    main()
