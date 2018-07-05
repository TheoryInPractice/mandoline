#!/usr/bin/env python3

from graph import Graph, load_graph
from collections import Counter
from pattern import PatternBuilder, Pattern
import sys, math, random
import argparse

import logging

log = logging.getLogger("mandoline")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Decomposes a given small graph H')

    parser.add_argument('H', help='Pattern graph H')
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

    # Load pattern
    H = load_graph(args.H)
    log.info("Loaded pattern graph with {} vertices and {} edges".format(len(H), H.num_edges()))
    log.info(H)
    
    count_patterns = 0
    count_pieces = 0
    count_degenerate_pieces = 0

    all_pieces = set()
    secondary_pieces = set()

    for pattern,indexmap in H.enum_patterns():
        log.debug("%s %d", pattern, hash(pattern))
        count_patterns += 1

        pieces = pattern.decompose()
        prim_leaves = set(pieces[-1].leaves)
        for piece in pieces:
            log.debug("    %s %d", piece, hash(piece))
            all_pieces.add(piece)
            count_pieces += 1

        for piece in pieces[:-1]:
            secondary_pieces.add(piece)

            if set(piece.leaves) == prim_leaves:
                count_degenerate_pieces += 1

    log.info("")
    log.info("{} patterns".format(count_patterns))
    log.info("{} pieces".format(count_pieces))

    k = len(all_pieces)
    log.info("  of which {} ({:.1f}%) are unique".format(k,k/count_pieces * 100))

    k = len(secondary_pieces) + count_patterns
    log.info("  of which {} ({:.1f}%) are primary or unique secondary".format(k,k/count_pieces * 100))

    k = count_degenerate_pieces
    log.info("  of which {} ({:.1f}%) or degenerate".format(k,k/count_pieces * 100))