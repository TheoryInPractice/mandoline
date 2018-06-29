#!/usr/bin/env python3

from graph import Graph, load_graph
from pattern import PatternBuilder, Pattern

import argparse
import itertools
import sys

from collections import defaultdict, Counter
from sortedcontainers import SortedSet
import bisect
import math, random
import cairo

from enumerate import find_matches

import logging

log = logging.getLogger("mandoline")

def find_matches_adh(LG, piece, adhesion):
    matches = defaultdict(SortedSet)
    for iu in LG:
        for match in LG.match(iu, piece):
            yield match

if __name__ == "__main__":
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

    minDeg = min(H.degree_sequence())
    G = G.compute_core(minDeg)
    # while True:
    #     smallDegree = set()
    #     for v in G:
    #         if G.degree(v) < minDeg:
    #             smallDegree.add(v)
    #     if len(smallDegree) == 0:
    #         break

    #     for v in smallDegree:
    #         G.remove_node(v)

    log.info("Reduced host graph to {} vertices and {} edges".format(len(G), G.num_edges()))

    LG, mapping = G.to_lgraph()
    LG.compute_wr(len(H)-1)

    for P,indexmap in H.enum_patterns():
        log.info("Current pattern {}".format(P))
        pieces = list(P.decompose())

        vertex_matches = defaultdict(lambda: defaultdict(SortedSet))
        vertex_leaf_colours = defaultdict(set)
        vertex_leaf_count = Counter()
        vertex_root_colours = defaultdict(set)        
        piece_membership = [SortedSet() for _ in range(len(P))]
        leaves = SortedSet()

        if len(pieces) == 1:
            log.info("(Skipping single-piece pattern)")
            continue

        for i,piece in enumerate(pieces):
            log.info("{} {}".format(i, piece))
            log.info("  Leaves: {}".format(piece.leaves))

            leaves.update(piece.leaves)

            for j in piece.leaves:
                piece_membership[j].add(i)
            piece_membership[piece.root].add(i)

            adh = list(sorted(set(piece.leaves) & set(pieces[-1].leaves)))
            # adhesions.append(adh)

            log.info("  Adhesion: {}".format(adh))
            log.info("  Matches:")
            for m in find_matches_adh(LG, piece, adh):
                for j,v in m.matched_vertices():
                    vertex_matches[v][j].add(i)
        
        log.info("Piece membership: {}".format(piece_membership))

        for v,matches in vertex_matches.items():
            for i in matches:
                if i in leaves:
                    # The recoreded matched indices must be a superset
                    # of the piece-membership (e.g which leaf appears
                    # in what pieces)
                    if matches[i] >= piece_membership[i]:
                        vertex_leaf_colours[v].add(i)
                        vertex_leaf_count[i] += 1
                else:
                    vertex_root_colours[v].add(i)

        n = len(LG)
        non_leaves = n - len(vertex_leaf_colours)
        non_roots = n - len(vertex_root_colours)
        unmatched = 0
        for iv in LG:
            if iv not in vertex_leaf_colours and iv not in vertex_root_colours:
                unmatched += 1

        log.info("Non-leaf vertices: {} ({:.1f}%)".format(non_leaves, 100.0*non_leaves/n ))
        log.info("Non-root vertices: {} ({:.1f}%)".format(non_roots, 100.0*non_roots/n ))
        log.info("Unmatched vertices: {} ({:.1f}%)".format(unmatched, 100.0*unmatched/n ))
        for l in leaves:
            log.info("{}:, {:.2f}%".format(l, 100*vertex_leaf_count[l]/len(G)))
        log.info("\n")