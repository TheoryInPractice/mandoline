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

    # Reduce to \delta(H)-core
    minDeg = min(H.degree_sequence())
    G = G.compute_core(minDeg)

    log.info("Reduced host graph to {} vertices and {} edges".format(len(G), G.num_edges()))

    LG, mapping = G.to_lgraph()
    LG.compute_wr(len(H)-1)

    marked = set() # Mark vertices that are useful for at least one patter

    for P,indexmap in H.enum_patterns():
        log.info("Current pattern {}".format(P))
        pieces = list(P.decompose())

        vertex_matches = defaultdict(lambda: defaultdict(SortedSet))
        vertex_leaf_colours = defaultdict(set)
        vertex_leaf_count = Counter()
        vertex_root_colours = defaultdict(set)        
        piece_membership = [SortedSet() for _ in range(len(P))]
        
        leaves = SortedSet() # Collect vertices that appear as leaves

        for i,piece in enumerate(pieces):
            log.info("{} {}".format(i, piece))
            log.info("  Leaves: {}".format(piece.leaves))

            # Update leaf-set            
            leaves.update(piece.leaves)

            # Update piece-lef-membership
            for j in piece.leaves:
                piece_membership[j].add(i)
            piece_membership[piece.root].add(i)

            # Compare to previous pieces. If they are equivalent up to
            # the index of the root, we can simply copy the information
            # from the previous piece since all matches will work out
            # exactly the same.
            equiv_piece = None
            for j,previous_piece in enumerate(pieces):
                if j >= i:
                    break
                if piece.root_equivalent(previous_piece):
                    equiv_piece = j
                    break

            if equiv_piece != None:
                log.info("Found equivalent piece {}!".format(j))

                # Copy info from equivalent piece
                for j in piece.leaves:
                    if equiv_piece in piece_membership[j]:
                        piece_membership[j].add(i)

                for iv in LG:
                    for index in leaves:
                        if equiv_piece in vertex_matches[iv][index]:
                            vertex_matches[iv][index].add(i)
            else:
                adh = list(sorted(set(piece.leaves) & set(pieces[-1].leaves)))

                log.info("  Adhesion: {}".format(adh))
                log.info("  Matches:")
                for m in find_matches_adh(LG, piece, adh):
                    for index,iv in m.matched_vertices():
                        # Record that v matches 'index' in piece i
                        vertex_matches[iv][index].add(i)
        

        # Determine which vertices can be matched to what pattern
        # by merging the information we obtained form each piece.
        # The idea here is that some vertice must appear as leaves in
        # multiple patterns in order to be mapped to a vertex of the 
        # whole pattern.
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
            else:
                marked.add(iv)

        log.info("Non-leaf vertices: {} ({:.1f}%)".format(non_leaves, 100.0*non_leaves/n ))
        log.info("Non-root vertices: {} ({:.1f}%)".format(non_roots, 100.0*non_roots/n ))
        log.info("Unmatched vertices: {} ({:.1f}%)".format(unmatched, 100.0*unmatched/n ))
        for l in leaves:
            log.info("{}:, {:.2f}%".format(l, 100*vertex_leaf_count[l]/len(G)))
        log.info("\n")


        # break # TEMP! For benchmark only

    log.info("{} of {} vertices marked as useful.".format(len(marked), len(LG)))