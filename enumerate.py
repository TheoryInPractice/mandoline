#!/usr/bin/env python3

from graph import Graph, load_graph
from pattern import PatternBuilder, Pattern

import argparse
import itertools
import sys

from collections import defaultdict
from sortedcontainers import SortedSet
import bisect
import math, random
import cairo

import logging

log = logging.getLogger(__name__)


def find_matches(LG, piece, adhesion):
    matches = defaultdict(SortedSet)
    for iu in LG:
        for match in LG.match(iu, piece):
            mapped_adhesion = match.restrict_to(adhesion)
            matches[mapped_adhesion].add(iu)
    return matches

def count_singleton_piece(LG, piece, truth):
    log.info("\nCounting singleton piece {}".format(piece))

    count = 0
    errors = 0
    matches = set()
    for iu in LG:
        log.debug("\n{} :".format(iu))

        uIN = set(LG.in_neighbours(iu))

        # Match primary piece
        for iumatch in LG.match(iu, piece):
            count += 1
            if truth != None:
                matches.add(iumatch)                
                if iumatch in truth:
                    log.debug(iumatch)
                else:
                    errors += 1
                    log.debug(">>> {} <<<".format(iumatch))
            else:
                log.debug(iumatch)     

    log.info("\nPattern count: {}\n".format(count))

    if truth:
        log.debug("False positives: {}".format(errors))

        missing = list(truth - matches)
        log.debug("Not found: {}".format(len(missing)))
        log.debug("Examples: ")
        log.debug(missing[:min(len(missing), 20)])      

        assert(len(missing) == 0)
        assert(errors == 0)          

    return count

def assemble_pieces(LG, pieces, truth):
    adhesions = []

    for i,piece in enumerate(pieces):
        log.info("{} {}".format(i, piece))
        log.info("  Leaves: {}".format(piece.leaves))

        adh = list(sorted(set(piece.leaves) & set(pieces[-1].leaves)))
        adhesions.append(adh)

        log.info("  Adhesion: {}".format(adh))

    log.debug("\nCounting secondary pieces:")
    secondary_matches = []
    for adh,piece in zip(adhesions[:-1], pieces[:-1]):
        matches = find_matches(LG, piece, adh)
        secondary_matches.append(matches)
        log.debug(piece)
        log.debug(matches)

    log.debug("\nAssembling with primary piece:")

    count = 0
    errors = 0
    matches = set()
    for iu, wreach in LG.wreach_iter():
        log.debug("\n{} :".format(iu))

        uIN = set(LG.in_neighbours(iu))

        # Match primary piece
        for iumatch in LG.match(iu, pieces[-1]):
            log.debug("Attempting to extend {}".format(iumatch))

            candidate_roots = []
            candidate_roots_indexed = []
            max_count = 1
            for i,(adh,piece) in enumerate(zip(adhesions[:-1], pieces[:-1])):
                mapped_adh = iumatch.restrict_to(adh)
                cands = secondary_matches[i][mapped_adh]
                
                # We can restrict ourselves to candidates that lie to
                # the left of iu 
                cands = cands[:cands.bisect_right(iu-1)] 

                log.debug("  piece {}: {}".format(i, piece))
                log.debug("  adhesion: {}".format(mapped_adh))
                log.debug("  candidate roots: {}".format(cands))

                candidate_roots.append(cands)
                candidate_roots_indexed.append(list(enumerate(cands)))
                max_count *= len(cands)

            assert len(candidate_roots) > 0 # Single-piece pattern case handled elsewhere 

            if max_count == 0: 
                # At least one candidate set was empty
                continue

            stack = [(0, 0, 0, iumatch)]
            while len(stack):
                log.debug("  STCK {}".format(stack))
                start_index, root_lower_bnd, piece_index, match = stack.pop()

                log.debug("  possible candidates: {}".format(candidate_roots[piece_index])) 

                lower_index = bisect.bisect_left(candidate_roots[piece_index], root_lower_bnd)
                lower_index = max(lower_index, start_index)

                log.debug("  restricted candidates: {}".format(candidate_roots[piece_index][lower_index:]))

                for i,iv in candidate_roots_indexed[piece_index][lower_index:]:
                    if iv in uIN:
                        continue # Abort: iu, iv are neighbours

                    if piece_index == len(pieces)-2:
                        # Last piece to be matched.  Every match here is
                        # a match for the whole pattern
                        for ivmatch in LG.match(iv, pieces[piece_index], partial_match=match):
                            count += 1
                            if truth != None:
                                matches.add(ivmatch)                                
                                if ivmatch in truth:
                                    log.debug(ivmatch)
                                else:
                                    errors += 1
                                    log.debug(">>> {} <<<".format(ivmatch))
                            else:
                                log.debug(ivmatch)             
                    else: 
                        candidates = []
                        for ivmatch in LG.match(iv, pieces[piece_index], partial_match=match):
                            candidates.append((0,iv,piece_index+1,ivmatch))

                        if len(candidates) > 0:
                            # Need to keep working one the `parent' match
                            stack.append((i+1,root_lower_bnd,piece_index,match))
                            stack = stack + candidates 
                            break
            log.debug("Done with extending {}\n".format(iumatch))

    log.info("\nPattern count: {}\n".format(count))

    if truth != None:
        log.debug("False positives: {}".format(errors))

        missing = list(truth - matches)
        log.debug("Not found: {}".format(len(missing)))
        log.debug("Examples:")
        log.debug(missing[:min(len(missing), 20)])

        assert(len(missing) == 0)
        assert(errors == 0)

    return count

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
    log.info("Loaded host graph with {} vertices".format(len(G)))

    LG, mapping = G.to_lgraph()
    LG.compute_wr(len(H)-1)

    if args.validate:
        log.debug("Using brute force to validate pattern count")

    # TODO: 
    # - pieces can be used as indices so adhesion/frequency should only 
    #   be computed once per piece and used in global data structure

    count = 0
    for P,indexmap in H.enum_patterns():
        log.debug("Searching pattern".format(P))
        truth = None
        if args.validate:
            truth = list(LG.brute_force_enumerate(P))
            log.debug("Found pattern {} times as ordered subgraph by brute force, e.g.".format(len(truth)))
            log.debug("{}\n".format(truth[:5]))
            truth = set(truth)

        pieces = list(P.decompose())

        if len(pieces) == 1:
            count += count_singleton_piece(LG, pieces[0], truth)
        else:
            count += assemble_pieces(LG, pieces, truth)

    log.info("\nTotal graph count: {}".format(count))