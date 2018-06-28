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

def find_matches(LG, piece, adhesion):
    matches = defaultdict(SortedSet)
    for iu in LG:
        for match in LG.match(iu, piece):
            mapped_adhesion = match.restrict_to(adhesion)
            matches[mapped_adhesion].add(iu)
    return matches

def count_singleton_piece(LG, piece, truth):
    print("\nCounting singleton piece")

    count = 0
    errors = 0
    matches = set()
    for iu in LG:
        print("\n")
        print(iu,":")

        uIN = set(LG.in_neighbours(iu))

        # Match primary piece
        for iumatch in LG.match(iu, piece):
            count += 1
            matches.add(iumatch)
            if truth != None:
                if iumatch in truth:
                    print(iumatch)
                else:
                    errors += 1
                    print(">>>", iumatch, "<<<")    
            else:
                print(iumatch)     

    print("\n\n")
    print("Total count:", count)

    if truth:
        print("False positives:", errors)

        missing = list(truth - matches)
        print("Not found:", len(missing))
        print("Examples:")
        print(missing[:min(len(missing), 20)])      

        assert(len(missing) == 0)
        assert(errors == 0)          

def assemble_pieces(LG, pieces, truth):
    adhesions = []

    for i,piece in enumerate(pieces):
        print(i, piece)
        print("  Leaves:", piece.leaves)

        adh = list(sorted(set(piece.leaves) & set(pieces[-1].leaves)))
        adhesions.append(adh)

        print("  Adhesion:", adh)

    print("\nCounting secondary pieces:")
    secondary_matches = []
    for adh,piece in zip(adhesions[:-1], pieces[:-1]):
        matches = find_matches(LG, piece, adh)
        secondary_matches.append(matches)
        print(piece)
        print(matches)

    print("\nAssembling with primary piece:")

    count = 0
    errors = 0
    matches = set()
    for iu, wreach in LG.wreach_iter():
        print("\n")
        print(iu,":")

        uIN = set(LG.in_neighbours(iu))

        # Match primary piece
        for iumatch in LG.match(iu, pieces[-1]):
            print("Attempting to extend", iumatch)

            candidate_roots = []
            candidate_roots_indexed = []
            max_count = 1
            for i,(adh,piece) in enumerate(zip(adhesions[:-1], pieces[:-1])):
                mapped_adh = iumatch.restrict_to(adh)
                cands = secondary_matches[i][mapped_adh]
                
                # We can restrict ourselves to candidates that lie to
                # the left of iu 
                cands = cands[:cands.bisect_right(iu-1)] 

                print("  piece {}: ".format(i), piece)
                print("  adhesion:", mapped_adh)
                print("  candidate roots:", cands)

                candidate_roots.append(cands)
                candidate_roots_indexed.append(list(enumerate(cands)))
                max_count *= len(cands)

            assert len(candidate_roots) > 0 # Single-piece pattern case handled elsewhere 

            if max_count == 0: 
                # At least one candidate set was empty
                continue

            stack = [(0, 0, 0, iumatch)]
            while len(stack):
                print("  STCK",stack)
                start_index, root_lower_bnd, piece_index, match = stack.pop()

                print("  possible candidates:", candidate_roots[piece_index]) 

                lower_index = bisect.bisect_left(candidate_roots[piece_index], root_lower_bnd)
                lower_index = max(lower_index, start_index)

                print("  restricted candidates:", candidate_roots[piece_index][lower_index:])

                for i,iv in candidate_roots_indexed[piece_index][lower_index:]:
                    if iv in uIN:
                        continue # Abort: iu, iv are neighbours

                    if piece_index == len(pieces)-2:
                        # Last piece to be matched.  Every match here is
                        # a match for the whole pattern
                        for ivmatch in LG.match(iv, pieces[piece_index], partial_match=match):
                            count += 1
                            matches.add(ivmatch)
                            if truth != None:
                                if ivmatch in truth:
                                    print(ivmatch)
                                else:
                                    errors += 1
                                    print(">>>", ivmatch, "<<<")       
                            else:
                                print(ivmatch)             
                    else: 
                        candidates = []
                        for ivmatch in LG.match(iv, pieces[piece_index], partial_match=match):
                            candidates.append((0,iv,piece_index+1,ivmatch))

                        if len(candidates) > 0:
                            # Need to keep working one the `parent' match
                            stack.append((i+1,root_lower_bnd,piece_index,match))
                            stack = stack + candidates 
                            break
            print("Done with extending", iumatch, "\n")

    print("\n\n")
    print("Total count:", count)

    if truth != None:
        print("False positives:", errors)

        missing = list(truth - matches)
        print("Not found:", len(missing))
        print("Examples:")
        print(missing[:min(len(missing), 20)])

        assert(len(missing) == 0)
        assert(errors == 0)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Enumerates H in G')

    parser.add_argument('H', help='Pattern graph H')
    parser.add_argument('G', help='Host graph G')
    parser.add_argument('--validate', action='store_true')

    args = parser.parse_args()

    H = load_graph(args.H)
    print("Loaded pattern graph with {} vertices and {} edges".format(len(H), H.num_edges()))
    print(H)

    G = load_graph(args.G)
    print("Loaded host graph with {} vertices".format(len(G)))

    LG, mapping = G.to_lgraph()
    LG.compute_wr(len(H)-1)

    # TODO: 
    # - pieces can be used as indices so adhesion/frequency should only 
    #   be computed once per piece and used in global data structure

    for P,indexmap in H.enum_patterns():
        print("Searching pattern", P)
        truth = None
        if args.validate:
            truth = list(LG.brute_force_enumerate(P))
            print("Found pattern {} times as ordered subgraph by brute force, e.g.".format(len(truth)))
            print(truth[:5], "\n")
            truth = set(truth)
        else:
            print("Graph to large to brute force")

        pieces = list(P.decompose())

        if len(pieces) == 1:
            count_singleton_piece(LG, pieces[0], truth)
        else:
            assemble_pieces(LG, pieces, truth)
