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

log = logging.getLogger("mandoline")

class MatchCounter:
    def __init__(self):
        self.count = 0

    def record(self, match):
        self.count += 1
        log.debug(match)

    def finalize(self):
        pass

    def __len__(self):
        return self.count

class MatchValidator:
    def __init__(self, truth):
        self.truth = truth
        self.errors = 0
        self.matches = set()

    def record(self, match):
        self.matches.add(match)                
        if match in self.truth:
            log.debug(match)
        else:
            self.errors += 1
            log.debug(">>> %s <<<", match)

    def finalize(self):
        log.debug("False positives: %d", self.errors)

        missing = list(self.truth - self.matches)
        log.debug("Not found: %d", len(missing))
        log.debug("Examples: ")
        log.debug(missing[:min(len(missing), 20)])      

        assert(len(missing) == 0)
        assert(self.errors == 0)      

    def __len__(self):
        return len(self.matches)    

def complete_rooted_match(LG, r, pieces, partial_match):
    """
        Complete a partial match in which all roots have
        been already matched.
    """
    if r >= len(pieces):
        assert(partial_match.is_complete())
        yield partial_match
        return

    for match in LG.complete(pieces[r], partial_match):
        for res in complete_rooted_match(LG, r+1, pieces, match):
            yield res

def independent_root_sets(LG, root_candidates):
    """
        Given a r lists of candidates (root_candidates) and a linear graph LG,
        compute all ordered independent sets of size r that pick exactly one 
        vertex from each set.
    """
    assert(len(root_candidates) >= 2)

    if len(root_candidates) == 2:
        for res in _independent_root_sets_binary(LG, *root_candidates):
            yield res
    elif len(root_candidates) == 3:
        for res in _independent_root_sets_ternary(LG, *root_candidates):
            yield res
    else:
        for iu in root_candidates[0]:
            for res in _independent_root_sets_rec(LG, [iu], 1, root_candidates):
                yield res

def _independent_root_sets_binary(LG, root_candsA, root_candsB):
    posB = 0
    nB = len(root_candsB)

    for vA in root_candsA:
        posB = bisect.bisect_right(root_candsB, vA, lo=posB)

        for vB in root_candsB[posB:]:
            if not LG.adjacent_ordered(vA, vB):
                yield [vA, vB]

def _independent_root_sets_ternary(LG, root_candsA, root_candsB, root_candsC):
    posB = 0
    nB = len(root_candsB)
    nC = len(root_candsC)

    for vA in root_candsA:
        posB = bisect.bisect_right(root_candsB, vA, lo=posB)
        posC = 0
        
        for vB in root_candsB[posB:]:
            if LG.adjacent_ordered(vA, vB):
               continue

            posC = bisect.bisect_right(root_candsC, vB, lo=posC)                
            if posC >= nC:
                break
            for vC in root_candsC[posC:]:
                if not LG.adjacent_ordered(vA, vC) and not LG.adjacent_ordered(vB, vC):
                    yield [vA, vB, vC]

def _independent_root_sets_rec(LG, selection, r, root_candidates):
    if r == len(root_candidates):
        yield selection
        return

    lower = bisect.bisect_right(root_candidates[r], selection[-1])
    for iu in root_candidates[r][lower:]:
        Nu = LG.in_neighbours(iu) # Selection so far lies to the left of iu
        if len(Nu & set(selection)) > 0:
            continue # Not independent
        for res in _independent_root_sets_rec(LG, selection+[iu], r+1, root_candidates):
            yield res

def find_matches(LG, piece, adhesion):
    matches = defaultdict(SortedSet)
    for iu in LG:
        for match in LG.match(iu, piece):
            mapped_adhesion = match.restrict_to(adhesion)
            matches[mapped_adhesion].add(iu)
    return matches

def count_singleton_piece(LG, piece, recorder):
    log.info("\nCounting singleton piece {}".format(piece))

    for iu in LG:
        log.debug("\n%d :", iu)

        uIN = set(LG.in_neighbours(iu))

        # Match primary piece
        for match in LG.match(iu, piece):
            recorder.record(match)

def assemble_pieces(LG, pieces, recorder):
    prim_piece = pieces[-1]
    sec_pieces = pieces[:-1]

    adhesions = []
    root_indices = []

    log.info("Primary {}".format(prim_piece))
    log.info("  Leaves: {}".format(prim_piece.leaves))

    # Compute adhesion sets and root indices  of every secondary piece
    for i,piece in enumerate(sec_pieces):
        log.info("{} {}".format(i, piece))
        log.info("  Leaves: {}".format(piece.leaves))

        adh = list(sorted(set(piece.leaves) & set(prim_piece.leaves)))
        adhesions.append(adh)
        root_indices.append(piece.root)

        log.info("  Adhesion: {}".format(adh))

    # Supplement with info from primary piece
    root_indices.append(prim_piece.root)

    log.info("Roots: %s", root_indices)

    # Collect boundaries for primary matches
    log.debug("Computing primary matches for piece %s", prim_piece)
    prim_matches = defaultdict(SortedSet)
    for iu in LG:
        for match in LG.match(iu, prim_piece):
            boundary = match.restrict_to(prim_piece.leaves)
            prim_matches[boundary].add(iu)
            log.debug("  Primary match %s", boundary)

    # Determine interesting secondary piece boundaries from
    # the collected primary boundaries and note for each vertex
    # which matches are `allowed' (meaning they might lead to 
    # a full match)
    sec_matches = defaultdict(SortedSet)
    filtered_leaves = set(prim_piece.leaves)
    allowed_matches = defaultdict(set)
    log.debug("Computing secondary bondaries")
    for prim_boundary in prim_matches:
        for index, iv in prim_boundary.matched_vertices():
            assert(index in filtered_leaves)
            allowed_matches[iv].add(index)
        for i, sec_adh in enumerate(adhesions):
            sec_boundary = prim_boundary.restrict_to(sec_adh)
            if sec_boundary not in sec_matches:
                log.debug("  Sec. boundary for piece %d: %s", i, sec_boundary)
            sec_matches[sec_boundary] # This adds the key to the defaultdict


    # Collect secondary matches
    for i, (piece, sec_adh) in enumerate(zip(sec_pieces,adhesions)):
        log.debug("Computing secondary matches for piece %s", piece)
        for iu in LG:
            for match in LG.match(iu, piece, filtered_leaves=filtered_leaves, allowed_matches=allowed_matches):
                boundary = match.restrict_to(sec_adh)
                if boundary in sec_matches:
                    log.debug("  Match: %s", boundary)
                    sec_matches[boundary].add(iu)
                else:
                    log.debug("  Dismissing match: %s, bondary %s not found", match, boundary)

    # Assemble!
    for prim_boundary in prim_matches:
        # Collect candidate roots
        cand_roots = []
        max_count = 1
        for i, (piece, sec_adh) in enumerate(zip(sec_pieces,adhesions)):
            sec_boundary = prim_boundary.restrict_to(sec_adh)
            cand_roots.append(sec_matches[sec_boundary])
            max_count *= len(cand_roots[-1]) 
        cand_roots.append(prim_matches[prim_boundary])
        max_count *= len(cand_roots[-1])

        # Early out: at least one piece cannot be matched
        if max_count == 0:
            continue

        log.debug("Matches for boundary %s", prim_boundary)
        for i,candidates in enumerate(cand_roots):
            log.debug("  (%d) %s",i, candidates)
        for root_set in independent_root_sets(LG, cand_roots):
            extension = list(zip(root_indices, root_set))
            partial_match = prim_boundary.extend_multiple(extension)

            if partial_match == None:
                continue

            if partial_match.is_complete():
                recorder.record(partial_match)
                log.debug(">>> %s %s", extension, partial_match)                
            else:
                # TODO!
                log.debug(">>> Completing partial match %s %s", extension, partial_match)  
                for match in complete_rooted_match(LG, 0, pieces, partial_match):
                    log.debug(">>> %s", match)
                    recorder.record(match)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Enumerates H in G')

    parser.add_argument('H', help='Pattern graph H')
    parser.add_argument('G', help='Host graph G')
    parser.add_argument('--validate', action='store_true')
    parser.add_argument('--debug', action='store_true')
    parser.add_argument('--no-reduction', action='store_true' )
    parser.add_argument('--quiet', action='store_true' )

    args = parser.parse_args()

    # Set up logging
    ch = logging.StreamHandler(sys.stdout)
    if args.quiet:
        # Mute on top level, don't add handler
        log.setLevel(logging.CRITICAL)
    elif args.debug:
        ch.setLevel(logging.DEBUG)
        log.setLevel(logging.DEBUG)
        log.addHandler(ch)
    else:
        ch.setLevel(logging.INFO)
        log.setLevel(logging.INFO)
        log.addHandler(ch)

    # Load pattern and graph
    H = load_graph(args.H)
    log.info("Loaded pattern graph with {} vertices and {} edges".format(len(H), H.num_edges()))
    log.info(H)

    G = load_graph(args.G)
    log.info("Loaded host graph with {} vertices and {} edges".format(len(G), G.num_edges()))

    if args.no_reduction:
        log.info("Skipping recution procedure because flag --no-reduction was set")
    else:
        mindeg = min(H.degree_sequence()) 
        log.info("Computing {}-core of host graph".format(mindeg))
        G = G.compute_core(mindeg)
        log.info("Reduced host graph to {} vertices and {} edges".format(len(G), G.num_edges()))

    log.info("Computing {}-wcol sets".format(len(H)-1))
    LG, mapping = G.to_lgraph()
    LG.compute_wr(len(H)-1)
    log.info("Done.")

    if args.validate:
        log.debug("Using brute force to validate pattern count")

    # TODO: 
    # - pieces can be used as indices so adhesion/frequency should only 
    #   be computed once per piece and used in global data structure

    count = 0
    for P,indexmap in H.enum_patterns():
        log.info("Searching pattern {}".format(P))

        recorder = MatchCounter()
        if args.validate:
            truth = list(LG.brute_force_enumerate(P))
            log.info("Found pattern {} times as ordered subgraph by brute force, e.g.".format(len(truth)))
            log.info("{}\n".format(truth[:5]))
            truth = set(truth)
            recorder = MatchValidator(truth)

        pieces = list(P.decompose())

        pattern_count = 0
        if len(pieces) == 1:
            count_singleton_piece(LG, pieces[0], recorder)
        else:
            assemble_pieces(LG, pieces, recorder)
        recorder.finalize()
        count += len(recorder)
        log.info("\nPattern count: %d\n", len(recorder))

    # Always print the final count, even in --quiet mode.
    print("\nTotal graph count: {}".format(count))