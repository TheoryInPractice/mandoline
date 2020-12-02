#!/usr/bin/env python3

from .graph import Graph, load_graph
from .pattern import PatternBuilder, Pattern
from colorama import *

import argparse
import itertools
import sys

from .helpers import CheckExt

from collections import defaultdict, Counter
from sortedcontainers import SortedSet
import bisect
import math, random

import logging

from enumerate import *
from count import *

log = logging.getLogger("mandoline")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Enumerates H in G')

    parser.add_argument('cdag', help='Counting DAG file', action=CheckExt({'dag'}))
    parser.add_argument('G', help='Host graph G')
    parser.add_argument('--debug', action='store_true')
    parser.add_argument('--no-reduction', action='store_true' )
    parser.add_argument('--quiet', action='store_true' )

    args = parser.parse_args()

    # Shut up logging.
    ch = logging.StreamHandler(sys.stdout)
    log.setLevel(logging.CRITICAL)

    # Load pattern and graph
    cdag = CDAG.load(args.cdag)
    H = cdag.target_graph()
    print("Loaded pattern graph with {} vertices and {} edges".format(len(H), H.num_edges()))
    print(H)

    G = load_graph(args.G)
    print("Loaded host graph with {} vertices and {} edges".format(len(G), G.num_edges()))

    if args.no_reduction:
        print("Skipping recution procedure because flag --no-reduction was set")
    else:
        mindeg = min(H.degree_sequence())
        G = G.compute_core(mindeg)

    LG, mapping = G.to_lgraph()
    LG.compute_wr(cdag.max_wreach)

    enum_count = 0
    enum_by_decomp = Counter()
    for P, indexmap in H.enum_patterns():
        td = TD.decompose(H, indexmap.order())
        recorder = MatchCounter()
        pieces = list(P.decompose())

        pattern_count = 0
        if len(pieces) == 1:
            count_singleton_piece(LG, pieces[0], recorder)
        else:
            assemble_pieces(LG, pieces, recorder)
        recorder.finalize()
        enum_by_decomp[td] += len(recorder)
        enum_count += len(recorder)

    count, by_decomp, _ = cdag.count(LG)

    print("\nEnumeration count: {}".format(enum_count))
    print("Count count: {}".format(count))

    print("\nBy decomposition:")
    print(f"                    Enum vs. Count {Style.RESET_ALL}")
    for td in enum_by_decomp:
        c_enum = enum_by_decomp[td]
        c = by_decomp[td]
        color = "" if c == c_enum else Fore.RED
        print(f" {color}{td.td_string():>18} {c_enum:4} vs. {c:4} {Style.RESET_ALL}")
