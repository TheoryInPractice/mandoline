#!/usr/bin/env python3

from graph import Graph, DAGError, load_graph
from datastructures import Bimap
from pattern import PatternBuilder, Pattern

import argparse
import itertools
import sys

from collections import defaultdict
from sortedcontainers import SortedSet
from itertools import permutations, product, combinations, chain
import bisect
import math, random
import cairo

from tree_decompose import TD, short_str

import logging

log = logging.getLogger("mandoline")

class CDAG:
    def __init__(self):
        self.graph = None
        self.index = None

        self.product_edges = {}
        self.subtract_edges = defaultdict(set)

    @staticmethod
    def load(filename):
        res = CDAG()

        pieces = {}
        decomps = {}
        base_decomps = {}

        product_edges = {}
        subtract_edges = defaultdict(set)
        
        def parse_base(line):
            pass

        def parse_composite(line):
            pass

        def parse_linear(line):
            pass

        def parse_edge(line):
            pass

        modes = {}
        modes['Base'] = parse_base
        modes['Composite'] = parse_composite
        modes['Linear'] = parse_linear
        modes['Edges'] = parse_edge
        with open(filename, 'r') as f:
            mode = None
            for line in f:
                line = line[:-1]
                if line[0] == '*':
                    mode = modes[line.split()[1]]
                    continue
                mode(line)

    def target_graph(self):
        return self.graph

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Counts subgraph in G according to counting DAG file')

    parser.add_argument('cdag', help='Counting DAG file')
    parser.add_argument('G', help='Host graph G')
    parser.add_argument('--debug', action='store_true')
    parser.add_argument('--quiet', action='store_true' )
    parser.add_argument('--no-reduction', action='store_true' )    
    parser.add_argument('--output', help='Output file for counting DAG' )

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
    cdag = CDAG.load(args.cdag)
    H = cdag.target_graph()

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

