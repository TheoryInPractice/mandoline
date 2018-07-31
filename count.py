#!/usr/bin/env python3

from graph import Graph, DiGraph, DAGError, load_graph
from datastructures import Bimap, Indexmap, Interval
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
        self.base_indices = None
        self.piece_indices = None
        self.inter_indices = None
        self.product_arcs = None
        self.dependency_dag = None

    @staticmethod
    def load(filename):
        res = CDAG()

        pieces = {}
        decomps = {}
        base_decomps = {}

        product_edges = {}
        subtract_edges = defaultdict(set)

        def parse_base(line):
            _id, td_str = line.split()
            _id = int(_id)
            base_decomps[_id] = TD.from_string(td_str)

        def parse_composite(line):
            _id, td_str = line.split()
            _id = int(_id)
            decomps[_id] = TD.from_string(td_str)

        def parse_linear(line):
            _id, td_str = line.split()
            _id = int(_id)
            pieces[_id] = TD.from_string(td_str)

        def parse_edge(line):
            source, left, right, *sub = list(map(int,line.split()))
            assert source not in product_edges
            assert source not in subtract_edges            
            product_edges[source] = (left, right)
            subtract_edges[source].update(sub)

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

        # Construct CDAG
        res.index = Indexmap(len(base_decomps)+len(decomps)+len(pieces))
        for i,td in chain(base_decomps.items(), decomps.items(), pieces.items()):
            res.index.put(i,td)
        res.graph = base_decomps[0].to_graph()

        res.base_indices = Interval(min(base_decomps), max(base_decomps))
        res.inter_indices = Interval(min(decomps), max(decomps))
        res.piece_indices = Interval(min(pieces), max(pieces))

        # Assert that indices are continuous and we can use intervals instead of sets
        assert len(res.base_indices) == len(base_decomps)
        assert len(res.inter_indices) == len(decomps)
        assert len(res.piece_indices) == len(pieces)

        res.dependency_dag = DiGraph()
        for s,(l,r) in product_edges.items():
            res.dependency_dag.add_arc(s, l)
            res.dependency_dag.add_arc(s, l)
        for s,N in subtract_edges.items():
            for t in N:
                res.dependency_dag.add_arc(s,t)            

        return res

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
    # log.info("Loaded counting dag graph with {} vertices and {} edges".format(len(G), G.num_edges()))
    log.info("Target graphh is {}".format(H))

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


    log.info("TODO: count")
