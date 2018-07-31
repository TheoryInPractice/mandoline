#!/usr/bin/env python3

from graph import Graph, DiGraph, DAGError, load_graph
from datastructures import Bimap, Indexmap, Interval
from pattern import PatternBuilder, Pattern
from sortedcontainers import SortedSet

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
        self.index = None # Global index for all decompositions

        # The three types of decompositions. 'base_decomps'
        # are those we want to count (via inclusion/exclusion);
        # 'decomps' are those intermediate decompositions which we 
        # have to count first to do so (also via inclusion/exclusion) and
        # 'pieces' are _linear_ TD decompositions that can actually be counted
        # directly in the graph.
        self.base_decomps = None
        self.decomps = None
        self.pieces = None

        self.adhesion_sizes = None

        self.product_arcs = None
        self.merge_operations = None
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
            assert source != left and source != right            
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

        res.dependency_dag = DiGraph()
        for s,(l,r) in product_edges.items():
            res.dependency_dag.add_arc(s, l)
            res.dependency_dag.add_arc(s, l)
        for s,N in subtract_edges.items():
            for t in N:
                res.dependency_dag.add_arc(s,t)           
        res.dependency_dag.remove_loops() # TODO: investigate why some base decomps have loops.

        res.base_decomps = base_decomps
        res.decomps = decomps
        res.pieces = pieces

        # Store product edges and reverse lookup of what the merging
        # of a td-pair results in
        res.product_edges = product_edges

        res.merge_operations = defaultdict(dict)
        for s,(l,r) in product_edges.items():
            adh_size = len(res.index[s]._sep)
            assert (l,r) not in res.merge_operations or adh_size not in res.merge_operations[(l,r)]
            res.merge_operations[(l,r)][adh_size] = s

        # Compute adhesion sizes. For base decompositions this is simply the
        # length of the root-path until the first branching vertex; we then
        # propagate this value downwards. The same holds for intermediate 
        # decompositions; hence every decomposition might need to be counted for
        # multiple adhesion sets (and thus sizes).
        res.adhesion_sizes = defaultdict(SortedSet)
        for i,td in res.base_decomps.items():
            if td.is_linear():
                res.adhesion_sizes[i].add(0) # Simply count
            else:
                res.adhesion_sizes[i].add(len(td._sep))
            print(i,td,len(td._sep), res.dependency_dag.in_neighbours(i))

        visited = set(res.base_decomps.keys())
        print(visited)
        frontier = res.dependency_dag.out_neighbours_set(visited)
        print(frontier)

        while len(frontier) != 0:
            print(frontier)
            for i in frontier:
                for parent in res.dependency_dag.in_neighbours(i): 
                    res.adhesion_sizes[i].update(res.adhesion_sizes[parent])
                td = res.index.vertex_at(i)
                if not td.is_linear():
                    res.adhesion_sizes[i].add(len(td._sep))

            visited |= frontier
            frontier = res.dependency_dag.out_neighbours_set(visited)

        for i in res.index:
            print(i, res.adhesion_sizes[i])

        return res

    def target_graph(self):
        return self.graph

    def count(self, LG):
        pass

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
    log.info("Target graph is {}".format(H))

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
