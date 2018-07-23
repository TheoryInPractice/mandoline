#!/usr/bin/env python3

from graph import Graph, load_graph
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

class Recorder:
    def __init__(self):
        self.pieces = set()
        self.decomps = set()

    def count_linear(self, td):
        assert(td.is_linear())
        if td in self.pieces:
            print("SEEN Piece", td.td_string())
        else:
            print("NEW Piece", td.td_string())
            self.pieces.add(td)

    def count_recursive(self, td):
        assert(not td.is_linear())
        if td in self.decomps:
            print("SEEN Decomp", td.td_string())
            return True
        else:
            print("NEW Decomp", td.td_string())
            self.decomps.add(td)
            return False

    def report(self):
        print("Recorded {} linear pieces and {} decompositions.".format(len(self.pieces), len(self.decomps)))

        pieces_str = sorted(map(lambda p: p.td_string(), self.pieces))
        for p in pieces_str:
            print(p)

def powerset_nonempty(iterable):
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(1, len(s)+1))

def suborder(order, vertices):
    return list(filter(lambda s: s in vertices, order))

def simulate_count(R, H, td):
    _simulate_count_rec(R, H, td, 0)

def _simulate_count_rec(R, H, td, depth):
    prefix = " "*(4*depth)
    # print(prefix+"We want to count", td, " for graph", list(H.edges()))
    print(prefix+"We want to count", td.td_string(), "(", td, ")")

    split_depth = td.adhesion_size()
    splits = list(td.split())
    if len(splits) == 1:
        print(prefix+"This decomposition is linear, we simply count it!")
        R.count_linear(td)
        return

    already_known = R.count_recursive(td)
    if already_known:
        return

    print(prefix+"The decomposition branches at depth", split_depth)

    order_prefix = td._sep[:split_depth]
    current = splits[0]
    print(prefix+"We first count the leftmost piece", current)
    _simulate_count_rec(R, H, current, depth+1)

    print(prefix+"Now we fold-count with the remaining pieces.")

    orders = list(td.orders())

    for td_next in splits[1:]:
        print(prefix+"The next piece is", td_next, "and we first count it.")
        _simulate_count_rec(R, H, td_next, depth+1)
        previous = current
        current = current.merge(td_next, split_depth)
        print(prefix+"The initial count of", current, "is the count of", previous, "times the count of", td_next)   

        old_nodes = previous.nodes()
        new_nodes = td_next.nodes()
        common_nodes = old_nodes & new_nodes
        old_nodes -= common_nodes
        new_nodes -= common_nodes
        joint_nodes = old_nodes | new_nodes | set(td._sep)

        print(prefix+"To account for non-induced instance, edges between", old_nodes, "and", new_nodes, "need to be considered" ) 
        potential_edges = list(product(old_nodes, new_nodes))

        print(prefix+"We subtract the results of the following counts:")
        seen = set()
        for oo in orders:
            o = suborder(oo, joint_nodes)
            for edges in powerset_nonempty(potential_edges):
                assert len(o) > 0
                HH = H.subgraph(joint_nodes)
                for u,v in edges:
                    assert u in HH and v in HH
                    HH.add_edge(u,v)
                print(prefix+"Decompositing",list(HH.edges()),"along order", o)                    
                tdHH = TD.decompose(HH, o)
                if tdHH in seen:
                    continue
                _simulate_count_rec(R, HH, tdHH, depth+1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Counts H in G')

    parser.add_argument('H', help='Pattern graph H')
    # parser.add_argument('G', help='Host graph G')
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

    # G = load_graph(args.G)
    # log.info("Loaded host graph with {} vertices and {} edges".format(len(G), G.num_edges()))

    # if args.no_reduction:
    #     log.info("Skipping recution procedure because flag --no-reduction was set")
    # else:
    #     mindeg = min(H.degree_sequence()) 
    #     log.info("Computing {}-core of host graph".format(mindeg))
    #     G = G.compute_core(mindeg)
    #     log.info("Reduced host graph to {} vertices and {} edges".format(len(G), G.num_edges()))

    # log.info("Computing {}-wcol sets".format(len(H)-1))
    # LG, mapping = G.to_lgraph()
    # LG.compute_wr(len(H)-1)
    # log.info("Done.")

    seen = set()
    R = Recorder()
    for order in permutations(H):
        tdH = TD.decompose(H, order)
        if tdH in seen:
            continue

        print("The decomposition {} represents the following orders:".format(tdH.order_string()))
        for o in tdH.orders():
            print("  ", short_str(o))

        seen.add(tdH)
        print()

        simulate_count(R, H, tdH)
        print()

    print("\n")
    print("Computed {} tree decompositions".format(len(seen)))

    R.report()