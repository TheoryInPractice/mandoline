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
        self.pieces = defaultdict(set)
        self.decomps = defaultdict(set)

    def count_linear(self, td, adh_size):
        assert(td.is_linear())
        if td not in self.pieces:
            log.info("Found new piece %s", td.td_string())
        self.pieces[td].add(adh_size)

    def count_recursive(self, td, adh_size):
        assert(not td.is_linear())
        if td in self.decomps and adh_size in self.decomps[td]:
            return True

        self.decomps[td].add(adh_size)
        log.info("Found new decomp/adh size combination %s, %s", td.td_string(), adh_size)
        return False

    def report(self):
        log.info("Recorded %d linear pieces and %d decompositions.", len(self.pieces), len(self.decomps))

def powerset_nonempty(iterable):
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(1, len(s)+1))

def suborder(order, vertices):
    return tuple(filter(lambda s: s in vertices, order))

def simulate_count(R, H, td):
    _simulate_count_rec(R, 0, H, td, 0)

def _simulate_count_rec(R, adh_size, H, td, depth):
    prefix = " "*(4*depth)
    log.debug("%sWe want to count %s (%s) for adhesion size %s", prefix, td.td_string(), str(td), adh_size)

    split_depth = td.adhesion_size()
    splits = list(td.split())
    if len(splits) == 1:
        log.debug("%sThis decomposition is linear, we simply count it!", prefix)
        R.count_linear(td, adh_size)
        return

    already_known = R.count_recursive(td, adh_size)
    if already_known:
        return

    log.debug("%sThe decomposition branches at depth %d", prefix, split_depth)

    order_prefix = td._sep[:split_depth]
    td_current = splits[0]
    log.debug("%sWe first count the leftmost piece %s", prefix, td_current)
    _simulate_count_rec(R, split_depth, H, td_current, depth+1)

    log.debug("%sNow we fold-count with the remaining pieces.", prefix)

    orders = list(td.orders())
    separator = frozenset(td._sep)

    for td_next in splits[1:]:
        log.debug("%sThe next piece is %s and we first count it.", prefix, td_next)
        _simulate_count_rec(R, split_depth, H, td_next, depth+1)
        td_previous = td_current
        td_current = td_current.merge(td_next, split_depth)
        log.debug("%sThe initial count of %s is the count of %s times the count of %s", prefix, td_current, td_previous, td_next)   

        for (HH, tdHH) in enumerate_defects(H, orders, separator, td_previous, td_next, depth):
            _simulate_count_rec(R, split_depth, HH, tdHH, depth+1)


def enumerate_defects(H, orders, separator, decompA, decompB, depth):
    """
        Enumerates all possible graphs with td decompositions that contain
        induced subgraph that decompose into decompA, decompB (with the joint
        separator 'separator').
    """
    prefix = " "*(4*depth)
    old_nodes = decompA.nodes()
    new_nodes = decompB.nodes()
    common_nodes = old_nodes & new_nodes
    old_nodes -= common_nodes
    new_nodes -= common_nodes
    joint_nodes = old_nodes | new_nodes | separator

    # TODO: enumerate _overlaps_

    log.debug("%sTo account for non-induced instance, edges between %s and %s need to be considered", prefix, old_nodes, new_nodes ) 
    potential_edges = list(product(old_nodes, new_nodes))

    log.debug("%sWe subtract the results of the following counts:", prefix)
    seen = set()
    seen_order = set() 
    for oo in orders:
        # Since we are constructing suborders, some might
        # appear twice.
        o = suborder(oo, joint_nodes)
        if o in seen_order:
            continue
        seen_order.add(o)

        for edges in powerset_nonempty(potential_edges):
            assert len(o) > 0
            HH = H.subgraph(joint_nodes)
            for u,v in edges:
                assert u in HH and v in HH
                HH.add_edge(u,v)
            log.debug("%sDecompositing %s along order %s", prefix,list(HH.edges()), o)                    
            tdHH = TD.decompose(HH, o)
            if tdHH in seen:
                continue
            yield (HH, tdHH)


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

        log.info("The decomposition %s represents the following orders:", tdH.order_string())
        for o in tdH.orders():
            log.info("  " + short_str(o))

        seen.add(tdH)
        log.info("")

        simulate_count(R, H, tdH)
        log.info("")

    log.info("\n")
    log.info("Computed %d tree decompositions", len(seen))

    R.report()