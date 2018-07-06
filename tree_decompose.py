#!/usr/bin/env python3

from graph import Graph, load_graph
from pattern import PatternBuilder, Pattern

import argparse
import itertools
import sys

from collections import defaultdict
from sortedcontainers import SortedSet
from itertools import permutations
import bisect
import math, random
import cairo

import logging

log = logging.getLogger("mandoline")

def short_str(l):
    if len(l) == 0:
        return '.'
    return ''.join(map(str,l))

class TD:
    @staticmethod
    def decompose(G, order):
        assert(G.is_connected())
       
        return TD._decompose_rec(G, G, order, [])

    @staticmethod
    def _decompose_rec(G, subG, order, sep):
        sep = list(sep)
        depth = len(sep)
        R = Graph.from_graph(subG)
        in_neighbours = []
        for v in order:
            if v not in R:
                continue

            Nv = sorted([sep.index(u) for u in G.neighbours(v) if u in sep])
            in_neighbours.append(tuple(Nv))
            sep.append(v)
            R.remove_node(v)
            if not R.is_connected():
                break

        children = []
        for CC in R.connected_components():
            children.append(TD._decompose_rec(G, CC, order, tuple(sep)))
        children.sort(key=lambda c: c.bag)

        res = TD(tuple(sep), in_neighbours, children, depth)
        return res

    def __init__(self, sep, in_neighbours, children, depth):
        self.parent = None
        self.sep = sep
        self.bag = sep[depth:]
        self.in_neighbours = in_neighbours
        self.depth = depth
        self.children = children
        for c in self.children:
            c.parent = self

    def __hash__(self):
        res = 219787954134**(self.depth+1)
        for N in self.in_neighbours:
            res += hash(N)
        for c in self.children:
            res += hash(c)
        return res

    def __eq__(self, other):
        if other == None:
            return False

        if self.depth != other.depth or len(self.children) != len(other.children):
            return False

        if self.in_neighbours != other.in_neighbours:
            return False

        for child, other_child in zip(self.children, other.children):
            if child != other_child:
                return False

        return True

    def __repr__(self):
        if len(self.children) == 0:
            return ','.join(map(lambda N: short_str(N), self.in_neighbours))
        else:
            return ','.join(map(lambda N: short_str(N), self.in_neighbours)) + '{' + ','.join(map(str, self.children)) + '}'       

    def order_string(self):
        if len(self.children) == 0:
            return ''.join(map(str, self.bag))
        else:
            return ''.join(map(str, self.bag)) + '{' + ','.join(map(lambda c: c.order_string(), self.children)) + '}'

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Exhaustively decomposes H into tree decompositions')

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

    seen = set()
    for order in permutations(H):
        tdH = TD.decompose(H, order)
        if tdH in seen:
            continue
        seen.add(tdH)

        print("\nOrder:", ''.join(map(str,order)))
        print(H.to_lgraph(order)[0])
        
        print(tdH.order_string())
        print(tdH)
    print("\n")
    print("Computed {} tree decompositions".format(len(seen)))