#!/usr/bin/env python3

from graph import Graph, DiGraph, DAGError, load_graph
from datastructures import Bimap, Indexmap, Interval
from pattern import PatternBuilder, Pattern
from sortedcontainers import SortedSet

import argparse
import itertools
import sys

from collections import defaultdict, Counter
from sortedcontainers import SortedSet
from itertools import permutations, product, combinations, chain
import bisect
import math, random

from tree_decompose import TD, short_str

import logging

log = logging.getLogger("mandoline")

class CDAG:
    def __init__(self):
        self.graph = None
        self.max_wreach = None
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

        self.product_edges = None
        self.subtract_edges = None
        self.dependency_dag = None

    def target_graph(self):
        return self.graph

    def count(self, LG):
        k = len(self.graph)
        counts = defaultdict(lambda: defaultdict(int))
        adhesions = defaultdict(set)

        # Traverse the topological order bases,inter,pieces in reverse
        # to ensure that dependencies are resolved.
        num_pieces = len(self.pieces)
        num_inter = len(self.decomps)
        num_bases = len(self.base_decomps)

        def count_linear(i):
            log.info("\nCounting %s %s (linear)", i, self.index[i].td_string())
            td = self.index[i]
            piece = td.to_piece(k)
            totalcount = 0
            for iu in LG:
                for match in LG.match(iu, piece):
                    # print(match)
                    for asize in self.adhesion_sizes[i]:
                        adh = match.get_adhesion(asize)
                        counts[adh][i] += 1
                        totalcount += 1
                        adhesions[i].add(adh)
            log.debug("Adhesions for {} are {}".format(i, adhesions[i]))

        def count_composite(i):
            log.info("\nCounting %s %s", i, self.index[i].td_string())
            left, right, auto_coeff = self.product_edges[i]

            td, td_left, td_right = self.index[i], self.index[left], self.index[right]
            adh_count_size = self.index[i].adhesion_size()
            for adh in (adhesions[left] & adhesions[right]):
                if len(adh) != adh_count_size:
                    continue

                c_left, c_right = counts[adh][left], counts[adh][right]
                c_subtract = sum([m*counts[adh][j] for j,m in self.subtract_edges[i].items()])
                debug_sub_str = ' - '.join(["{}*{}".format(m, counts[adh][j]) for j,m in self.subtract_edges[i].items()])
                c = (c_left * c_right) - c_subtract
                if c < 0:
                    print(f"Count for {i} bugged, negative count")
                    print(f"{c_left} x {c_right} - {debug_sub_str} = {c}")
                assert c >= 0

                if c % auto_coeff != 0:
                    print(f"Count for {i} bugged")
                    print(f"{c_left} x {c_right} - {c_subtract} = {c} not divisible by {auto_coeff}")
                assert c % auto_coeff == 0
                c //= auto_coeff

                log.info("  {} = ({} * {} - {})/{}".format(c, c_left, c_right, c_subtract, auto_coeff))

                if c == 0:
                    continue

                log.info("  Counted {} instances of {} on {}".format(c,td.td_string(),adh))
                for asize in self.adhesion_sizes[i]:
                    assert asize <= len(adh)
                    counts[adh[:asize]][i] += c
                    adhesions[i].add(adh[:asize])
                    log.info("  > Adding to {}, summing to {}".format(adh[:asize], counts[adh[:asize]][i]))


        # (num_bases+num_inter+num_pieces-1)..num_bases+num_inter
        for i in reversed(range(num_bases+num_inter,num_bases+num_inter+num_pieces)):
            assert i in self.pieces
            count_linear(i)

        # (num_bases+num_inter-1)..num_bases
        for i in reversed(range(num_bases,num_bases+num_inter)):
            assert i in self.decomps
            count_composite(i)

        # (num_bases-1)..0
        for i in reversed(range(num_bases)):
            assert i in self.base_decomps
            if self.dependency_dag.out_degree(i) == 0:
                assert self.index[i].is_linear()
                count_linear(i)
            else:
                count_composite(i)

        log.info("\nLinear counts:")
        for i in range(num_bases+num_inter,num_bases+num_inter+num_pieces):
            # We can compute the total count by fixing one adhesion size (say, the smallest) and
            # sum the counts for those adhesions only.
            adh_size = min(self.adhesion_sizes[i])
            log.info("%i %s %i", i, self.index[i].td_string(), sum([counts[adh][i] for adh in adhesions[i] if len(adh) == adh_size]))
            for adh in adhesions[i]:
                log.debug("  %s %i", adh, counts[adh][i])

        log.info("\nIntermediate counts:")
        for i in range(num_bases,num_bases+num_inter):
            # We can compute the total count by fixing one adhesion size (say, the smallest) and
            # sum the counts for those adhesions only.
            adh_size = min(self.adhesion_sizes[i])
            log.info("%i %s %i", i, self.index[i].td_string(), sum([counts[adh][i] for adh in adhesions[i] if len(adh) == adh_size]))
            for adh in adhesions[i]:
                log.debug("  %s %i", adh, counts[adh][i])

        log.info("\nFinal counts:")
        total = 0
        by_decomposition = Counter()
        for i in range(num_bases):
            pattern_count = counts[tuple()][i]
            total += pattern_count
            by_decomposition[self.index[i]] = pattern_count
            log.info("Pattern {} {} counted {} times".format(i, self.index[i].td_string(), pattern_count))
        log.info("Counted target graph {} times in host graph".format(total))
        return total, by_decomposition, counts

    @staticmethod
    def load(filename):
        res = CDAG()

        pieces = {}
        decomps = {}
        base_decomps = {}

        product_edges = {}
        subtract_edges = defaultdict(dict)
        adhesion_sizes = defaultdict(SortedSet)

        def parse_graph(line):
            var, *rest = line.split()
            if var == 'nodes':
                res.graph = Graph()
                res.graph.add_nodes(map(int, rest))
            elif var == 'edges':
                assert res.graph != None
                edges = [s.split('|') for s in rest]
                edges = [(int(x),int(y)) for x,y in edges]
                res.graph.add_edges(edges)
            elif var == 'wreach':
                res.max_wreach = int(rest[0])

        def parse_base(line):
            _id, td_str, *adh = line.split()
            _id = int(_id)
            base_decomps[_id] = TD.from_string(td_str)
            adhesion_sizes[_id].update(map(int, adh))

        def parse_composite(line):
            _id, td_str, *adh = line.split()
            _id = int(_id)
            decomps[_id] = TD.from_string(td_str)
            adhesion_sizes[_id].update(map(int, adh))

        def parse_linear(line):
            _id, td_str, *adh = line.split()
            _id = int(_id)
            pieces[_id] = TD.from_string(td_str)
            adhesion_sizes[_id].update(map(int, adh))

        def parse_edge(line):
            source, lr_string, *sub = list(line.split())
            lr_string, auto_coeff = lr_string.split('|')
            auto_coeff = int(auto_coeff)
            left, right = lr_string.split('x')

            source, left, right = int(source), int(left), int(right)

            sub = map(lambda s: tuple(s.split('|')), sub)
            sub = [(int(a), int(b)) for a,b in sub]

            assert source not in product_edges
            assert source not in subtract_edges
            assert source != left and source != right
            product_edges[source] = (left, right, auto_coeff)
            subtract_edges[source].update(sub)

        modes = {}
        modes['Base'] = parse_base
        modes['Composite'] = parse_composite
        modes['Linear'] = parse_linear
        modes['Edges'] = parse_edge
        modes['Graph'] = parse_graph
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
        res.dependency_dag = DiGraph()
        for i,td in chain(base_decomps.items(), decomps.items(), pieces.items()):
            res.index.put(i,td)
            res.dependency_dag.add_node(i)

        for s,(l,r,_) in product_edges.items():
            res.dependency_dag.add_arc(s,l)
            res.dependency_dag.add_arc(s,r)
        for s,N in subtract_edges.items():
            for t in N:
                res.dependency_dag.add_arc(s,t)
        res.dependency_dag.remove_loops() # TODO: investigate why some base decomps have loops.

        # Sanity checks: base_decomps should be sources, pieces should be sinks
        # and the indices provide a topological embedding for the graph (hence proving
        # that it is indeed a DAG).
        for i in base_decomps:
            assert res.dependency_dag.in_degree(i) == 0
            if res.dependency_dag.out_degree(i) == 0:
                assert  res.index[i].is_linear()
        for i in decomps:
            assert res.dependency_dag.in_degree(i) > 0, 'Decomp ({}) {} has in-degree zero'.format(i, res.index[i].td_string())
            assert res.dependency_dag.out_degree(i) > 0
        for i in pieces:
            assert res.dependency_dag.out_degree(i) == 0
        for i,j in res.dependency_dag.arcs():
            assert i < j

        res.base_decomps = base_decomps
        res.decomps = decomps
        res.pieces = pieces

        log.info("Base decompositions: {}--{}".format(min(base_decomps), max(base_decomps)))
        if len(decomps) > 0:
            log.info("Interm. decompositions: {}--{}".format(min(decomps), max(decomps)))
        log.info("Linear. decompositions: {}--{}".format(min(pieces), max(pieces)))

        # Store product edges and reverse lookup of what the merging
        # of a td-pair results in
        res.product_edges = product_edges
        res.subtract_edges = subtract_edges
        res.adhesion_sizes = adhesion_sizes

        return res


if __name__ == "__main__":
    from helpers import CheckExt
    parser = argparse.ArgumentParser(description='Counts subgraph in G according to counting DAG file')

    parser.add_argument('cdag', help='Counting DAG file', action=CheckExt({'dag'}))
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

    max_wreach = cdag.max_wreach
    log.info("Computing {}-wcol sets".format(max_wreach))
    LG, mapping = G.to_lgraph()
    LG.compute_wr(max_wreach)
    log.info("Done.")

    cdag.count(LG)
