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

        self.product_edges = None
        self.subtract_edges = None
        self.merge_operations = None
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
            print()
            print("Counting", i, self.index[i].td_string(), "(linear)")
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
            print("Adhesions for {} are {}".format(i, adhesions[i]))

        def count_composite(i):
            print()
            print("Counting", i, self.index[i].td_string())
            left, right = self.product_edges[i]

            if len(self.subtract_edges[i]) > 0:
                subtract_ids, subtract_multis = zip(*self.subtract_edges[i])
            else:
                subtract_ids, subtract_multis = [], []

            td, td_left, td_right = self.index[i], self.index[left], self.index[right]
            td_subtract = self.index.vertices_at(subtract_ids)
            # td_subtract_str = '; '.join(map(lambda t: t.td_string(), td_subtract))
            td_subtract_str = ' '.join([ '- {}#{}'.format(m,self.index.vertex_at(x).td_string()) for x,m in zip(subtract_ids, subtract_multis)])
            subtract_str = ' '.join([ '- {}#{}'.format(m,x) for x,m in zip(subtract_ids, subtract_multis)])
            print("#{} = #{} x #{} {}".format(i, left, right, subtract_str))
            print("#{} = #{} x #{} {}".format(td.td_string(), td_left.td_string(), td_right.td_string(), td_subtract_str))

            print("Left pattern found for {}".format(adhesions[left]))
            print("Right pattern found for {}".format(adhesions[right]))
            print("  Intersection contains {}".format(adhesions[left] & adhesions[right]))

            adh_count_size = self.index[i].adhesion_size()
            for adh in (adhesions[left] & adhesions[right]):
                if len(adh) != adh_count_size:
                    continue

                c_left, c_right = counts[adh][left], counts[adh][right]
                c_subtract = sum([m*counts[adh][j] for j,m in self.subtract_edges[i]])
                c = c_left * c_right - c_subtract

                if left == right:
                    assert c % 2 == 0
                    c //= 2
                assert c >= 0

                print("  {} = {} * {} - {}".format(c, c_left, c_right, c_subtract))

                if c == 0:
                    continue

                print("  Counted {} instances of {} on {}".format(c,td.td_string(),adh))
                for asize in self.adhesion_sizes[i]:
                    assert asize <= len(adh)
                    counts[adh[:asize]][i] += c
                    adhesions[i].add(adh[:asize])
                    print("  > Adding to {}, summing to {}".format(adh[:asize], counts[adh[:asize]][i]))


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

        print()
        print("Linear counts:")
        for i in range(num_bases+num_inter,num_bases+num_inter+num_pieces):
            # We can compute the total count by fixing one adhesion size (say, the smallest) and
            # sum the counts for those adhesions only.
            adh_size = min(self.adhesion_sizes[i])
            print(i, self.index[i].td_string(), sum([counts[adh][i] for adh in adhesions[i] if len(adh) == adh_size]))
            # for adh in adhesions[i]:
            #     print("  ", adh, counts[adh][i])

        print()
        print("Intermediate counts:")
        for i in range(num_bases,num_bases+num_inter):
            # We can compute the total count by fixing one adhesion size (say, the smallest) and
            # sum the counts for those adhesions only.
            adh_size = min(self.adhesion_sizes[i])
            print(i, self.index[i].td_string(), sum([counts[adh][i] for adh in adhesions[i] if len(adh) == adh_size]))
            for adh in adhesions[i]:
                print("  ", adh, counts[adh][i])

        print()
        print("Final counts:")
        total = 0
        for i in range(num_bases):
            pattern_count = counts[tuple()][i]
            total += pattern_count
            print("Pattern {} {} counted {} times".format(i, self.index[i].td_string(), pattern_count))
        print("Counted target graph {} times in host graph".format(total))

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
            source, left, right, *sub = list(line.split())
            source, left, right = int(source), int(left), int(right)

            sub = map(lambda s: tuple(s.split('|')), sub)
            sub = [(int(a), int(b)) for a,b in sub]

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
        res.dependency_dag = DiGraph()
        for i,td in chain(base_decomps.items(), decomps.items(), pieces.items()):
            res.index.put(i,td)
            res.dependency_dag.add_node(i)
        res.graph = base_decomps[0].to_graph()

        conflicts = {}
        for s,(l,r) in product_edges.items():
            res.dependency_dag.add_arc(s, l)
            res.dependency_dag.add_arc(s, r)
            if (l,r) in conflicts:
                td_left, td_right = res.index[l], res.index[r]
                td_res, td_conf = res.index[s], res.index[conflicts[(l,r)]]
                print(f"Conflict: {td_left} + {td_right} = {td_res} AND {td_conf}")
                print(f"          {td_left.td_string()} + {td_right.td_string()} = {td_res.td_string()} AND {td_conf.td_string()}")
                assert False
            conflicts[(l,r)] = conflicts[(r,l)] = s
        for s,N in subtract_edges.items():
            for t,_ in N:
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

        print("Base decompositions: {}--{}".format(min(base_decomps), max(base_decomps)))
        print("Interm. decompositions: {}--{}".format(min(decomps), max(decomps)))
        print("Linear. decompositions: {}--{}".format(min(pieces), max(pieces)))

        # Store product edges and reverse lookup of what the merging
        # of a td-pair results in
        res.product_edges = product_edges
        res.subtract_edges = subtract_edges

        res.merge_operations = defaultdict(dict)
        for s,(l,r) in product_edges.items():
            adh_size = len(res.index[s]._sep)
            assert (l,r) not in res.merge_operations or adh_size not in res.merge_operations[(l,r)]
            res.merge_operations[(l,r)][adh_size] = s

        # TODO: at this point we could also compute and propagate the wcol-distances
        # of decompositions, thus pruning the search space for couting pieces.
        # TODO: this should probably all be done during construction and stored in
        # the .dag file.

        # Compute adhesion sizes. For base decompositions this is simply zero, meaning
        # the global counter. For all other decompositions, the values are derived from
        # the products they are involved in: if decompositions A and B merge together
        # into decomposition C with a root-path of length r, then A and B need to be
        # counted for adhesions of length r.
        res.adhesion_sizes = defaultdict(SortedSet)
        for i,td in res.base_decomps.items():
            res.adhesion_sizes[i].add(0) # Simply count

        # Now compute adhesion sizes for all other decompositions
        visited = set(res.base_decomps)
        frontier = res.dependency_dag.out_neighbours_set(visited)

        while len(frontier) != 0:
            for i in frontier:
                for parent in res.dependency_dag.in_neighbours(i):
                    parent_adhesion = res.index[parent].adhesion_size()
                    if parent_adhesion > res.index[i].adhesion_size():
                        # This issue should be solved since commit a8202d9.
                        print("Decomp: ", res.index[i].td_string())
                        print("Parent: ", res.index[parent].td_string())
                        assert False
                    res.adhesion_sizes[i].add(parent_adhesion)


            visited |= frontier
            frontier = res.dependency_dag.out_neighbours_set(visited)

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

    log.info("Computing {}-wcol sets".format(len(H)-1))
    LG, mapping = G.to_lgraph()
    LG.compute_wr(len(H)-1)
    log.info("Done.")

    cdag.count(LG)
