#!/usr/bin/env python3

from graph import Graph, DiGraph, DAGError, load_graph
from datastructures import Bimap
from pattern import PatternBuilder, Pattern

import argparse
import itertools
import sys

from collections import defaultdict, Counter
from sortedcontainers import SortedSet
from itertools import permutations, product, combinations, chain

from helpers import CheckExt, suborder
from tree_decompose import TD, short_str

import logging

log = logging.getLogger("mandoline")

class Recorder:
    def __init__(self):
        self.pieces = set()
        self.decomps = set()
        self.base_decomps = set()

        self.product_edges = dict()
        self.product_edges_count = 0
        self.subtract_edges = defaultdict(set)
        self.subtract_edges_count = 0

    def add_base_decomp(self, td):
        if td in self.base_decomps:
            return True
        log.info("Found new base decomp %s", td.td_string())
        self.base_decomps.add(td)
        return False

    def count_linear(self, td):
        assert(td.is_linear())
        if td not in self.pieces:
            log.info("Found new piece %s", td.td_string())
        self.pieces.add(td)

    def count_recursive(self, td):
        assert(not td.is_linear())
        if td in self.decomps:
            return True

        self.decomps.add(td)
        log.info("Found new decomp %s", td.td_string())
        return False

    def count_product(self, td_left, td_right, td_result):
        if td_result in self.product_edges:
            assert self.product_edges[td_result] == (td_left, td_right)
        else:
            self.product_edges_count += 1
        self.product_edges[td_result] = (td_left, td_right)

    def count_subtract(self, td_count, td_subtract, multi):
        if td_subtract not in self.subtract_edges[td_count]:
            self.subtract_edges_count += 1
        self.subtract_edges[td_count].add((td_subtract, multi))

    def report(self):
        log.info("Recorded %d linear pieces, %d decompositions of which %d are the basis.", len(self.pieces), len(self.decomps), len(self.base_decomps))
        log.info("  We have %d product-count and %d subtract-count edges", self.product_edges_count, self.subtract_edges_count)

    def output(self, filename):
        log.info("Writing counting dag to %s", filename)

        # Construct dependency DAG for intermediate decompositions
        # in order to output decompositions in a topological order.
        # We know that 'base_decomps' are sources and 'pieces' are sinks,
        # hence we only need to compute an order for the intermediate
        # decompositions in 'decomps'.
        dag = DiGraph()
        for td in self.decomps:
            if td in self.base_decomps:
                continue
            dag.add_node(td)
        nodes = set(dag)
        for td in nodes:
            left, right =  self.product_edges[td]
            out_neighbours = set([left, right])
            out_neighbours.update([t for (t,_) in self.subtract_edges[td]])
            out_neighbours = out_neighbours & nodes
            for other in out_neighbours:
                dag.add_arc(td, other)
        dag.remove_loops()

        _, imap = dag.normalize()
        decomp_order = imap.order() # Iterator!

        index = dict()
        curr_index = 0
        with open(filename, 'w') as f:
            f.write('* Base\n')
            for td in self.base_decomps:
                index[td] = curr_index
                f.write('{} {}\n'.format(curr_index, td.td_string()))
                curr_index += 1
            f.write('* Composite\n')
            for td in decomp_order:
                if td in index:
                    assert td in self.base_decomps
                    continue
                index[td] = curr_index
                f.write('{} {}\n'.format(curr_index, td.td_string()))
                curr_index += 1
            f.write('* Linear\n')
            for td in self.pieces:
                if td in index:
                    assert td in self.base_decomps
                    continue
                index[td] = curr_index
                f.write('{} {}\n'.format(curr_index, td.td_string()))
                curr_index += 1
            f.write('* Edges\n')

            edges_rows = []
            for td, (left, right) in self.product_edges.items():
                assert td in index
                assert left in index
                assert right in index
                if index[left] > index[right]:
                    left, right = right, left # For consistency

                subtr_sorted = sorted([(index[td_sub],multi) for td_sub,multi in self.subtract_edges[td]])
                subtr_tokens = [str(x)+"|"+str(multi) for x, multi in subtr_sorted]
                edges_rows.append( (index[td], index[left], index[right], ' '.join(subtr_tokens)) )

            for e in sorted(edges_rows):
                f.write('{} {} {} {}\n'.format(*e))
        pass

def powerset(iterable):
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(0, len(s)+1))

def powerset_nonempty(iterable):
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(1, len(s)+1))

def simulate_count(R, H, td):
    already_known = R.add_base_decomp(td)
    if already_known:
        return
    _simulate_count_rec(R, H, td, 0)

def _simulate_count_rec(R, H, td, depth):
    prefix = " "*(4*depth)
    log.debug("%sWe want to count %s (%s) for adhesion size %s", prefix, td.td_string(), str(td), len(td._sep))

    split_depth = td.adhesion_size()
    splits = list(td.split())
    if len(splits) == 1:
        log.debug("%sThis decomposition is linear, we simply count it!", prefix)
        R.count_linear(td)
        return

    already_known = R.count_recursive(td)
    if already_known:
        return

    # Merge splits from left to right
    merged = [splits[0]]
    for s in splits[1:]:
        merged.append(merged[-1].merge(s, split_depth))

    assert merged[-1] == td


    log.debug("%sThe decomposition branches at depth %d", prefix, split_depth)

    td_current = splits[0]
    log.debug("%sWe first count the leftmost piece %s", prefix, td_current)
    _simulate_count_rec(R, H, td_current, depth+1)

    log.debug("%sNow we fold-count with the remaining pieces.", prefix)

    for tdB, result, tdA in zip(splits[1:], merged[1:], merged):
        log.debug("%sThe next piece is %s and we first count it.", prefix, tdB)
        _simulate_count_rec(R, H, tdB, depth+1)

        log.debug("%sThe initial count of %s is the count of %s times the count of %s", prefix, result, tdA, tdB)

        nodesA, nodesB, _ = td_overlap(tdA, tdB)

        subisosAintoB = 0
        subisosBintoA = 0
        for (H1, tdH1, mapping) in enumerate_merges(tdA, tdB):
            # Note: The resulting merge is labelled with vertices from nodesB,
            #   e.g. if x \in nodesA is mapped onto y \in nodesB, the resulting
            #   node has the label y.

            # Compute remaining vertices of nodesA, nodesB between which
            # edges are still allowed (the merge makes some edges unavailable!).
            # For this, we remove all nodes _above_ nodes that participate in the
            # merge.
            closure = tdH1.upwards_closure(mapping.target()) | mapping.source() # mapping.source() are labels not found in tdH1
            nodesA1 = nodesA - closure
            nodesB1 = nodesB - closure

            # Enumerate additional cases with additional edges
            for (H2, tdH2, edges) in enumerate_edge_faults(H1, tdH1, nodesA1, nodesB1, depth):
                if len(edges) == 0 and len(mapping) == 0:
                    assert tdH2 == result, "{} != {}".format(tdH2.td_string(), result.td_string())
                    continue # No merge, no edge addition. This case is included in the iteration for convenience.

                _simulate_count_rec(R, H2, tdH2, depth+1)
                coeff = compute_coefficient(result, nodesA, nodesB, tdH2)
                R.count_subtract(result, tdH2, coeff)

        R.count_product(tdA, tdB, result)

        # Make sure the resulting decomposition is noted
        R.count_recursive(result)

def td_overlap(decompA, decompB):
    """
        Computes three sets of nodes: nodes exclusive to decompA,
        nodes exclusive to decompB and all nodes of the joint
        decomposition (including their joint root-path).
    """
    rootPath = set(decompA._sep) & set(decompB._sep)
    nodesA = decompA.nodes()
    nodesB = decompB.nodes()
    common_nodes = nodesA & nodesB
    nodesA -= common_nodes
    nodesB -= common_nodes
    nodesAll = nodesA | nodesB | rootPath
    return (nodesA, nodesB, nodesAll)

def enumerate_merges(decompA, decompB):
    """
        Enumerates all graphs and td decompositions that
        contain decompA, decompB as subgraph/subdecompositions
        while overlapping in some subset of vertices.
    """
    rootPathSet = set(decompA._sep) & set(decompB._sep)
    rootPath = decompA._sep[:len(rootPathSet)]
    rootPathEdges = list(zip(rootPath[:-1],rootPath[1:]))
    assert len(rootPath) > 0 # Paranoia: decompositions must have common root-path.
    assert rootPath == decompB._sep[:len(rootPathSet)]

    nodesA, nodesB, nodesAll = td_overlap(decompA, decompB)

    decompJoint = decompA.merge(decompB, len(rootPath))
    # print("{} + {} = {}".format(decompA, decompB, decompJoint))
    dagJoint = decompJoint.to_ditree()
    graphJoint = decompJoint.to_graph()

    # print("TD Tree",dagJoint)
    # print("Graph", graphJoint)

    # Make a list of which nodes in nodesA can potentially be
    # mapped onto nodes in nodesB. The first constraint enforced here
    # is that the must have the same in-neighbourhood w.r.t the
    # (joint) root-path of the decomposition.
    candidates = {}
    for u in nodesA:
        candidates[u] = set()
        rp_neighboursA = decompA.in_neighbours(u) & rootPathSet
        for v in nodesB:
            rp_neighboursB = decompB.in_neighbours(v) & rootPathSet
            if rp_neighboursA == rp_neighboursB:
                candidates[u].add(v)

    # print(candidates)

    # Now try all subsets of nodesA, including the empty set
    seen_tdM = set()
    for sourceA in powerset(nodesA):
        candidate_sets = [candidates[x] for x in sourceA]
        subgraphA = graphJoint.subgraph(sourceA)
        for targetB in product(*candidate_sets):
            if len(set(targetB)) != len(targetB):
                continue # Mapping must be onto

            mapping = Bimap()
            mapping.put_all(zip(sourceA, targetB))

            # print("  Mapping {} -> {} = {}".format(sourceA, targetB, mapping))

            # Check that subgraphs induced by 'sourceA' and 'sourceB' are the same
            # under the mapping, e.g. that the mapping is a isomorphism.
            subgraphB = graphJoint.subgraph(targetB)

            if subgraphA.relabel(mapping) != subgraphB:
                continue # Not an isomorphism

            # print("  Isomorphism! {} == {}".format(subgraphA.relabel(mapping), subgraphB))

            dagMerged = dagJoint.copy()
            dagMerged.merge_pairs(mapping.items())
            # print("  Dag contraction {} --{}--> {}".format(dagJoint, mapping, dagMerged))

            # Check whether resulting decomposition graph is a DAG, otherwise
            # this mapping is incompatible with the orderings associated with decompA, decompB.
            if not dagMerged.is_acyclic():
                # print("  > not a DAG")
                continue

            # This is going somewhere, so we can finally construct the
            # resulting graph
            graphMerged = graphJoint.copy()
            graphMerged.merge_pairs(mapping.items())

            # print("  Graph contraction {} --{}--> {}".format(graphJoint, mapping, graphMerged))

            # Decompose resulting graph according to DAG embeddings
            # print("  TD decompositions of {}:".format(graphMerged))
            for o in dagMerged.embeddings():
                # We decompose 'graphMerged' with the additional constraint that
                # the root-path must stay a prefix of the resulting decomposition
                # (by using virtual edges along that path)
                tdMerged = TD.decompose(graphMerged, o, rootPathEdges)

                # Some embeddings will produce the same td decomposition, but in this
                # inner loop we consider every decomposition only once.
                if tdMerged in seen_tdM:
                    continue
                yield graphMerged, tdMerged, mapping

def enumerate_edge_faults(H, tdH, nodesA, nodesB, depth):
    """
        Enumerates all td decompositions that can be constructed from tdH by
        adding edges between the sets nodesA, nodesB. We enforce that the root-path
        of tdH must be a prefix of the resulting decompositions' root-path.
    """
    prefix = " "*(4*depth)
    log.debug("%sTo account for non-induced instance, edges between %s and %s need to be considered", prefix, nodesA, nodesB )
    potential_edges = list(product(nodesA, nodesB))

    # We add virtual root-path edges in order to preserve the root-path
    # in the resulting decompositions.
    rootPath = tuple(tdH._sep)
    rootPathEdges = list(zip(rootPath[:-1],rootPath[1:]))

    # Return original decomposition without edges added
    yield H, tdH, tuple()

    if len(potential_edges) == 0:
        return

    log.debug("%sWe subtract the results of the following counts:", prefix)
    for edges in powerset_nonempty(potential_edges):
        seen_decomp = set()
        for o in tdH.suborders(H):
            assert len(o) > 0
            HH = H.copy()
            HH.add_edges(edges)
            log.debug("%sDecompositing %s along order %s", prefix,list(HH.edges()), o)
            tdHH = TD.decompose(HH, o, rootPathEdges)
            if tdHH in seen_decomp:
                continue
            assert tuple(tdHH._sep[:len(rootPath)]) == tuple(rootPath)
            seen_decomp.add(tdHH)
            yield HH, tdHH, edges

def count_automorphisms(td, nodes):
    """
        Counts all automorphisms of the subgraph induced
        by (root-path + nodes) in which the root-path is mapped by
        the identity.
    """
    return sum([1 for _ in compute_automorphisms(td, nodes)])

def compute_automorphisms(td, nodes):
    """
        Computes all depth-preserving automorphisms of the subgraph induced
        by (root-path + nodes). Since the depths of nodes are preserved, the
        root-path will always be mapped to itself and is _not_ returned
    """
    rootPath = tuple(td._sep)
    rootPathSet = set(rootPath)
    subgraph = td.to_graph().subgraph(nodes)

    # Partition nodes according to a) their depth and b)
    # their root-path neighbourhood
    partition = defaultdict(list)
    for s in nodes:
        rpNeighS = td.in_neighbours(s) & rootPathSet
        N = suborder(rootPath, rpNeighS)
        d = td.depth_of(s)
        partition[(d,N)].append(s)
    partition = list(partition.values())

    # Combine all permutations of each block and check whether the
    # resulting mapping is an automorphism
    for perms in product(*[permutations(p) for p in partition]):
        mapping = dict()
        for orig, perm in zip(partition, perms):
            mapping.update(zip(orig, perm))
        if subgraph.relabel(mapping) == subgraph:
            yield mapping


def compute_coefficient(td, nodesA, nodesB, defect):
    """
        Returns how many mappings from 'td' to 'defect' there
        are such that the two pieces induced by (root-path + nodesA) and (root-path + nodesB)
        are individually preserved (but if td != defect they will necessarily either
        intersect or be connected by an unwanted edge). Note that we only need
        to consider mappings that are surjective, e.g. all nodes of 'defect' must be hit.
    """
    rootPath = tuple(td._sep)
    rootPathSet = set(rootPath)
    nodesA, nodesB = list(nodesA), list(nodesB) # We need a consistent iteration order
    defectNodes = set(defect.nodes()) - rootPathSet
    graph = td.to_graph()
    graphDefect = defect.to_graph()
    subgraphA = graph.subgraph(set(nodesA))
    subgraphB = graph.subgraph(set(nodesB))

    assert rootPath == tuple(defect._sep[:len(rootPath)]) # Ensure decompositions agree on labelling of root-path

    # Make list of potential mapping candidates; the constraint here
    # is that source and target vertex must agree on the root-path neighbourhood.
    candidates = defaultdict(set)
    for s in chain(nodesA, nodesB):
        rpNeighS = td.in_neighbours(s) & rootPathSet
        for t in (defect.nodes() - rootPathSet):
            rpNeighT = defect.in_neighbours(t) & rootPathSet
            if rpNeighS == rpNeighT:
                candidates[s].add(t)
    # print("Candidates:", dict(candidates))
    # print("Need to cover nodes", defectNodes)

    candidateSetsA = [candidates[x] for x in nodesA]
    candidateSetsB = [candidates[x] for x in nodesB]
    count = 0
    for choicesA in product(*candidateSetsA):
        if len(set(choicesA)) != len(nodesA):
            continue # Not a bijection

        mappingA = dict(zip(nodesA, choicesA))
        mappingArev = dict(zip(choicesA, nodesA))

        # Ensure that 'choicesA' induce a subgraph in 'graphDefect' that is
        # isomorphic the subgraph induced by 'nodesA' in 'graph'. Note that edges
        # towards the root-path are not tested here, those are taken care of by
        # our selection of candidates.
        if subgraphA.relabel(mappingA) != graphDefect.subgraph(choicesA):
            continue

        # Check whether mapping is order-compatible; meaning that valid orderings
        # of the target nodes according the 'defect' decomposition must all be
        # compatible with the original decomposition.
        # TODO: This could probably done more efficiently by considering the induced
        #       posets.
        orderCompatible = True
        for o in defect.suborders(choicesA):
            oo = [mappingArev[x] for x in o]
            if not td.compatible_with(oo):
                orderCompatible = False
                break

        if not orderCompatible:
            continue

        for choicesB in product(*candidateSetsB):
            if len(set(choicesB)) != len(nodesB):
                continue # Not a bijection

            mappingB = dict(zip(nodesB, choicesB))
            mappingBrev = dict(zip(choicesB, nodesB))
            targets = set(choicesA) | set(choicesB)
            if targets != defectNodes:
                continue

            # Similar to above, but for our mapping of 'choicesB'.
            if subgraphB.relabel(mappingB) != graphDefect.subgraph(choicesB):
                continue

            orderCompatible = True
            for o in defect.suborders(choicesB):
                oo = [mappingBrev[x] for x in o]
                if not td.compatible_with(oo):
                    orderCompatible = False
                    break

            if not orderCompatible:
                continue

            # Whatever remains is a valid mapping!
            count += 1

    autoA = count_automorphisms(td, nodesA)
    autoB = count_automorphisms(td, nodesB)
    if count % (autoA * autoB) != 0:
        print(td.td_string())
        print(td)
        print(nodesA)
        print(nodesB)
        print(defect.td_string())
    assert count % (autoA * autoB) == 0, "{} not divisible by {} = {} x {}".format(count, autoA*autoB, autoA, autoB)
    count //= autoA * autoB

    return count

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Counts H in G')

    parser.add_argument('H', help='Pattern graph H')
    # parser.add_argument('G', help='Host graph G')
    parser.add_argument('--debug', action='store_true')
    parser.add_argument('--quiet', action='store_true' )
    parser.add_argument('--output', '-o', help='Output file for counting DAG', action=CheckExt({'dag'}))

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

    if args.output:
        R.output(args.output)
