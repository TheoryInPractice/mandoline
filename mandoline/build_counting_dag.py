#!/usr/bin/env python3

import argparse
import sys

from collections import defaultdict
from sortedcontainers import SortedSet
from itertools import permutations, product, chain
from operator import itemgetter

import logging

from .graph import DiGraph, load_graph
from .datastructures import Bimap

from .helpers import CheckExt, suborder, powerset_nonempty, powerset
from .helpers import short_str
from .tree_decompose import TD


def join_dicts(dictA, dictB, f):
    """
        Joins to dictionaries. The value of elements sharing the same key
        is determined by the function f, that is, if `dictA` and `dictB` both
        contain the key `k` then the resulting dictionary will contain the
        element `f(dictA[k], dictB[k])` for `k`.

        Values with non-shared keys are simply copied over.
    """
    res = {}
    res.update(dictA)
    res.update(dictB)
    for k in (dictA.keys() &  dictB.keys()):
        res[k] = f(dictA[k], dictB[k])
    return res


class TDIndex:
    def __init__(self):
        self.td_to_ix = {} 
        self.ix_to_td = []
        self.curr_index = 0

    def append(self, td):
        assert td not in self.td_to_ix
        self.td_to_ix[td] = self.curr_index
        self.ix_to_td.append(td)
        assert len(self.ix_to_td) == self.curr_index+1
        self.curr_index += 1

    def __contains__(self, td):
        return td in self.td_to_ix

    def __getitem__(self, key):
        if isinstance(key, int):
            return self.ix_to_td[key]
        else:
            return self.td_to_ix[key]

    def __len__(self):
        return len(self.ix_to_td)


log = logging.getLogger("mandoline")

class Recorder:
    def __init__(self, graph):
        self.graph = graph
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

    def count_product(self, td_left, td_right, td_result, coeff):
        if td_result in self.product_edges:
            assert self.product_edges[td_result] == (td_left, td_right, coeff)
        else:
            self.product_edges_count += 1
        self.product_edges[td_result] = (td_left, td_right, coeff)

    def count_subtract(self, td_count, td_subtract, multi):
        if td_subtract not in self.subtract_edges[td_count]:
            self.subtract_edges_count += 1
        self.subtract_edges[td_count].add((td_subtract, multi))

    def report(self):
        log.info("Recorded %d linear pieces, %d decompositions of which %d are the basis.", len(self.pieces), len(self.decomps), len(self.base_decomps))
        log.info("  We have %d product-count and %d subtract-count edges", self.product_edges_count, self.subtract_edges_count)

    def compute_decomp_order(self):
        """
            Construct dependency DAG for intermediate decompositions
            in order to output decompositions in a topological order.
            We know that 'base_decomps' are sources and 'pieces' are sinks,
            hence we only need to compute an order for the intermediate
            decompositions in 'decomps'.
        """
        dag = DiGraph()
        for td in self.decomps:
            if td in self.base_decomps:
                continue
            dag.add_node(td)
        nodes = set(dag)
        for td in nodes:
            left, right, _ =  self.product_edges[td]
            out_neighbours = set([left, right])
            out_neighbours.update([t for (t,_) in self.subtract_edges[td]])
            out_neighbours = out_neighbours & nodes
            for other in out_neighbours:
                dag.add_arc(td, other)
        dag.remove_loops()

        _, imap = dag.normalize()
        return imap.order() # Iterator

    def compute_index(self):
        # Compute order in which to compute intermediate decompositions,
        # e.g. a topological ordering of the dependency-dag.
        decomp_order = self.compute_decomp_order()

        # Build index for td decompositions. There is a slight complication here
        # since a decomposition might appear both as a 'base' decomposition and
        # as a 'decomp' or a 'piece'. Hence the loop unrolling.
        td_index = TDIndex()
        for td in self.base_decomps:
            td_index.append(td)

        index_base = td_index.curr_index
        for td in decomp_order:
            if td in td_index:
                assert td in self.base_decomps
                continue
            td_index.append(td)

        index_decomp = td_index.curr_index
        for td in self.pieces:
            if td in td_index:
                assert td in self.base_decomps
                continue
            td_index.append(td)
        index_pieces = td_index.curr_index

        return td_index, (index_base, index_decomp, index_pieces)

    def compute_nroot_hints(self, td_index):
        """
            An /nroot/ is a node in a linear decomposition where it looks like the 
            connectvitiy is broken, i.e. if we remove all preceding vertices up to the
            nroot we find that the remaining graph has several connected components. 
            However, we _know_ that we are only interested in linear occurences of this
            pattern where these vertices are connected. This connectivity happens
            in another linear piece. 

            It is therefore enough to search for the nroot vertex within the wreach
            set of some subsequent vertex, which we need to idenfify. To do that, we consider
            all "product parents", that is, all parent decompositions T which are decomposed
            into linear pieces L_1 + L_2. For example, let's say that vertex u is an
            nroot in L_1 and u has distance 5 from v (which comes after u in L_1) in
            T. Then we know that it is sufficient to search for u in W^5(v) and we call
            v a /hint/ vertex.

            Of course, our goal is to find the best hint vertices, meaning those that
            minimize the distance to u because this globally minimises the wreach number
            we have to compute.

            To simplify this algorithm, we want to find a single hint vertex for each nroot
            (instead of one hint vertex per every "product parent"). Therefore we collect all
            candidates across all product parents and then choose one that works in every
            situation.
        """

        product_parents = defaultdict(set)

        for parent, (left, right, _) in self.product_edges.items():
            iparent = td_index[parent]
            product_parents[td_index[left]].add(iparent)
            product_parents[td_index[right]].add(iparent)

        nroot_hints = defaultdict(list)


        def find_distances(i, nroots, depth=0):
            res = dict([(x,{}) for x in nroots])
            for p in product_parents[i]:
                td_parent = td_index[p]
                dists = td_parent.to_graph().all_distances()
                for r in nroots:
                    if r not in dists:
                        rec = find_distances(p, set([r]), depth+1)
                        res[r] = join_dicts(res[r], rec[r], max)
                    else:
                        # print(" "*(2*depth+1), r, dists[r])
                        res[r] = join_dicts(res[r], dists[r], max)
            return res

        for i in range(len(td_index)):
            td = td_index[i]

            if not td.is_linear():
                continue

            graph = td.to_graph()
            piece = td.to_piece(len(graph))
            nroots = SortedSet(piece.nroots) # These are indices, not vertices!
            nroots.remove(piece.depth()-1) # Remove roots of piece, it's alway an nroot

            if len(nroots) == 0:
                continue

            order = next(td.orders()) # Piece is linear, so there is only one order

            nroots = SortedSet([order[r] for r in nroots]) # Map indices to vertices

            # print(i, td.td_string(), td, nroots, order);
            nroot_dists = find_distances(i, nroots)
            for r in nroots:
                r_index = order.index(r)
                # Pair nroot index, hint index with distances
                cands = [(r_index, i, nroot_dists[r][x]) for (i,x) in enumerate(order)]
                cands = cands[r_index+1:] # Only allow vertices to the right of r
                nroot_hints[i].append(min(cands, key=itemgetter(2)))

        return nroot_hints

    def output(self, filename):
        log.info("Writing counting dag to %s", filename)

        # Compute index for td decompositions
        td_index, boundaries = self.compute_index()
        index_base, index_decomp, index_pieces = boundaries

        nroot_hints = self.compute_nroot_hints(td_index)

        # Determine maximum wreach needed
        max_wreach = 1
        for hints in nroot_hints.values():
            for (_,_,dist) in hints:
                max_wreach = max(max_wreach, dist)

        # Write to file
        with open(filename, 'w') as f:
            # Write preamble
            f.write('* Graph\n')
            node_str = ' '.join(map(str, self.graph.nodes))
            f.write('nodes ' + node_str + '\n')
            edge_str = ' '.join(map(lambda x: '{}|{}'.format(x[0],x[1]), self.graph.edges()))
            f.write('edges ' + edge_str + '\n')
            f.write('wreach {}\n'.format(max_wreach))

            # Write decompositions
            f.write('* Base\n')
            for i in range(index_base):
                assert(i not in nroot_hints)
                f.write('{} {}\n'.format(i, td_index[i].td_string()))
            f.write('* Composite\n')
            for i in range(index_base, index_decomp):
                assert(i not in nroot_hints)
                f.write('{} {}\n'.format(i, td_index[i].td_string()))
            f.write('* Linear\n')
            for i in range(index_decomp, index_pieces):
                if i in nroot_hints:
                    # nroot_hints[i] contains tuples of the form (nroot_index, hint_index, wreach_dist)
                    nroot_str = ' '.join(map(lambda e: '{}|{}|{}'.format(*e) , nroot_hints[i]))
                    f.write('{} {} {}\n'.format(i, td_index[i].td_string(), nroot_str))
                else:
                    f.write('{} {}\n'.format(i, td_index[i].td_string()))

            # Write 'edges' of counting-DAG
            f.write('* Edges\n')
            edges_rows = []
            for td, (td_left, td_right, auto_coeff) in self.product_edges.items():
                assert td in td_index
                assert td_left in td_index
                assert td_right in td_index

                left, right = td_index[td_left], td_index[td_right]
                if left > right:
                    left, right = right, left # For consistency
                    td_left, td_right = td_right, td_left

                subtract = dict([(td_index[td_sub],m) for td_sub,m in self.subtract_edges[td]])

                subtr_sorted = sorted(subtract.items())
                subtr_tokens = [str(x)+"|"+str(multi) for x, multi in subtr_sorted]
                edges_rows.append( (td_index[td], left, right, auto_coeff, ' '.join(subtr_tokens)) )

            for e in sorted(edges_rows):
                f.write('{} {}x{}|{} {}\n'.format(*e))


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

        nodesA, nodesB, _ = tdA.overlap_with(tdB)

        for (coeff, defect_graph, defect_td) in enumerate_defects(result, tdA, tdB, depth):
            _simulate_count_rec(R, defect_graph, defect_td, depth+1)
            R.count_subtract(result, defect_td, coeff)

        # Compute automorphism coefficient: assume we find the graph `result`,
        # how many pairs of embeddings of `tdA`, `tdB` are that that
        # cover `result`?
        coeff = result.count_relaxed_automorphisms(nodesA, nodesB)
        R.count_product(tdA, tdB, result, coeff)

        # Make sure the resulting decomposition is noted
        R.count_recursive(result)

def enumerate_defects(result, tdA, tdB, depth):
    nodesA, nodesB, _ = tdA.overlap_with(tdB)

    for (H1, tdH1, mapping) in enumerate_merges(tdA, tdB):
        # Note: The resulting merge is labelled with vertices from nodesB,
        #   e.g. if x \in nodesA is mapped onto y \in nodesB, the resulting
        #   node has the label y.

        # Instead of finding out which edge-pairs are still allowed, we
        # enumerate all possible edge-pairs and let compute_coefficient below
        # figure out whether tdA, tdB can still be found in the resulting graph.
        nodesA1 = nodesA - mapping.source()
        nodesB1 = nodesB

        # Enumerate additional cases with additional edges
        for (H2, tdH2, edges) in enumerate_edge_faults(H1, tdH1, nodesA1, nodesB1, depth):
            if len(edges) == 0 and len(mapping) == 0:
                assert tdH2 == result, "{} != {}".format(tdH2.td_string(), result.td_string())
                continue # No merge, no edge addition. This case is included in the iteration for convenience.

            coeff = result.count_relaxed_embeddings(nodesA, nodesB, tdH2)
            if coeff == 0:
                continue # Could not find a valid tdA, tdB mapping
            
            yield (coeff, H2, tdH2)

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

    nodesA, nodesB, _ = decompA.overlap_with(decompA)

    decompJoint = decompA.merge(decompB, len(rootPath))
    dagJoint = decompJoint.to_ditree()
    graphJoint = decompJoint.to_graph()

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

            # Check that subgraphs induced by 'sourceA' and 'sourceB' are the same
            # under the mapping, e.g. that the mapping is a isomorphism.
            subgraphB = graphJoint.subgraph(targetB)

            if subgraphA.relabel(mapping) != subgraphB:
                continue # Not an isomorphism

            dagMerged = dagJoint.copy()
            dagMerged.merge_pairs(mapping.items())

            # Check whether resulting decomposition graph is a DAG, otherwise
            # this mapping is incompatible with the orderings associated with decompA, decompB.
            if not dagMerged.is_acyclic():
                continue

            # This is going somewhere, so we can finally construct the
            # resulting graph
            graphMerged = graphJoint.copy()
            graphMerged.merge_pairs(mapping.items())

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
            log.debug("%sDecomposing %s along order %s", prefix,list(HH.edges()), o)
            tdHH = TD.decompose(HH, o, rootPathEdges)
            if tdHH in seen_decomp:
                continue
            assert tuple(tdHH._sep[:len(rootPath)]) == tuple(rootPath)
            seen_decomp.add(tdHH)
            yield HH, tdHH, edges

def decompose(args):
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
    H = load_graph(args.pattern)
    log.info("Loaded pattern graph with {} vertices and {} edges".format(len(H), H.num_edges()))
    log.info(H)

    seen = set()
    R = Recorder(H)
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
    else:
        log.info("No output was written. To write the cdag, use the --output argument.")
