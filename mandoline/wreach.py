#!/usr/bin/env python3

import sys
from graph import Graph, load_graph
from pattern import Pattern
import argparse


parser = argparse.ArgumentParser(
			description='Computes wreach sets of a graph and prints them to the terminal.')

parser.add_argument('G', help='Host graph G')
parser.add_argument('wreach', type=int,
						help='Up which depth the wreach sets are computed')
parser.add_argument('--natural-order', action='store_true', help='Use order as given by node ids' )

args = parser.parse_args()

G = load_graph(sys.argv[1])

if args.natural_order:
    LG, _ = G.to_lgraph(sorted(G))
else:
    LG, _ = G.to_lgraph()
LG.compute_wr(args.wreach+1)

print(';'.join(["wr{}".format(d) for d in range(LG.depth())]))
for i in LG:
	print(	';'.join([str(list(LG.wr[d][i])) for d in range(LG.depth())]))
