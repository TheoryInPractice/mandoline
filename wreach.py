#!/usr/bin/env python3

import sys
from graph import Graph, load_graph
from pattern import Pattern

G = load_graph(sys.argv[1])

LG = G.to_lgraph()
LG.compute_wr(int(sys.argv[2]))

print(';'.join(["wr{}".format(d) for d in range(LG.depth())]))
for i in LG:
	print(	';'.join([str(LG.wr[d][i]) for d in range(LG.depth())]))

