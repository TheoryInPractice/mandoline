#!/usr/bin/env python3

from graph import Graph, load_graph
from pattern import PatternBuilder, Pattern
import math, random

if __name__ == '__main__':
	H = Graph()
	H.add_edge(0,1)
	H.add_edge(1,2)
	H.add_edge(2,3)
	H.add_edge(3,4)

	for pattern,indexmap in H.enum_patterns():
		print(pattern)
		for piece in pattern.decompose():
			print("   ", piece)