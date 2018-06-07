#!/usr/bin/env python3

from graph import Graph, load_graph
from collections import Counter
from pattern import PatternBuilder, Pattern
import math, random

if __name__ == '__main__':
	H = Graph()
	H.add_edge(0,1)
	H.add_edge(1,2)
	H.add_edge(2,3)
	H.add_edge(3,4)
	H.add_edge(4,5)
	H.add_edge(5,6)
	H.add_edge(7,8)

	count_patterns = 0
	count_pieces = 0

	all_pieces = set()
	secondary_pieces = set()

	for pattern,indexmap in H.enum_patterns():
		# print(pattern, hash(pattern))
		count_patterns += 1

		patterns = pattern.decompose()
		for piece in patterns:
			# print("   ", piece, hash(piece))
			all_pieces.add(piece)
			count_pieces += 1

		for piece in patterns[:-1]:
			secondary_pieces.add(piece)

	print("")
	print("{} patterns".format(count_patterns))
	print("{} pieces".format(count_pieces))

	k = len(all_pieces)
	print("  of which {} ({:.1f}%) are unique".format(k,k/count_pieces * 100))

	k = len(secondary_pieces) + count_patterns
	print("  of which {} ({:.1f}%) are primary or unique secondary".format(k,k/count_pieces * 100))