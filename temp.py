#!/usr/bin/env python3

import sys, os

total = 0
leaves = 0
dist = None

with open(sys.argv[1], 'r') as f:
    for _ in range(3):
        next(f)
    l = next(f)
    dist = int(l.split()[1])
    for l in f:
        l = l.rstrip()

        if l == '* Edges':
            break
        if l[0] == '*':
            continue
        total += 1
        if '{' not in l:
            leaves += 1
name = os.path.basename(sys.argv[1])
name = name.split('.')[0]
print("{} & {} & {} & {} \\\\".format(name, total, leaves, dist))

