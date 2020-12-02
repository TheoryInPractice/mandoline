#!/usr/bin/env python3


def stats(file):
    total = 0
    leaves = 0
    edges = 0
    dist = None

    with open(file, 'r') as f:
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
        for l in f:
            l = l.rstrip()
            if len(l) == 0:
                continue
            edges += len(l.split()) - 1
    return total, leaves, edges, dist


def main():
    graphs = ['P{}'.format(i) for i in range(3,6)]
    graphs += ['C{}'.format(i) for i in range(4,7)]
    graphs += ['S{}'.format(i) for i in range(3,6)]
    graphs += ['W{}'.format(i) for i in range(3,7)]
    graphs += ['K{},{}'.format(i,i) for i in range(3,6)]

    # for H in graphs:
    #     name = H[0] + "_{" + H[1:] + "}"
    #     total, leaves, edges, dist = stats('example-graphs/{}.dag'.format(H))
    #     print("${}$ & {} ({}) & {} & {} \\\\".format(name, total, leaves, edges, dist))

    n = 4
    for H in ['diamond', 'paw']:
        total, leaves, edges, dist = stats('example-graphs/{}.dag'.format(H))
        print("{} & {} & {} ({}) & {} & {} \\\\".format(n, H, total, leaves, edges, dist))
        n = ''

if __name__ == "__main__":
    main()
