from .build_counting_dag import decompose
from .count import count

import sys
import argparse

from .helpers import CheckExt

def main():
    parser = argparse.ArgumentParser(description='Graph counting based on weak colouring orders')
    subparsers = parser.add_subparsers(dest="mode", help='Program mode')
    subparsers.required = True

    parser.add_argument('--debug', action='store_true')
    parser.add_argument('--quiet', action='store_true' )

    decomp_parser = subparsers.add_parser("decompose", description="Decompose a pattern graph H and produce a counting DAG as output.")
    decomp_parser.add_argument('pattern', help='Pattern graph file')
    decomp_parser.add_argument('--output', '-o', help='Output file for counting DAG', action=CheckExt({'dag'}))
    decomp_parser.set_defaults(func=decompose)

    count_parser = subparsers.add_parser("count", description="Counts number of times a pattern graph appears induced in a host graph")

    count_parser.add_argument('cdag', help='Counting DAG file of the pattern graph', action=CheckExt({'dag'}))
    count_parser.add_argument('host', help='The host graph')
    count_parser.add_argument('--no-reduction', action='store_true' )
    count_parser.add_argument('--natural-order', action='store_true', help='Use order as given by node ids' )
    count_parser.set_defaults(func=count)

    args = parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    sys.exit(main())