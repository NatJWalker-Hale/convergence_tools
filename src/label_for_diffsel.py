#! /usr/bin/python3

import sys
import argparse
import tree_reader


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("treeFile", help="tree file in newick format with \
                        matching taxa")
    parser.add_argument("-n", "--nCond", type=int, default=2,
                        choices=[2, 3], help="number of conditions. If n=2 \
                        [default], nodes including and descending the \
                        foreground branch will be labelled 1, and all \
                        others 0. Else, the script will expect a #1 and #2 \
                        label in the tree. Nodes including and descending \
                        #1 will be treated as above, and same for #2, while \
                        all others will be labelled 0")
    args = parser.parse_args()

    with open(args.treeFile, "r") as inf:
        nwkString = inf.readline().strip()
        curroot = tree_reader.read_tree_string(nwkString)

    # first, strip any other node labels present
    for n in curroot.iternodes(order="preorder"):
        if n.istip:
            continue
        if "#" not in n.label:
            n.label = ""

    if args.nCond == 3:
        label1 = False
        label2 = False
        for n in curroot.iternodes(order="preorder"):
            if n.label == "#1":
                if not label1:
                    label1 = True
                    n.label = ":1"
                    for c in n.children:
                        if c.istip:
                            c.label += ":1"
                        else:
                            c.label = ":1"
            elif n.label == "#2":
                if not label2:
                    label2 = True
                    n.label = ":2"
                    for c in n.children:
                        if c.istip:
                            c.label += ":2"
                        else:
                            c.label = ":2"
            else:
                if n.istip:
                    n.label += ":0"
                else:
                    n.label = ":0"
        if not label2:
            sys.stderr.write("expected both #1 and #2 in tree, but only \
                             #1 found")
            sys.exit()

    print(curroot.get_newick_repr())
