#! /usr/bin/python3

import sys
import re
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

    # need to do lengths for diffsel, stupidly
    visited = []
    if args.nCond == 3:
        label1 = False
        label2 = False
        for n in curroot.iternodes(order="preorder"):  # start at root
            if n.istip:
                if n.label.endswith("#1"):
                    n.length = 1.0
                    n.label = n.label[:-2]
                    visited.append(n)
                elif n.label.endswith("#2"):
                    n.length = 2.0
                    n.label = n.label[:-2]
                    visited.append(n)
            if n.label == "#1":
                label1 = True
                n.label = ""
                n.length = 1.0
                visited.append(n)
                for c in n.iternodes(order="preorder"):
                    c.length = 1.0
                    visited.append(c)
            elif n.label == "#2":
                n.label = ""
                label2 = True
                visited.append(n)
                for c in n.iternodes(order="preorder"):
                    c.length = 2.0
                    visited.append(c)
            else:
                if n not in visited:
                    n.length = 0.0
        if not label2:
            sys.stderr.write("expected a class 2 label \
                             but only 1 found. Check input.\n")
    else:
        label1 = False
        for n in curroot.iternodes(order="preorder"):  # start at root
            if n.istip:
                if n.label.endswith("#1"):
                    n.label = n.label[:-2]
                    n.length = 1.0
                    visited.append(n)
            if n.label == "#1":
                label1 = True
                n.label = ""
                n.length = 1.0
                visited.append(n)
                for c in n.iternodes(order="preorder"):
                    c.length = 1.0
                    visited.append(c)
            if n.label == "#2":  # in case formatted for both
                n.label = ""
            else:
                if n not in visited:
                    n.length = 0.0

    newNwkString = curroot.get_newick_repr(showbl=True)
    # now we need to strip the branch lengths to integers
    p = re.compile(":[0-9].[0-9]*")

    def repl(m):
        n = m.group(0).split(".")[0]
        return n

    newNwkString = p.sub(repl, newNwkString)

    print(newNwkString+";")
