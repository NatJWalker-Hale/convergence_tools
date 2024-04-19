#! /usr/bin/python3

import sys
import argparse
import tree_reader


def parse_scenario(path):
    scenarios = {}
    scenarioStr = {}
    with open(path, "r") as inf:
        nCond = 1
        for s in inf:
            if s != "\n":
                scenarios[nCond] = []
                scenarioStr[nCond] = s.strip()
                for i in s.strip().split("/"):
                    scenarios[nCond] += [int(x) for x in i.split(",")]
                nCond += 1
    return scenarios, nCond, scenarioStr


def label_for_nhx(root, scenario, scenStr, pelican=False):
    root.number_tree()
    #print(root.get_newick_repr(shownum=True)+";")
    if pelican:
        for n in root.iternodes(order="preorder"):
            if n.number in scenario:
                n.note += "&&NHX:Trait=foreground"
            else:
                n.note += "&&NHX:Trait=background"
    else:
        for n in root.iternodes(order="preorder"):
            if n.number in scenario:
                n.note += "&&NHX:Condition=1:"
            else:
                n.note += "&&NHX:Condition=0:"
        trNodes = []
        for i in scenStr.split("/"):
            trNodes.append(int(i.split(",")[0]))
        for n in root.iternodes(order="preorder"):
            if n.number in trNodes:
                n.note += "Transition=1"
            else:
                n.note += "Transition=0"
    


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("tree", help="newick tree to tag")
    parser.add_argument("scenario", help="PCOC-formatted scenario string \
                        in text file")
    parser.add_argument("-p", "--pelican", help="label for PELICAN", action="store_true")
    args = parser.parse_args()

    with open(args.tree, "r") as inf:
        nwkString = inf.readline().strip()
        curroot = tree_reader.read_tree_string(nwkString)

    scenarios, _, scenStr = parse_scenario(args.scenario)

    if args.pelican:
        label_for_nhx(curroot, scenarios[1], scenStr[1], pelican=True)
    else:
        label_for_nhx(curroot, scenarios[1], scenStr[1])

    print(curroot.get_newick_repr_ete(showbl=True, shownum=False) + ";")
