#! /usr/bin/python3

import sys
import argparse
import subprocess
import tree_reader
from copy import deepcopy
from label_for_diffsel import label_for_diffsel
from label_for_tdg09 import label_for_tdg09
from prep_for_msd import prep_for_msd


def unroot_tree(path):
    # relies on pxrr from phyx - v hacky, figure out
    # how to do myself later
    cmd = ["pxrr", "-t", path, "-u"]
    res = subprocess.run(cmd, capture_output=True)
    return res.stdout


def parse_scenario(path):
    scenarios = {}
    with open(path, "r") as inf:
        nCond = 1
        for s in inf:
            if s != "\n":
                scenarios[nCond] = []
                for i in s.strip().split("/"):
                    scenarios[nCond] += [int(x) for x in i.split(",")]
                nCond += 1
    return scenarios


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("tree", help="tree on which to conduct convergence \
                        analyses. Should be rooted and have branch lengths \
                        optimised under an appropriate amino acid model \
                        (e.g. JTT, WAG, LG)")
    parser.add_argument("CDS_aln", help="CDS alignment of sequences matching \
                        tree")
    parser.add_argument("AA_aln", help="amino acid alignment of sequences \
                        matching tree")
    parser.add_argument("scenario", help="PCOC-formatted scenario string \
                        in text file")
    parser.add_argument("-m", "--methods", nargs="+", default=["diffsel",
                        "msd", "PCOC", "tdg09", "topo"],
                        help="space-separated detection methods to run, \
                        default all")
    parser.add_argument("-nt", "--num_threads", type=int, default=6,
                        help="number of threads to use. Has no impact on \
                        individual jobs, which will each run on 1 thread \
                        but will impact how jobs are parallelised. Defaults \
                        mean all 5 methods will be run in parallel")
    args = parser.parse_args()
    #print(args.methods)

    with open(args.tree, "r") as inf:
        nwkString = inf.readline().strip()
        curroot = tree_reader.read_tree_string(nwkString)

    curroot.number_tree()

    if curroot.is_rooted():
        # create unrooted tree for diffsel
        curroot_unroot = deepcopy(curroot)
        curroot_unroot.unroot()
        curroot_unroot.number_tree()
    else:
        sys.stdout.write("a rooted tree is required!")
        sys.exit()

    print(curroot.get_newick_repr(True) + ";")
    print(curroot_unroot.get_newick_repr(True) + ";")