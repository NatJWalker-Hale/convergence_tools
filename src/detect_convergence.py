#! /usr/bin/python3

import os
import sys
import argparse
import subprocess
import tree_reader
from copy import deepcopy
from parse_fasta import parse_fasta
from label_for_diffsel import label_for_diffsel
from label_for_tdg09 import label_for_tdg09
from prep_for_msd import prep_for_msd


def seqDict_to_phylip(seqDict, outPath):
    with open(outPath, "w") as outf:
        for k, v in seqDict.items():
            outf.write(k + "\t")
            outf.write(v + "\n")


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
    return scenarios, nCond


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
                        help="number of threads to use. Impacts how jobs \
                        are parallelised. Defaults mean all 5 methods will \
                        be run in parallel on 1 thread each, while nt < 6 \
                        means jobs will be run sequentially, with nt threads \
                        each, where possible.")
    args = parser.parse_args()
    # print(args.methods)

    # aa_seqs = dict([x for x in parse_fasta(args.AA_aln)])
    cds_seqs = dict([x for x in parse_fasta(args.CDS_aln)])

    with open(args.tree, "r") as inf:
        nwkString = inf.readline().strip()
        curroot = tree_reader.read_tree_string(nwkString)

    curroot.number_tree()

    if curroot.is_rooted():
        # create unrooted tree for diffsel
        if "diffsel" in args.methods:
            curroot_unroot = deepcopy(curroot)
            curroot_unroot.unroot()
            curroot_unroot.number_tree()
    else:
        sys.stdout.write("a rooted tree is required!\n")
        sys.exit()

    scenarios, conds = parse_scenario(args.scenario)

    cmds = {}
    if "diffsel" in args.methods:
        # make files
        os.mkdir("diffsel")
        seqDict_to_phylip(cds_seqs, "diffsel/diffsel_aln.phy")
        label_for_diffsel(curroot_unroot, scenarios, conds,
                          "diffsel/diffsel_tree.nwk")
        cmd1 = ["diffsel",
                "-d",
                "diffsel/diffsel_aln.phy",
                "-t",
                "diffsel/diffsel_tree.nwk",
                "-ncond",
                str(conds),
                "-x",
                "1 3000",
                "diffsel"+str(conds)+"cond/run1"]
        cmd2 = ["diffsel",
                "-d",
                "diffsel/diffsel_aln.phy",
                "-t",
                "diffsel/diffsel_tree.nwk",
                "-ncond",
                conds,
                "-x",
                "1 3000",
                "diffsel/"+str(conds)+"cond/run2"]
    cmds["diffsel1"] = cmd1
    cmds["diffsel2"] = cmd2

    print(cmds["diffsel1"])
    subprocess.run(cmds["diffsel1"])
