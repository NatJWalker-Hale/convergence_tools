#!/usr/bin/env python3


import os
import sys
import glob
import argparse
import subprocess
import pandas as pd
import tree_reader
import calc_topological
from copy import deepcopy
from parse_fasta import parse_fasta
from label_for_diffsel import label_for_diffsel
from label_for_tdg09 import label_for_tdg09
from optimise_raxmlng import optimise, root_and_clean
from prep_for_msd import prep_for_msd


def seqDict_to_phylip(seqDict, outPath):
    with open(outPath, "w") as outf:
        outf.write(str(len(seqDict)) + "\t" +
                   str(len(list(seqDict.values())[0])) + "\n")
        for k, v in seqDict.items():
            outf.write(k + "\t")
            outf.write(v + "\n")


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
    parser.add_argument("-c0", "--cond0", help="background condition for \
                        tdg09. Default D2", type=str, default="D2")
    parser.add_argument("-c1", "--cond1", help="foreground condition for \
                        tdg09. Default D1", type=str, default="D1")
    parser.add_argument("-nt", "--num_threads", type=int, default=6,
                        help="number of threads to use. Impacts how jobs \
                        are parallelised. Defaults mean all 5 methods will \
                        be run in parallel on 1 thread each, while nt < 6 \
                        means jobs will be run sequentially, with nt threads \
                        each, where possible.")
    args = parser.parse_args()
    print(args.methods)
    # need to sanitise inputs at some point so that shell injections cannot
    # occur. Could do this by writing tmp aln and tree files so I know
    # exactly what they are

    TOP = os.getcwd()
    print(TOP)

    aa_seqs = dict([x for x in parse_fasta(args.AA_aln)])
    cds_seqs = dict([x for x in parse_fasta(args.CDS_aln)])

    with open(args.tree, "r") as inf:
        nwkString = inf.readline().strip()
        curroot = tree_reader.read_tree_string(nwkString)

    curroot.number_tree()

    if curroot.is_rooted():
        # determine og. Assumes child with fewest taxa is og
        ch1, ch2 = curroot.children
        ch1Lvs = ch1.lvsnms()
        ch2Lvs = ch2.lvsnms()
        if len(ch1Lvs) > len(ch2Lvs):
            og = ch2Lvs
        else:
            og = ch1Lvs
        print(og)
        # create unrooted tree for diffsel
        if "diffsel" in args.methods:
            curroot_unroot = deepcopy(curroot)
            curroot_unroot.unroot()
            curroot_unroot.number_tree()
            print(curroot_unroot.get_newick_repr()+";")
    else:
        sys.stderr.write("a rooted tree is required!\n")
        sys.exit()

    scenarios, conds, scenStr = parse_scenario(args.scenario)
    print(scenStr)

    cmds = {}

    if "diffsel" in args.methods:
        # make files
        try:
            os.mkdir("diffsel")
        except FileExistsError:
            sys.exit("A diffsel run already exists. "
                     "Delete to restart.")
        seqDict_to_phylip(cds_seqs, "diffsel/diffsel_aln.phy")
        label_for_diffsel(curroot_unroot, scenarios, conds,
                          "diffsel/diffsel_tree.nwk")
        cmd1 = ["diffsel",
                "-d",
                "diffsel_aln.phy",
                "-t",
                "diffsel_tree.nwk",
                "-ncond",
                str(conds),
                "-x",
                "1 3000",
                "run1"
                ]
        cmd2 = ["diffsel",
                "-d",
                "diffsel_aln.phy",
                "-t",
                "diffsel_tree.nwk",
                "-ncond",
                str(conds),
                "-x",
                "1 3000",
                "run2"
                ]
        cmds["diffsel1"] = cmd1
        cmds["diffsel2"] = cmd2

    if "msd" in args.methods:
        # make files
        try:
            os.mkdir("msd")
        except FileExistsError:
            sys.exit("An msd run already exists. "
                     "Delete to restart.")
        prep_for_msd(curroot, scenarios, "msd/msd_tree.nwk",
                     "msd/msd_states.tsv")
        cmd = ["msd",
               "-o",
               "msd/output.txt",
               "msd/msd_tree.nwk",
               "msd/msd_states.tsv",
               args.AA_aln
               ]
        cmds["msd"] = cmd

    if "PCOC" in args.methods:
        # make files
        try:
            os.mkdir("PCOC")
        except FileExistsError:
            sys.exit("A PCOC run already exists. "
                     "Delete to restart.")
        dockerStr = "docker run -e LOCAL_USER_ID=`id -u $USER` --rm -v \
                     $PWD:$PWD -e CWD=$PWD carinerey/pcoc:v1.0.1"
        cmd = [dockerStr,  # hard-coding a lot here
               "pcoc_det.py",
               "-t",
               args.tree,
               "-aa",
               args.AA_aln,
               "-m",
               scenStr[1],
               "-o",
               "PCOC",
               "--plot",
               "--svg",
               "-CATX_est",
               "60",
               "--gamma",
               "--max_gap_allowed",
               "10",
               "--max_gap_allowed_in_conv_leaves",
               "10",
               "--no_cleanup"
               ]
        # tried to run this with shell=False but too hard to get docker to work
        # userID, _ = subprocess.Popen(["id", "-u", os.environ["USER"]],
        #                              stdout=subprocess.PIPE,
        #                              universal_newlines=True).communicate()
        # print(str(userID))
        # cmd = ["docker",
        #        "run",
        #        "-e",
        #        "LOCAL_USER_ID=" + str(userID),
        #        "--rm",
        #        "-v",
        #        TOP + ":" + TOP,
        #        "-e",
        #        "CWD=" + TOP,
        #        "carinerey/pcoc",
        #        "pcoc_det.py",
        #        "-t",
        #        args.tree,
        #        "-aa",
        #        args.AA_aln,
        #        "-m",
        #        scenStr[1],
        #        "-o",
        #        "PCOC",
        #        "--plot",
        #        "--svg",
        #        "-CATX_est",
        #        "60",
        #        "--gamma",
        #        "--max_gap_allowed",
        #        "10",
        #        "--max_gap_allowed_in_conv_leaves",
        #        "10",
        #        ]
        cmds["PCOC"] = cmd

    if "topo" in args.methods:
        if "PCOC" not in args.methods:
            sys.exit("topo must be run with PCOC!")
        # make files
        try:
            os.mkdir("topo")
        except FileExistsError:
            sys.exit("A topo run already exists. "
                     "Delete to restart.")
        # all we need to do here is prep dir, since we
        # steal topos from PCOC

    if "tdg09" in args.methods:
        # make files
        try:
            os.mkdir("tdg09")
        except FileExistsError:
            sys.exit("A tdg09 run already exists. "
                     "Delete to restart.")
        # optimise with WAG+G for tdg09 site model
        optimise(args.tree, args.AA_aln, "WAG+G")
        root_and_clean("WAG+G_tree.raxml.bestTree", outgroup=og)
        with open("WAG+G_tree.raxml.bestTree.rr.cltr", "r") as inf:
            nwkString = inf.readline().strip()
            tdg09root = tree_reader.read_tree_string(nwkString)
        tdg09root.number_tree()
        label_for_tdg09(tdg09root, aa_seqs, scenarios,
                        args.cond0, args.cond1,
                        "tdg09/tdg09_tree.nwk",
                        "tdg09/tdg09_aln.phy")
        tdg09Str = "java -cp \
                    /home/nat/Applications/tdg09/tdg09-1.1.2/dist/tdg09.jar \
                    tdg09.Analyse"
        # cmd = [tdg09Str,
        #        "-alignment",
        #        "tdg09_aln.phy",
        #        "-groups",
        #        args.cond0,
        #        args.cond1,
        #        "-tree",
        #        "tdg09_tree.nwk"
        #        ]
        cmd = ["java",
               "-cp",
               "/home/nat/Applications/tdg09/tdg09-1.1.2/dist/tdg09.jar",
               "tdg09.Analyse",
               "-alignment",
               "tdg09_aln.phy",
               "-groups",
               args.cond0,
               args.cond1,
               "-tree",
               "tdg09_tree.nwk"
               ]
        cmds["tdg09"] = cmd

    for m in args.methods:
        if m == "diffsel":
            subprocess.Popen(" ".join(cmds["diffsel1"]), shell=True,
                             cwd=TOP + "/" + "diffsel",
                             stdout=subprocess.DEVNULL,
                             stderr=subprocess.DEVNULL
                             )
            # try:
            #     diff1.wait(timeout=5)
            # except subprocess.TimeoutExpired:
            #     print("")
            subprocess.Popen(" ".join(cmds["diffsel2"]), shell=True,
                             cwd=TOP + "/" + "diffsel",
                             stdout=subprocess.DEVNULL,
                             stderr=subprocess.DEVNULL
                             )
            # try:
            #     diff2.wait(timeout=5)
            # except subprocess.TimeoutExpired:
            #     print("")
        if m == "msd":
            subprocess.Popen(cmds["msd"], shell=False,
                             stderr=subprocess.DEVNULL,
                             stdout=subprocess.DEVNULL
                             )
        if m == "PCOC":
            print(" ".join(cmds["PCOC"]))
            pcoc = subprocess.Popen(" ".join(cmds["PCOC"]), shell=True,
                                    stderr=subprocess.DEVNULL,
                                    stdout=subprocess.DEVNULL
                                    )
            try:
                pcoc.wait(5)
            except subprocess.TimeoutExpired:
                print("")
        if m == "tdg09":
            f = open("tdg09/tdg09_out.txt", "w")
            subprocess.Popen(cmds["tdg09"], shell=False,
                             cwd=TOP + "/" + "tdg09",
                             stderr=subprocess.DEVNULL,
                             stdout=f)
        if m == "topo":
            with open(glob.glob("PCOC/RUN_*/Trees/tree.nhx")[0], "r") as treF:
                nwkString = treF.readline().strip()
                non_conv_root = tree_reader.read_tree_string(nwkString)
            for n in non_conv_root.iternodes(order="preorder"):
                n.note = ""
            with open("topo/tree.nwk", "w") as outTre:
                outTre.write(non_conv_root.get_newick_repr() + ";\n")

            with open(
                      glob.glob("PCOC/RUN_*/Trees/tree_conv.nhx")[0], "r"
                      ) as treF:
                nwkString = treF.readline().strip()
                conv_root = tree_reader.read_tree_string(nwkString)
            for n in conv_root.iternodes(order="preorder"):
                n.note = ""
            with open("topo/tree_conv.nwk", "w") as outTre:
                outTre.write(conv_root.get_newick_repr() + ";\n")

            os.chdir("topo")
            m = calc_topological.optimize_model(TOP + "/" + args.AA_aln,
                                                "tree.nwk",
                                                "non_conv",
                                                model=None)
            calc_topological.optimize_model(TOP + "/" + args.AA_aln,
                                            "tree_conv.nwk",
                                            "conv",
                                            model=m)
            df_non_conv = calc_topological.tp_to_dict("non_conv.sitelh")
            df_conv = calc_topological.tp_to_dict("conv.sitelh")
            df_final = pd.merge(df_non_conv[["sites", "logl"]],
                                df_conv[["sites", "logl"]],
                                on="sites",
                                suffixes=['_non_conv', '_conv'])
            df_final["topological"] = list(map(calc_topological.prob_ap,
                                               df_final["logl_conv"],
                                               df_final["logl_non_conv"]))
            df_final = df_final[["sites", "topological"]]
            df_final.to_csv("topo_output.txt", sep=",", index=False)
            os.chdir(TOP)
