#! /usr/bin/python3

import sys
import argparse
from parse_fasta import parse_fasta
from itertools import permutations, chain


def read_subs(path):
    """takes as input a tab-delimited substitution file produced by
    summarise_asr_over_tree.py
    Consists of parent  descendant  subs[comma-separated, format A265V]"""
    subsDict = {}
    with open(path, "r") as inf:
        next(inf)  # skip header
        for line in inf.readlines():
            line = line.strip().split("\t")
            parentDescendant = (line[0], line[1])
            try:
                subsStr = line[2].split(",")
            except IndexError:  # no subs
                continue
            subs = [(x[0], int(x[1:-1]), x[-1]) for x in subsStr]
            subsDict[parentDescendant] = subs

    return subsDict


def get_ref_pos_dict(seqDict):
    allCorrDict = {}
    for k in seqDict.keys():
        corrDict = {}
        aln_pos = 1
        seq_pos = 1
        for i in seqDict[k]:
            if i == "-":
                corrDict[aln_pos] = 0
                aln_pos += 1
            else:
                corrDict[aln_pos] = seq_pos
                aln_pos += 1
                seq_pos += 1
        allCorrDict[k] = corrDict
    return allCorrDict


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("subsfile", help="tsv with subs produced by \
                        summarise_asr_over_tree.py")
    parser.add_argument("branches", help="space-separated parent-daughter \
                        comparisons in format parent,daughter",
                        type=str, nargs="+")
    parser.add_argument("-g", "--groups", help="space-separated groups of \
                        branches that represent particular origins, e.g. \
                        N7,N8,N9 N171,N224,N232", type=str, nargs="+")
    parser.add_argument("-a", "--atleast", help="sub must occur in at least \
                        n branches (default 2)", type=int, default=2)
    parser.add_argument("-i", "--isin", help="for overlap mode: list \
                        substitutions in branch parent,child. Can be \
                        space-separated list", type=str, nargs="+")
    parser.add_argument("-n", "--notin", help="for overlap mode: list \
                        substitutions in --in, but not in branch \
                        parent,child", type=str, nargs="+")
    parser.add_argument("-t", "--type", help="classify substitutions as \
                        convergent or coincident", action="store_true")
    parser.add_argument("-rd", "--ref_desc", help="additionally print \
                        substitutions positions relative to the descendant \
                        sequence in each comparison, instead of alignment \
                        position", action="store_true")
    parser.add_argument("-rp", "--ref_par", help="additionally print \
                        substitutions positions relative to the parent \
                        sequence in each comparison, instead of alignment \
                        position", action="store_true")
    parser.add_argument("-aln", "--alignment", help="if reference: \
                        FASTA-formatted alignment to extract position \
                        correspondences", type=str)
    parser.add_argument("-nc", "--nonconv", help="also write \
                        non-site-overlapping changes on each branch",
                        action="store_true")
    # parser.add_argument("-d", "--degree", help="number of origins required \
    #                     to call a convergent substitution (minimum 2, \
    #                     default)", type=int, default=2)
    args = parser.parse_args()

    branches = [(x.split(",")[0], x.split(",")[1]) for x in args.branches]
    if len(branches) == 1:
        args.isin = args.branches

    subs = read_subs(args.subsfile)

    if args.ref_desc or args.ref_par:
        if not args.alignment:
            sys.stderr.write("must supply alignment for reference mode\n")
            sys.exit()
        seqs = dict([x for x in parse_fasta(args.alignment)])
        corres = get_ref_pos_dict(seqs)

    if args.isin and not args.notin:
        inNodes = [(x.split(",")[0], x.split(",")[1]) for x in args.isin]
        allSubsPos = []
        for b in inNodes:
            allSubsPos.append(set([x[1] for x in subs[b]]))
        allSubsPosUnion = sorted(list(set.union(*allSubsPos)))
        nodes = inNodes
    elif args.isin and args.notin:
        inNodes = [(x.split(",")[0], x.split(",")[1]) for x in args.isin]
        notInNodes = [(x.split(",")[0], x.split(",")[1]) for x in args.notin]
        inSubsPos = []
        for b in inNodes:
            inSubsPos += [x[1] for x in subs[b]]
        notInSubsPos = []
        for b in notInNodes:
            notInSubsPos += [x[1] for x in subs[b]]
        allSubsPosUnion = sorted(list(set(inSubsPos) - set(notInSubsPos)))
    else:
        branchCombs = permutations(branches, args.atleast)

        allSubsPos = []
        for b in branchCombs:
            if args.groups is not None:
                contained = False
                for g in args.groups:
                    g = [n for n in g.split(",")]
                    if set(chain.from_iterable(b)).issubset(set(g)):
                        sys.stderr.write("skipping %s contained in %s\n"
                                         % (b, g))
                        contained = True
                if contained:
                    continue
            subsPos = []
            for n in b:
                subsPos.append(set([x[1] for x in subs[n]]))
            subsInAll = set.intersection(*subsPos)
            allSubsPos.append(subsInAll)

        allSubsPosUnion = sorted(list(set.union(*allSubsPos)))

    if args.ref_desc:
        print("\t".join(["pos_aln",
                         "pos_desc",
                         "par",
                         "desc",
                         "sub",
                         "type"]))
    if args.ref_par:
        print("\t".join(["pos_aln",
                         "pos_par",
                         "par",
                         "desc",
                         "sub",
                         "type"]))
    for pos in allSubsPosUnion:
        if args.type:
            endStates = []
            for b in branches:
                endStates += [s[2] for s in subs[b] if pos == s[1]]
                if len(endStates) == 0:
                    continue
            for b in branches:
                subsInPos = [s[0] + str(s[1]) + s[2] for s in subs[b]
                             if pos == s[1]]
                if len(subsInPos) == 0:
                    continue
                # if args.groups is not None:
                #     uniq = len(branches)
                #     for g in args.groups:
                #         g = [n for n in g.split(",")]
                #         if set(chain.from_iterable(b)).issubset(set(g)):
                #             uniq -= 1
                #     if uniq < args.atleast:
                #         print("Uniq less than at least")
                #         continue
                endState = [s[2] for s in subs[b] if pos == s[1]][0]
                if endStates.count(endState) >= args.atleast:
                    """make changes here in future to add requirements of 2-,
                    3-way overlap etc."""
                    if args.ref_desc:
                        print("\t".join([str(pos),
                                         str(corres[b[1]][pos]),
                                         b[0],
                                         b[1],
                                         ",".join(subsInPos), "CONV"]))
                    elif args.ref_par:
                        print("\t".join([str(pos),
                                         str(corres[b[0]][pos]),
                                         b[0],
                                         b[1],
                                         ",".join(subsInPos), "CONV"]))
                    else:
                        print("\t".join([str(pos), b[0], b[1],
                                        ",".join(subsInPos), "CONV"]))
                else:
                    if args.ref_desc:
                        print("\t".join([str(pos),
                                         str(corres[b[1]][pos]),
                                         b[0],
                                         b[1],
                                         ",".join(subsInPos), "COIN"]))
                    elif args.ref_par:
                        print("\t".join([str(pos),
                                         str(corres[b[0]][pos]),
                                         b[0],
                                         b[1],
                                         ",".join(subsInPos), "COIN"]))
                    else:
                        print("\t".join([str(pos), b[0], b[1],
                                        ",".join(subsInPos), "COIN"]))
        else:
            for b in branches:
                subsInPos = [s[0] + str(s[1]) + s[2] for s in subs[b]
                             if pos == s[1]]
                if len(subsInPos) == 0:
                    continue
                else:
                    if args.ref_desc:
                        print("\t".join([str(pos),
                                         str(corres[b[1]][pos]),
                                         b[0],
                                         b[1],
                                         ",".join(subsInPos)]))
                    elif args.ref_par:
                        print("\t".join([str(pos),
                                         str(corres[n[0]][pos]),
                                         b[0],
                                         b[1],
                                         ",".join(subsInPos)]))
                    else:
                        print("\t".join([str(pos), b[0], b[1],
                                        ",".join(subsInPos)]))
    if args.nonconv:
        for b in branches:
            for s in subs[b]:
                pos = s[1]
                if pos not in allSubsPosUnion:
                    if args.ref_desc:
                        print("\t".join([str(pos),
                                         str(corres[b[1]][pos]),
                                         b[0],
                                         b[1],
                                         f"{s[0]}{s[1]}{s[2]}", "NONC"]))
                    elif args.ref_par:
                        print("\t".join([str(pos),
                                         str(corres[b[0]][pos]),
                                         b[0],
                                         b[1],
                                         f"{s[0]}{s[1]}{s[2]}", "NONC"]))
                    else:
                        print("\t".join([str(pos), b[0], b[1],
                                        f"{s[0]}{s[1]}{s[2]}", "NONC"]))
    # subsPos = []
    # for n in nodes:
    #     subsPos.append(set([x[1] for x in subs[n]]))

    # subsInAll = set.intersection(*subsPos)

    # print(subsInAll)
