#!/usr/bin/env python3


"""
writes inferred convergent, divergent and other substitutions from ASR per site
"""


import sys
import argparse
from itertools import combinations
from calc_expected_conv import get_good_branch_combs
from count_conv_subs_from_asr import add_subs
import newick as nwk
import sequence as sq


def read_subs(path) -> dict[tuple: list[tuple[str, int, str]]]:
    """
    takes as input a tab-delimited substitution file produced by summarise_asr_over_tree.py
    Consists of parent  descendant  subs[comma-separated, format A265V]
    """
    subs_dict = {}
    with open(path, "r", encoding="utf-8") as inf:
        next(inf)  # skip header
        for line in inf.readlines():
            line = line.strip().split("\t")
            par_desc = (line[0], line[1])
            try:
                subs_str = line[2].split(",")
            except IndexError:  # no subs
                continue
            subsl = [(x[0], int(x[1:-1]), x[-1]) for x in subs_str]
            subs_dict[par_desc] = subsl

    return subs_dict


def get_ref_pos_dict(seq_dict: dict[str: str]) -> dict[str: dict[int: int]]:
    """
    returns a dictionary of {aln_pos: seq_pos} for all sequences in seq_dict
    """
    all_corr_dict = {}
    for header in seq_dict:
        corr_dict = {}
        aln_pos = 1
        seq_pos = 1
        for i in seq_dict[header]:
            if i == "-":
                corr_dict[aln_pos] = 0
                aln_pos += 1
            else:
                corr_dict[aln_pos] = seq_pos
                aln_pos += 1
                seq_pos += 1
        all_corr_dict[header] = corr_dict
    return all_corr_dict


# def get_pos_overlap(branch_comb: tuple[tuple[str, str]],
#                     subs_dict: dict[tuple: list[tuple[str, int, str]]]) -> set:
#     print(branch_comb)
#     nodes = [branch for branch in branch_comb]
#     print(nodes)



if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("alignment", help="FASTA-formatted alignment including ancestral seqs")
    parser.add_argument("tree", help="tree with node labels matching headers in alignment")
    parser.add_argument("branches", help="space-separated parent-daughter comparisons in format \
                        parent,daughter (can be multibranch lineages)", type=str, nargs="+")
    parser.add_argument("-a", "--atleast", help="sub must occur in at least n branches \
                        (default 2)", type=int, default=2)
    parser.add_argument("-i", "--isin", help="for overlap mode: list substitutions in branch \
                        parent,child. Can be space-separated list", type=str, nargs="+")
    parser.add_argument("-n", "--notin", help="for overlap mode: list substitutions in --in, but \
                        not in branch parent,child", type=str, nargs="+")
    parser.add_argument("-rd", "--ref_desc", help="additionally print substitutions positions \
                        relative to the descendant sequence in each comparison",
                        action="store_true")
    parser.add_argument("-rp", "--ref_par", help="additionally print substitutions positions \
                        relative to the parent sequence in each comparison", action="store_true")
    # parser.add_argument("-d", "--degree", help="number of origins required \
    #                     to call a convergent substitution (minimum 2, \
    #                     default)", type=int, default=2)
    args = parser.parse_args()

    alignment = dict(sq.parse_fasta(args.alignment))

    curroot = nwk.parse_from_file(args.tree)

    # make dict of tree nodes indexed by label
    node_dict = {}
    for n in curroot.iternodes():
        node_dict[n.label] = n

    branches = [(x.split(",")[0], x.split(",")[1]) for x in args.branches]
    print(branches)
    if len(branches) == 1:
        args.isin = args.branches

    subs = {br: [] for br in branches}

    add_subs(subs, alignment)

    corres = get_ref_pos_dict(alignment)

    if args.isin and not args.notin:
        inNodes = [(x.split(",")[0], x.split(",")[1]) for x in args.isin]
        allSubsPos = []
        for b in inNodes:
            allSubsPos.append({x[1] for x in subs[b]})
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
        combs = combinations(branches, 2)
        good_combs = get_good_branch_combs(combs, curroot)

        allSubsPos = []
        for br_comb in good_combs:
            subsPos = []
            for n in br_comb:
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
                # make changes here in future to add requirements of 2-,
                # 3-way overlap etc.
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
