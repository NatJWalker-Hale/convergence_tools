#!/usr/bin/env python3


import sys
import argparse
from collections import Counter
from itertools import combinations
import numpy as np
import tree_reader
from treenode import Node
from parse_fasta import parse_fasta
from discrete_models import Discrete_model
from optimise_raxmlng import get_optimised_freqs


AAIDX = {"A": 0, "R": 1, "N": 2, "D": 3, "C": 4, "Q": 5, "E": 6, "G": 7,
         "H": 8, "I": 9, "L": 10, "K": 11, "M": 12, "F": 13, "P": 14, "S": 15,
         "T": 16, "W": 17, "Y": 18, "V": 19}
AA = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F',
      'P', 'S', 'T', 'W', 'Y', 'V']


def get_columns(seq_dict: dict) -> dict:
    """takes a dictionary of an alignment (key: name, value: sequence),
    and returns a dictionary of columns 
    (key: position, value: dict{name: state})"""
    col_dict = {}
    pos = 0
    for header, seq in seq_dict.items():
        for char in seq:
            try:
                col_dict[pos][header] = char
            except KeyError:
                col_dict[pos] = {}
                col_dict[pos][header] = char
            pos += 1
        pos = 0
    return col_dict


def get_site_specific_frequencies(col_dict: dict) -> dict:
    """takes a dictionary from get_columns and returns a dictionary of 
    site-specific frequencies"""
    freq_dict = {}  # key is pos, value is np array of freqs
    for pos, column in col_dict.items():
        freqs = []
        tot = 0
        counts = Counter([c for c in column.values()])
        tot += sum(x[1] for x in counts.items() if x[0] in AA)
        for state in AA:
            freqs.append(counts[state] / tot)
        freqs = np.array(freqs)
        freq_dict[pos] = freqs
    return freq_dict


def get_opt_site_specific_frequencies(col_dict: dict, tree: str) -> dict:
    """takes a dictionary from get_columns, writes individual column FASTAs,
    optimises ML frequencies with raxml-ng and returns the frequencies as a 
    dictionary"""
    out = {}
    for pos, col in col_dict.items():
        with open(f"tmp_{pos}.aln", "w", encoding='utf-8') as sitef:
            for header, char in col.items():
                sitef.write(f">{header}\n")
                sitef.write(f"{char}\n")
        freqs = get_optimised_freqs(tree, f"tmp_{pos}.aln", "JTT+FO")
        out[pos] = freqs
    return out


def parse_site_frequencies_file(inf: str) -> dict:
    """parse a TSV of pos\tfreqs, where freqs is comma-separated"""
    freqs = {}
    with open(inf, "r", encoding="utf-8") as sff:
        for line in sff:
            line = line.strip().split("\t")
            freqs[int(line[0])] = np.array([float(x) for
                                            x in line[1].split(",")])
    return freqs


def parse_paml_rates(path: str) -> dict:
    """Takes a paml 'rates' file containing estimated posterior mean rate per
    site and returns a dictionary of {site: rate}"""
    out = {}
    with open(path, "r", encoding='utf-8') as inf:
        going = False
        reading = False
        for line in inf:
            if line.strip().startswith("Site"):
                going = True
                continue
            if line == "\n" and going and not reading:
                reading = True
                continue
            if line == "\n" and going and reading:
                break

            if going and reading:
                line = line.strip().split()
                out[int(line[0])] = float(line[-2])
    return out


def get_mrca(tree: Node, node1: str, node2: str) -> Node:
    """returns the node that is the mrca for two labelled nodes"""
    for n in tree.iternodes():
        if n.label == node1:
            going = True
            while going:
                p = n
                if node2 in [i.label for i in p.iternodes()]:
                    going = False
                else:
                    n = p.parent
    return p


def map_state_array(aa: str) -> np.array:
    """returns a dim20 array of zeros with idx aa = 1"""
    idx = AAIDX[aa]
    array = np.zeros(20)
    array[idx] = 1.
    return array


def get_path_length_mrca(tree: Node, node1: str, node2: str) -> tuple:
    """returns the sum of the branch lengths between two labelled nodes and
    their mrca"""
    mrca = None
    lengths = []
    for n in tree.iternodes():
        if n.label == node1:
            going = True
            brlen = 0.
            while going:
                if node2 in [i.label for i in n.iternodes()]:  # node IS mrca
                    going = False
                    mrca = n
                else:
                    brlen += n.length
                    n = n.parent
    lengths.append(brlen)
    for n in tree.iternodes():
        if n.label == node2:
            going = True
            brlen = 0.
            while going:
                if node1 in [i.label for i in n.iternodes()]:  # node IS mrca
                    going = False
                    mrca = n
                else:
                    brlen += n.length
                    n = n.parent
    lengths.append(brlen)
    return lengths, mrca


def get_path_length_root(tree: Node, node1: str):
    """returns the sum of branch lengths between the node and the root"""
    for n in tree.iternodes():
        if n.label == node1:
            going = True
            brlen = 0.
            while going:
                p = n
                brlen += n.length
                if p.parent is None:  # root
                    going = False
                else:
                    n = p.parent
    return brlen


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("tree", help="Tree file containing newick-formatted \
                        tree to run analysis on")
    parser.add_argument("ancestors", help="FASTA-formatted alignment of MAP \
                        sequences for nodes, with names matching labels \
                        produced by --label")
    parser.add_argument("branches", help="space-separated parent-daughter \
                        comparisons in format parent,daughter (at least two)",
                        type=str, nargs="+")
    group = parser.add_mutually_exclusive_group(required = False)
    group.add_argument("-f", "--empirical_frequencies", help="Use \
                        frequencies calculated from extant sequences. \
                        'ancestors' should contain sequences for both \
                        nodes and tips", action="store_true")
    group.add_argument("-sf", "--site_frequencies", help="Use per-site \
                        frequencies calculated from extant sequences. \
                        'ancestors' should contain sequences for both \
                        nodes and tips", action="store_true")
    parser.add_argument("-sff", "--site_frequencies_file", help="if using \
                        site frequencies, read from file. Format: tab \
                        separated, first column position (0-indexed) \
                        second column comma separated frequencies (20)")
    parser.add_argument("-r", "--rates", help="File containing estimated \
                        posterior mean rate for each site from PAML")
    parser.add_argument("-l", "--label", help="Print labelled tree and \
                        quit", action="store_true")
    args = parser.parse_args()
    if args.site_frequencies_file:
        args.site_frequencies = True

    curroot = next(tree_reader.read_tree_file_iter(args.tree))

    curroot.number_nodes()
    # make dict of tree nodes indexed by label
    node_dict = {}
    for n in curroot.iternodes():
        node_dict[n.label] = n

    if args.label:
        print(f"{curroot.get_newick_repr(False, False)};")
        sys.exit()

    if args.rates:
        rates_dict = parse_paml_rates(args.rates)

    ancs = dict(parse_fasta(args.ancestors))
    anc_cols = get_columns(ancs)

    # set up model
    modJTT = Discrete_model()
    modJTT.set_rate_JTT()
    if args.empirical_frequencies:
        extants = {}
        for n in curroot.iternodes():
            if n.istip:
                try:
                    extants[n.label] = ancs[n.label]
                except KeyError:
                    sys.stderr.write(f"{n.label} not in provided FASTA - "
                                     f"extant sequences need to be included\n")
                    sys.exit()
        modJTT.calc_empirical_freqs(extants)
        modJTT.set_empirical_freqs()
    else:
        modJTT.set_model_freqs()
    modJTT.scale_rate_matrix()

    if args.site_frequencies:
        if args.site_frequencies_file:
            site_freqs = parse_site_frequencies_file(args.site_frequencies_file)
        else:
            extants = {}
            for n in curroot.iternodes():
                if n.istip:
                    try:
                        extants[n.label] = ancs[n.label]
                    except KeyError:
                        sys.stderr.write(f"{n.label} not in provided FASTA - "
                                         f"extant sequences need to be "
                                         f"included\n")
                        sys.exit()
            extant_cols = get_columns(extants)
            site_freqs = get_opt_site_specific_frequencies(extant_cols,
                                                        args.tree)
            with open("site_frequencies.tsv", "w", encoding="utf-8") as siteff:
                for k, v in site_freqs.items():
                    siteff.write(f"{k}\t{','.join([str(f) for f in v])}\n")

    branches = [(x.split(",")[0], x.split(",")[1]) for x in args.branches]
    if len(branches) == 1:
        print("More than one branch required")
        sys.exit()

    if len(branches) == 2:
        branch_combs = [tuple(branches)]
    else:
        branch_combs = list(combinations(branches, 2))

    good_combs = []
    for b in branch_combs:
        remove = False
        if b[0][0] == b[1][0]:  # sisters
            sys.stderr.write(f"skipping combination {b} as {b[0]} is "
                             f"sister to {b[1]}\n")
            remove = True
        # start arbitrarily with left parent
        n = node_dict[b[0][0]]
        going = True
        while going:
            p = n.parent
            if p.parent is None:  # reach root without encountering
                going = False
            if p.label == b[1][1]:  # right descendant is in parents
                sys.stderr.write(f"skipping combination {b} as {b[1]} is "
                                 f"direct ancestor of {b[0]}\n")
                remove = True
            else:
                n = p

        n = node_dict[b[0][1]]  # left descendant
        for c in n.iternodes():
            if c.label == b[1][0]:  # right parent is in descendants
                sys.stderr.write(f"skipping combination {b} as {b[1]} is "
                                 f"direct descendant of {b[0]}\n")
                remove = True

        if not remove:
            good_combs.append(b)


    for b in good_combs:
        conv_prob_sum = 0.
        div_prob_sum = 0.
        lengths, mrca = get_path_length_mrca(curroot, b[0][0], b[1][0])
        if mrca.label in (b[0][0], b[1][0]):
            # if either of the nodes is the MRCA, we calculate probs for
            # this node by conditioning on the state one node back, to
            # avoid fixed states
            lengths = [x + mrca.length for x in lengths]
            mrca = mrca.parent
        # get lengths between node and MRCA
        sites = 0
        if not args.rates:  # do these once and save time
            p0 = modJTT.get_P(lengths[0])
            p1 = modJTT.get_P(lengths[1])
            p01 = modJTT.get_P(node_dict[b[0][1]].length)
            p11 = modJTT.get_P(node_dict[b[1][1]].length)
        for k, v in anc_cols.items():   # going by site
            map_state = v[mrca.label]  # get MAP state of MRCA
            if map_state == "-":
                continue  # skip cases where mrca MAP is gap
            nodes = [i for j in b for i in j]
            states = [v[n] for n in nodes]
            if "-" in states:
                continue  # skip cases where any MAPs are inferred as gaps
            if args.site_frequencies:
                # set up model for site specific freqs
                modJTT = Discrete_model()
                modJTT.set_rate_JTT()
                modJTT.set_frequencies(site_freqs[k])
                with np.errstate(divide="raise"):
                    # account for cases where 0 freqs cause problems
                    try:
                        modJTT.scale_rate_matrix()
                    except FloatingPointError:
                        sys.stderr.write(f"skipping site {k} where "
                                            f"sparse frequencies caused a "
                                            f"problem\n")
                        continue
            if args.rates:
                p0 = modJTT.get_P(lengths[0], rates_dict[k])
                p1 = modJTT.get_P(lengths[1], rates_dict[k])
                p01 = modJTT.get_P(node_dict[b[0][1]].length,
                                    rates_dict[k])
                p11 = modJTT.get_P(node_dict[b[1][1]].length,
                                    rates_dict[k])
            else:
                p0 = modJTT.get_P(lengths[0])
                p1 = modJTT.get_P(lengths[1])
                p01 = modJTT.get_P(node_dict[b[0][1]].length)
                p11 = modJTT.get_P(node_dict[b[1][1]].length)
            # dict for probs
            sites += 1
            node_probs = {}
            # probs for parents
            map_array = map_state_array(map_state)
            node_probs[b[0][0]] = np.matmul(map_array, p0)  #1x20
            node_probs[b[1][0]] = np.matmul(map_array, p1)  #1x20
            # conditional probs for children
            conditionals_l = np.zeros((20,20))
            conditionals_r = np.zeros((20,20))
            for i in range(20):
                map_array = np.zeros(20)
                map_array[i] = 1.
                prob_l = np.matmul(map_array, p01)  # 1x20
                # if node_dict[b[0][1]].istip:
                #     prob_l = (prob_l *
                #                 map_state_array(v[node_dict[b[0][1]].label]))
                conditionals_l[i] = prob_l
                prob_r = np.matmul(map_array, p11)  # 1x20
                # if node_dict[b[1][1]].istip:
                #     prob_r = (prob_r *
                #                 map_state_array(v[node_dict[b[1][1]].label]))
                conditionals_r[i] = prob_r
                # conditionals now holds the probability of desc being
                # in state j (col) given starting in state i (row)
                joint_probs_l = np.matmul(np.diag(node_probs[b[0][0]]),
                                            conditionals_l)
                joint_probs_r = np.matmul(np.diag(node_probs[b[1][0]]),
                                            conditionals_r)
            # if node_dict[b[1][1]].istip:
            #     print(conditionals_r)
            #     print(joint_probs_r)
            for i in range(20):  # parent 1
                for j in range(20):  # child 1
                    if i != j:  # only bother enumerating if these are
                                # different
                        for k in range(20):  # parent 2
                            for l in range(20):  # child 2
                                if l != k and j == l:
                                # print(f"Branch 1: {AA[i]} -> {AA[j]}")
                                # print(f"Branch 2: {AA[k]} -> {AA[j]}")
                                    # conv_prob_sum += (node_probs[b[0][0]][i] *
                                    #                   conditionals_l[i, j] *
                                    #                   node_probs[b[1][0]][k] *
                                    #                   conditionals_r[k, l])
                                    conv_prob_sum += (joint_probs_l[i, j] *
                                                        joint_probs_r[k, l])
                                elif l not in (k, j):
                                    # div_prob_sum += (node_probs[b[0][0]][i] *
                                    #                  conditionals_l[i, j] *
                                    #                  node_probs[b[1][0]][k] *
                                    #                  conditionals_r[k, l])
                                    div_prob_sum += (joint_probs_l[i, j] *
                                                     joint_probs_r[k, l])
            # for i in range(20):
            #     for j in range(20):
            #         for k in range(20):
            #             for l in range(20):
            #                 if j == l and i != j and k != l:
            #                     # no restriction of i and k - sum over
            #                     # parallel and convergent
            #                     # print(f"{b[0][0]} is in {AA[i]}\n"
            #                     #       f"descendent {b[0][1]} probability"
            #                     #       f"of {AA[j]} given {b[0][0]} is in"
            #                     #       f" {AA[i]}\n"
            #                     #       f"{b[1][0]} is in {AA[k]}\n"
            #                     #       f"descendent {b[1][1]} probability"
            #                     #       f"of {AA[l]} given {b[1][0]} is in"
            #                     #       f" {AA[k]}\n")
            #                     conv_prob_sum += (node_probs[b[0][0]][i] *
            #                                       conditionals_l[i, j] *
            #                                       node_probs[b[1][0]][k] *
            #                                       conditionals_r[k, l])
            #                 elif j != l and i != j and k != l:
            #                     div_prob_sum += (node_probs[b[0][0]][i] *
            #                                       conditionals_l[i, j] *
            #                                       node_probs[b[1][0]][k] *
            #                                       conditionals_r[k, l])

        print(b)
        print(f"{conv_prob_sum}\t{div_prob_sum}")
        print(sites)
