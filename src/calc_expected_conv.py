#!/usr/bin/env python3


import sys
import argparse
import tree_reader
import numpy as np
from treenode import Node
from collections import Counter
from itertools import combinations
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
    for k, v in seq_dict.items():
        for i in v:
            try:
                col_dict[pos][k] = i
            except KeyError:
                col_dict[pos] = {}
                col_dict[pos][k] = i
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
        tot += sum([x[1] for x in counts.items() if x[0] in AA])
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
    for k, v in col_dict.items():
        with open(f"tmp_{k}.aln", "w") as sitef:
            for header, char in v.items():
                sitef.write(f">{header}\n")
                sitef.write(f"{char}\n")
        freqs = get_optimised_freqs(tree, f"tmp_{k}.aln", "JTT+FO")
        out[k] = freqs
    return out


def parse_site_frequencies_file(inf: str) -> dict:
    


def parse_paml_rates(path: str) -> dict:
    """Takes a paml 'rates' file containing estimated posterior mean rate per
    site and returns a dictionary of {site: rate}"""
    out = {}
    with open(path, "r") as inf:
        going = False
        reading = False
        for line in inf:
            if line.strip().startswith("Site"):
                going = True
                continue
            elif line == "\n" and going and not reading:
                reading = True
                continue
            elif line == "\n" and going and reading:
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
    idx = AAIDX[aa]
    array = np.zeros(20)
    array[idx] = 1.
    return array


def get_path_length_mrca(tree: Node, node1: str, node2: str):
    """returns the sum of the branch lengths between two labelled nodes and
    their mrca"""
    mrca = None
    lengths = []
    for n in tree.iternodes():
        if n.label == node1:
            going = True
            brlen = 0.0
            while going:
                p = n
                brlen += n.length
                if node2 in [i.label for i in p.iternodes()]:
                    going = False
                    mrca = p
                else:
                    n = p.parent
    lengths.append(brlen)
    for n in tree.iternodes():
        if n.label == node2:
            going = True
            brlen = 0.0
            while going:
                p = n
                brlen += n.length
                if node1 in [i.label for i in p.iternodes()]:
                    going = False
                else:
                    n = p.parent
    lengths.append(brlen)
    return lengths, mrca


def get_path_length_root(tree: Node, node1: str):
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

    ancs = dict([x for x in parse_fasta(args.ancestors)])
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
                                    f"extant sequences need to be included")
                    sys.exit()
        modJTT.calc_empirical_freqs(extants)
        modJTT.set_empirical_freqs()
    else:
        modJTT.set_model_freqs()
    modJTT.scale_rate_matrix()

    if args.site_frequencies:
        if args.site_frequencies_file:

        extants = {}
        for n in curroot.iternodes():
            if n.istip:
                try:
                    extants[n.label] = ancs[n.label]
                except KeyError:
                    sys.stderr.write(f"{n.label} not in provided FASTA - "
                                    f"extant sequences need to be included")
                    sys.exit()
        extant_cols = get_columns(extants)
        site_freqs = get_opt_site_specific_frequencies(extant_cols,
                                                       args.tree)
        with open("site_frequencies.tsv", "w") as siteff:
            for k, v in site_freqs.items():
                siteff.write(f"{k}\t{','.join([str(f) for f in v])}\n")

    branches = [(x.split(",")[0], x.split(",")[1]) for x in args.branches]
    if len(branches) == 1:
        print("More than one branch required")
        sys.exit()

    if len(branches) == 2:
        branch_combs = [tuple(branches)]
    else: 
        branch_combs = [x for x in combinations(branches, 2)]

    good_combs = []
    for b in branch_combs:
        remove = False
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
        conv_prob_sum = 0.0
        if b[0][0] == b[1][0]:
            # print("this happens")
            continue  # place holder for the special case
        elif node_dict[b[0][1]].istip or node_dict[b[1][1]].istip:
            # in this case, the conditional probabilities that node 2 or 4 
            # are in any state will be 0 unless the tip is observed to be in 
            # that state
            continue
        else:
            lengths, mrca = get_path_length_mrca(curroot, b[0][0], b[1][0])
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
                            continue
                    if not args.rates:  # do these once and save time
                        p0 = modJTT.get_P(lengths[0])
                        p1 = modJTT.get_P(lengths[1])
                        p01 = modJTT.get_P(node_dict[b[0][1]].length)
                        p11 = modJTT.get_P(node_dict[b[1][1]].length)
                # dict for probs
                sites += 1
                node_probs = {}
                # probs for parents
                map_array = map_state_array(map_state)
                if args.rates:
                    node_probs[b[0][0]] = np.matmul(map_array,
                                                    modJTT.get_P(lengths[0],
                                                                 rates_dict[k])
                                                   )  #1x20
                    node_probs[b[1][0]] = np.matmul(map_array,
                                                    modJTT.get_P(lengths[1],
                                                                 rates_dict[k])
                                                   )  #1x20
                else:
                    node_probs[b[0][0]] = np.matmul(map_array, p0)
                    node_probs[b[1][0]] = np.matmul(map_array, p1)
                # conditional probs for children
                conditionals_l = np.zeros((20,20))
                conditionals_r = np.zeros((20,20))
                for i in range(20):
                    map_array = np.zeros(20)
                    map_array[i] = 1.
                    if args.rates:
                        prob_l = np.matmul(map_array,
                                           modJTT.get_P(node_dict[b[0]
                                                                  [1]].length,
                                                        rates_dict[k]))  # 1x20
                        conditionals_l[i] = prob_l
                        prob_r = np.matmul(map_array,
                                           modJTT.get_P(node_dict[b[1]
                                                                  [1]].length,
                                                        rates_dict[k]))  # 1x20
                        conditionals_r[i] = prob_r
                    else:
                        prob_l = np.matmul(map_array, p01)  # 1 x 20
                        conditionals_l[i] = prob_l
                        prob_r = np.matmul(map_array, p11)  # 1 x 20
                        conditionals_r[i] = prob_r
                        # conditionals now holds the probability of desc being 
                        # in state j (col) given starting in state i (row)
                for i in range(20):
                    for j in range(20):
                        for k in range(20):
                            for l in range(20):
                                if j == l and i != j and k != l:
                                    # no restriction of i and k - sum over
                                    # parallel and convergent
                                    # print(f"{b[0][0]} is in {AA[i]}\n"
                                    #       f"descendent {b[0][1]} probability"
                                    #       f"of {AA[j]} given {b[0][0]} is in"
                                    #       f" {AA[i]}\n"
                                    #       f"{b[1][0]} is in {AA[k]}\n"
                                    #       f"descendent {b[1][1]} probability"
                                    #       f"of {AA[l]} given {b[1][0]} is in" 
                                    #       f" {AA[k]}\n")
                                    conv_prob_sum += (node_probs[b[0][0]][i] * 
                                                      conditionals_l[i, j] * 
                                                      node_probs[b[1][0]][k] *
                                                      conditionals_r[k, l])
                
        print(b)
        print(conv_prob_sum)
        print(sites)


                
                
                

