#!/usr/bin/env python3


"""given an input file of site-specific frequencies, use seq-gen to repeatedly
simulate single sites and concatenate"""


import sys
import argparse
import subprocess
import numpy as np
from treenode import Node
from parse_fasta import parse_fasta, parse_fasta_str
from calc_expected_conv import 


def parse_site_frequencies_file(inf: str) -> dict:
    """parse a TSV of pos\tfreqs, where freqs is comma-separated"""
    freqs = {}
    with open(inf, "r", encoding="utf-8") as sff:
        for line in sff:
            line = line.strip().split("\t")
            freqs[int(line[0])] = np.array([float(x) for
                                            x in line[1].split(",")])
    return freqs


def sim_sites(treef: str, freqs: np.array, model:str="JTT", use_anc:int=None,
             rate:float=1.0, sites:int=1) -> dict:
    """
    Simulates a single site with seq-gen given a treefile, which may
    optionally include sequences to serve as ancestral state, in which case
    treef should be a file with the following formatting:

    2 4
    taxon1 ATGC
    taxon2 ATGC
    1
    (taxon1:0.1,taxon2:0.1);

    rate information (e.g. from PAML posterior mean rates) can be incorporated
    using a branch length scalar
    """
    models = ["JTT", "WAG", "PAM", "BLOSUM", "MTREV", "CPREV45",
              "MTART", "LG", "GENERAL"]
    if model not in models:
        sys.stderr.write(f"{model} not in {models}\n")
        sys.exit()
    cmd = [
            "seq-gen",
            "-l",
            f"{sites}",
            "-m"
            f"{model}",
            "-f",
            f"{' '.join([str(x) for x in freqs])}"
            "-of",
            "-wa"
        ]
    if use_anc is not None:
        cmd += ["-k", f"{use_anc}"]
    if rate != 1.0:
        cmd += ["-s", f"{rate}"]
    cmd += [treef]
    print(subprocess.list2cmdline(cmd))
    res = subprocess.run(cmd, shell=False, capture_output=True, text=True,
                         check=False)
    out = dict(parse_fasta_str(res.stdout))
    return out


def concat_seqs(seq_dict1: dict, seq_dict2: dict) -> dict:
    """concatenate contents of two dictionaries with matching keys into a
    third"""
    out = {}
    for header, seq in seq_dict1.items():
        try:
            out[header] = seq + seq_dict2[header]
        except KeyError:
            sys.stderr.write(f"{header} not in {seq_dict2}\n")
            sys.exit()
    return out


def write_seq_gen_input(tree: Node, seq_dict:dict=None):
    """write input tree (and optionally, sequence) file for seq-gen"""
    with open("seq_gen_input.txt", "w", encoding="utf-8") as outf:
        if seq_dict is not None:
            outf.write(f"{len(seq_dict)} {len(seq_dict.values()[0])}\n")
            for header, seq in seq.items():
                outf.write(f"{header} {seq}\n")
            outf.write("1\n")
        outf.write(f"{tree.get_newick_repr(showbl=True, shownum=False)};")


def sim_alignment_gene_freqs(tree: Node, seq_dict:dict,
                             freqs=np.array, use_anc:int=None) -> dict:
    