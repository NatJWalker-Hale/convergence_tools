#!/usr/bin/env python3


"""given an input file of site-specific frequencies, use seq-gen to repeatedly
simulate single sites and concatenate"""


import sys
import argparse
import subprocess
from collections import Counter
import numpy as np
from phylo import Node
from parse_phylip import parse_phylip_str
from calc_expected_conv import get_columns, get_site_specific_frequencies


def check_aligned(seq_dict: dict) -> bool:
    """
    checks if sequences in a sequence dictionary are the same length
    """
    if len(set(len(s) for s in seq_dict.values())) > 1:
        return False
    return True


def calc_gene_frequencies(seq_dict: dict) -> np.array:
    """
    calculates the frequencies of amino acid character states in sequence dictionary, ignoring gaps
    """
    aas = "ARNDCQEGHILKMFPSTWYV"
    counts = Counter(char for seq in seq_dict.values() for char in seq if
                     char in set(aas))
    tot = sum(counts.values())
    frequencies = np.fromiter((counts[aa] for aa in aas), dtype=float) / tot
    return frequencies


def parse_site_frequencies_file(inf: str) -> dict:
    """
    parse a TSV of pos\tfreqs, where freqs is comma-separated
    """
    frequencies = {}
    with open(inf, "r", encoding="utf-8") as sff:
        for line in sff:
            line = line.strip().split("\t")
            frequencies[int(line[0])] = np.fromiter(line[1].split(","),
                                                    dtype=float)
    return frequencies


def sim_sites(treef: str, frequencies: np.array, model: str="JTT",
              gamma: tuple[int,float]=None, use_anc: int=None, rate: float=1.0,
              sites: int=1) -> dict:
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
    using a branch length scalar, or (for sites > 1) by using seq-gen directly,
    in which case arg gamma is a tuple of (ncat,alpha).
    """
    models = ["JTT", "WAG", "PAM", "BLOSUM", "MTREV", "CPREV45",
              "MTART", "LG", "GENERAL"]
    if model not in models:
        raise ValueError(f"{model} not in {models}")
    cmd = [
            "seq-gen",
            "-l",
            f"{sites}",
            "-m"
            f"{model}",
            "-f",
            f"{' '.join([str(x) for x in frequencies])}"
            "-of",
            "-wa"
        ]
    if use_anc is not None:
        cmd += ["-k", f"{use_anc}"]
    if rate != 1.0:
        if gamma is not None:
            raise ValueError("only one of rate or gamma allowed")
        cmd += ["-s", f"{rate}"]
    if gamma is not None:
        if rate != 1.0:
            raise ValueError("only one of rate or gamma allowed")
        cmd += ["-g", f"{gamma[0]}", "-a", f"{gamma[1]}"]
    cmd += [treef]
    print(subprocess.list2cmdline(cmd))
    res = subprocess.run(cmd, shell=False, capture_output=True, text=True,
                         check=False)
    out = dict(parse_phylip_str(res.stdout))
    return out


def concat_seqs(seq_dict1: dict, seq_dict2: dict) -> dict:
    """
    concatenate contents of two dictionaries with matching keys into a
    third
    """
    out = {}
    for header, seq in seq_dict1.items():
        try:
            out[header] = seq + seq_dict2[header]
        except KeyError:
            sys.stderr.write(f"{header} not in {seq_dict2}\n")
            sys.exit()
    return out


def write_seq_gen_input(tree: Node, seq_dict:dict=None):
    """
    write input tree (and optionally, sequence) file for seq-gen
    """
    with open("seq_gen_input.txt", "w", encoding="utf-8") as outf:
        if seq_dict is not None:
            outf.write(f"{len(seq_dict)} {len(seq_dict.values()[0])}\n")
            for header, seq in seq.items():
                outf.write(f"{header} {seq}\n")
            outf.write("1\n")
        outf.write(f"{tree.to_string()};")


def sim_alignment_gene_freqs(tree: Node, seq_dict: dict, model: str="JTT",
                             gamma: tuple[int,float]=None,
                             use_anc: int=None) -> dict:
    """
    simulates an alignment with length matching the input alignment and 
    optionally starting from a specified sequence as ancestor. Returns a
    sequence dictionary {header: sequence}
    """
    sites = len(seq_dict.values[0])
    freqs = calc_gene_frequencies(seq_dict)
    if use_anc is not None:
        write_seq_gen_input(tree, seq_dict)
        sim_aln = sim_sites(treef="seq_gen_input.txt", frequencies=freqs,
                            model=model, gamma=gamma, use_anc=use_anc,
                            sites=sites)
    else:
        write_seq_gen_input(tree)
        sim_aln = sim_sites(treef="seq_gen_input.txt", frequencies=freqs,
                            model=model, gamma=gamma, sites=sites)
    return sim_aln


def sim_alignment_site_freqs(tree: Node, seq_dict: dict, model: str="JTT",
                             gamma: tuple[int,float]=None, use_anc: int=None) -> dict:
    """
    simulates an alignment with length matching the input alignment and optionally starting from a
    specified sequence as ancestor, using site-specific frequencies for each site. Returns a
    sequence dictionary {header: sequence}
    """
    columns = get_columns(seq_dict)
    site_freqs = get_site_specific_frequencies(columns)
    ancs = {n.label: "" for n in tree.iternodes() if not n.istip}
    out = ancs | {header: "" for header in seq_dict}
    for pos, site in columns.items():
        if use_anc is not None:
            write_seq_gen_input(tree, site)
            sim = sim_sites(treef="seq_gen_input.txt", frequencies=site_freqs[pos],
                            model=model, gamma=gamma, use_anc=use_anc)
            
            out = concat_seqs(out, sim)
        
        