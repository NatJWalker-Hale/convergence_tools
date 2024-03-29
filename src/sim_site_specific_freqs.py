#!/usr/bin/env python3


"""given an input file of site-specific frequencies, use seq-gen to repeatedly
simulate single sites and concatenate"""


import sys
import argparse
import subprocess
from copy import deepcopy
import numpy as np
import sequence as sq
import newick as nwk
from phylo import Node


def sim_sites(treef: str, frequencies: np.array, model: str="JTT",
              gamma: tuple[int,float]=None, use_anc: int=None, rate: float=1.0,
              sites: int=1) -> dict:
    """
    Simulates a single site with seq-gen given a treefile, which may optionally include sequences
    to serve as ancestral state, in which case treef should be a file with the following formatting:

    2 4
    taxon1 ATGC
    taxon2 ATGC
    1
    (taxon1:0.1,taxon2:0.1);

    rate information (e.g. from PAML posterior mean rates) can be incorporated using a branch length
    scalar, or (for sites > 1) by using seq-gen directly, in which case arg gamma is a tuple of 
    (ncat,alpha).
    """
    models = ["JTT", "WAG", "PAM", "BLOSUM", "MTREV", "CPREV45",
              "MTART", "LG", "GENERAL"]
    if model not in models:
        raise ValueError(f"{model} not in {models}")
    cmd = [
            "seq-gen",
            "-l",
            f"{sites}",
            "-m",
            f"{model}",
            "-f"
    ]
    cmd += [f"{s:.4f}" for s in frequencies]
    cmd += ["-of", "-wa"]
    # cmd += ["-z", "12345"]
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
    sys.stderr.write(f"{subprocess.list2cmdline(cmd)}\n")
    res = subprocess.run(cmd, shell=False, capture_output=True, text=True,
                         check=False)
    out = dict(sq.parse_phylip_str(res.stdout))
    return out


def sim_sites_indelible(tree: Node, seq_dict: dict, model: str="JTT", rates: dict=None,
                        site_freqs: dict=None, reps: int=100) -> dict[dict]:
    """
    mimics an alignment with site-specific frequencies (and optionally, rates) with indelible. Site
    rates are assumed to have mean 1.0 and are implemented as a branch length scalar. Returns a
    dictionary where key is replicate and value is seq_dict of replicate
    """
    extants = {n.label: seq_dict[n.label] for n in tree.iternodes() if n.istip}
    if site_freqs is None:
        site_freqs = sq.get_site_specific_frequencies(extants)
    with open("control.txt", "w", encoding="utf-8") as inf:
        inf.write("[TYPE] AMINOACID 1\n\n")
        inf.write("[SETTINGS]\n")
        inf.write("  [ancestralprint] SAME\n")
        inf.write("  [output] FASTA\n")
        inf.write("  [fastaextension] fa\n")
        inf.write("  [fileperrep] TRUE\n\n")
        for pos in site_freqs:  # iterate over sites
            inf.write(f"[MODEL] m{pos}\n")
            inf.write(f"  [submodel] {model}\n")
            inf.write("  [indelrate] 0\n")
            inf.write("  [rates] 0 0 0\n")
            inf.write("  [statefreq] " + " ".join([f"{s:.4f}" for s in site_freqs[pos]]) + "\n\n")
        for pos in site_freqs:
            inf.write(f"[TREE] t{pos} ")
            if rates is not None:
                scaled = deepcopy(tree)
                for n in scaled.iternodes():
                    n.length *= rates[pos]
                inf.write(f"{scaled.to_string(label=False)};\n")
            else:
                inf.write(f"{tree.to_string(label=False)};\n")
        inf.write("\n[PARTITIONS] sites\n")
        for pos in site_freqs:
            inf.write(f"  [t{pos}  m{pos}  1]\n")
        inf.write("[EVOLVE]\n")
        inf.write(f"  sites {reps} indelible_out\n")
    cmd = ["indelible", "control.txt"]
    sys.stderr.write(f"{subprocess.list2cmdline(cmd)}\n")
    subprocess.run(cmd, shell=False, check=False)
    reps_dict = {i: {} for i in range(reps)}
    for i in reps_dict:
        sim_rep = dict(sq.parse_fasta(f"indelible_out_TRUE_{i+1}.fa"))
        sim_rep = rename_ancestral_seqs_indelible(tree, sim_rep)
        sim_rep = sq.insert_gaps_by_seq(seq_dict, sim_rep)
        sq.write_fasta(sim_rep, f"indelible_out_TRUE_{i+1}_rn_gap.fa")
        reps_dict[i] = sim_rep

    return reps_dict


def concat_seqs(seq_dict1: dict, seq_dict2: dict) -> dict:
    """
    concatenate contents of two dictionaries with matching keys into a third
    """
    out = {}
    for header, seq in seq_dict1.items():
        try:
            out[header] = seq + seq_dict2[header]
        except KeyError as e:
            raise KeyError(f"{header} not in {seq_dict2}") from e
    return out


def write_seq_gen_input(tree: Node, seq_dict:dict=None):
    """
    write input tree (and optionally, sequence) file for seq-gen
    """
    with open("seq_gen_input.txt", "w", encoding="utf-8") as outf:
        if seq_dict is not None:
            outf.write(sq.get_phylip_str(seq_dict))
            outf.write("1\n")
        # note that we strip node labels in the seq-gen input because it causes the program to bug
        outf.write(f"{tree.to_string(label=False)};")


def rename_ancestral_seqs(tree: Node, seq_dict: dict) -> dict:
    """
    takes a seq_dict output including ancestral sequences as written directly from seq-gen and
    parsed, and a tree, and renames sequences in the seq_dict that are not tips in the tree 
    according to internal node labels of the tree
    """
    tips = len(tree.lvsnms())
    name = tips + 1
    out = {}
    for n in tree.iternodes(order=0):
        if not n.istip:
            out[n.label] = seq_dict[str(name)]
            name += 1
        else:
            out[n.label] = seq_dict[n.label]
    return out


def rename_ancestral_seqs_indelible(tree: Node, seq_dict: dict):
    """
    rename ancestral sequences from indelible
    """
    tips = len(tree.lvsnms())
    name = tips + 1
    out = {}
    for n in tree.iternodes():
        if not n.istip:
            if n.parent is None:
                out[n.label] = seq_dict["ROOT"]
            else:
                out[n.label] = seq_dict[f"N{name}"]
                name += 1
        else:
            out[n.label] = seq_dict[n.label]
    return out


def replace_gaps_by_random(frequencies: np.array=np.array([1.]*20)/20) -> str:
    """
    returns a random amino acid character in proportion to the probabilities specified in 
    frequencies
    """
    return np.random.choice(np.array(list("ARNDCQEGHILKMFPSTWYV"), dtype="str"),
                            size=1, p=frequencies)[0]


def sim_alignment_gene_freqs(tree: Node, seq_dict: dict, model: str="JTT",
                             gamma: tuple[int,float]=None,
                             use_anc: int=None) -> dict:
    """
    simulates an alignment with length matching the input alignment and optionally starting from a 
    specified sequence as ancestor. Returns a sequence dictionary {header: sequence}
    """
    sites = [len(v) for v in seq_dict.values()][0]
    freqs = sq.calc_gene_frequencies(seq_dict)
    if use_anc is not None:
        single_seq = {"N1": ""}
        for i in range(len(seq_dict["N1"])):
            if seq_dict["N1"][i] == "-":
                single_seq["N1"] += replace_gaps_by_random(freqs)
            else:
                single_seq["N1"] += seq_dict["N1"][i]
        write_seq_gen_input(tree, single_seq)
    else:
        write_seq_gen_input(tree)
    sim_aln = sim_sites(treef="seq_gen_input.txt", frequencies=freqs,
                        model=model, gamma=gamma, use_anc=use_anc,
                        sites=sites)
    sim_rn = rename_ancestral_seqs(tree, sim_aln)
    return sim_rn


def sim_alignment_site_freqs(tree: Node, seq_dict: dict, model: str="JTT", use_anc: int=None,
                             rates: dict=None) -> dict:
    """
    simulates an alignment with length matching the input alignment and optionally starting from a
    specified sequence as ancestor, using site-specific frequencies for each site. Returns a
    sequence dictionary {header: sequence}
    """
    extants = {n.label: seq_dict[n.label] for n in tree.iternodes() if n.istip}
    site_freqs = sq.get_site_specific_frequencies(extants)
    all_columns = sq.get_columns(seq_dict)
    out = {header: "" for header in seq_dict}
    for pos, site in all_columns.items():
        if use_anc is not None:
            single_site = {}
            if site["N1"] == "-":
                single_site["N1"] = replace_gaps_by_random(site_freqs[pos])
            else:
                single_site["N1"] = site["N1"]
            write_seq_gen_input(tree, single_site)
        else:
            write_seq_gen_input(tree)
        if rates is not None:
            rate = rates[pos]
        else:
            rate = 1.
        sim = sim_sites(treef="seq_gen_input.txt", frequencies=site_freqs[pos],
                        model=model, use_anc=use_anc, rate=rate)
        sim_rn = rename_ancestral_seqs(tree, sim)
        out = concat_seqs(out, sim_rn)
    return out


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    parser = argparse.ArgumentParser()
    parser.add_argument("aln", help="FASTA-formatted alignment to mimic with simulations")
    parser.add_argument("tree", help="tree to simulate sequences along, with node labels")
    parser.add_argument("-r", "--rates", help="PAML 'rates' file containing BEB estimates of \
                        site-specific rates to use per-site")
    parser.add_argument("-g", "--gapped", help="clone gaps from input alignment into simulated",
                        action="store_true")
    args = parser.parse_args()

    curroot = nwk.parse_from_file(args.tree)
    seqs = dict(sq.parse_fasta(args.aln))
    if args.rates:
        rates_dict = sq.parse_paml_rates(args.rates)
    else:
        rates_dict = None
    # sim_gf = sim_alignment_gene_freqs(tree=curroot, seq_dict=seqs, use_anc=1, gamma=(4, 0.8))
    # print(seqs)
    # if args.gapped:
    #     sim_gf = sq.insert_gaps_by_seq(seqs, sim_gf)
    # print(sq.get_fasta_str(sim_gf))

    # sim_sf = sim_alignment_site_freqs(tree=curroot, seq_dict=seqs, use_anc=1, rates=rates_dict)
    # if args.gapped:
    #     sim_sf = sq.insert_gaps_by_seq(seqs, sim_sf)
    # print(sq.get_fasta_str(sim_sf))

    # sim_sites_indelible(curroot, seqs, model="JTT", rates=rates_dict, reps=100)