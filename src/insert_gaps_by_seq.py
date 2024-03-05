#!/usr/bin/env python3


import sys
import argparse
from parse_fasta import parse_fasta, get_fasta_str


def insert_gaps_by_seq(ref:dict, query:dict):
    """clones gaps from a reference sequence into a query sequence, for example
    to insert gaps in simulated sequences to match an empirical alignment"""
    out = {}
    for header, seq in query.items():
        try:
            for n, _ in enumerate(seq):
                if ref[header][n] == "-":
                    out[header]= query[header][:n] + "-" + query[header][n+1:]
        except KeyError:
            sys.stderr.write(header + " not in simulated data, skipping\n")
    return out


if __name__ == "__main__":
    if len(sys.argv[1:]) == 0:
        sys.argv.append("-h")

    ap = argparse.ArgumentParser()
    ap.add_argument("-a", "--alignment", help="Empirical alignment, \
                    in FASTA format.")
    ap.add_argument("-s", "--simulated", help="Simulated alignment, \
                    in FASTA format.")
    args = ap.parse_args()

    emp_seqs = dict(parse_fasta(args.alignment))
    sim_seqs = dict(parse_fasta(args.simulated))

    gaps_clone = insert_gaps_by_seq(emp_seqs, sim_seqs)

    print(get_fasta_str(gaps_clone))