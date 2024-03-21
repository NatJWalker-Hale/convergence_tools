#!/usr/bin/env python3


"""
utilties for parsing, writing and processing common sequence file formats
"""


from collections import Counter
import numpy as np


def check_aligned(seq_dict: dict) -> bool:
    """
    checks if sequences in a sequence dictionary are the same length
    """
    if len(set(len(s) for s in seq_dict.values())) > 1:
        return False
    return True


def parse_phylip(path: str):
    """
    parse relaxed PHYLIP from file. Will fail with interleaved
    """
    with open(path, "r", encoding="utf-8") as inf:
        phy_str = inf.read()
        parse_phylip_str(phy_str)


def parse_phylip_str(phy_str: str):
    """
    parse relaxed PHYLIP from string. Will fail with interleaved
    """
    lines = phy_str.strip().split("\n")
    for line in lines[1:]:
        line = line.strip().split()
        yield line[0], line[1]


def get_phylip_str(seq_dict) -> str:
    """
    writes a PHYLIP-formatted string from an aligned sequence dictionary {header: sequence}. First
    allows a particular sequence to be placed at the start
    """
    if not check_aligned(seq_dict):
        raise ValueError("sequences are not aligned, write to FASTA instead")
    nseq = len(seq_dict)
    seql = set(len(v) for v in seq_dict.values()).pop()
    out = ""
    out += f" {nseq} {seql}\n"
    for header, seq in seq_dict.items():
        out += f"{header}\t{seq}\n"
    return out


def parse_fasta(path):  # courtesy of Jonathan Chang https://gist.github.com/jonchang/6471846
    """
    given a path tries to parse a fasta file. Returns an iterator which yields a (name, sequence) 
    tuple
    """
    with open(path, "r", encoding="utf-8") as handle:
        name = sequence = ""
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    yield name, sequence
                name = line[1:]
                sequence = ""
                continue
            sequence += line
        # yield the last sequence
        if name and sequence:
            yield name, sequence


def parse_fasta_str(fa_str: str):
    """
    Given a str representing the content of a FASTA-formatted file, return an iterator yielding 
    a (name, sequence tuple)
    """
    name = sequence = ""
    for line in fa_str.splitlines():
        line = line.strip()
        if line.startswith(">"):
            if name:
                yield name, sequence
            name = line[1:]
            sequence = ""
            continue
        sequence += line
    if name and sequence:
        yield name, sequence


def get_fasta_str(seq_dict: dict):
    """
    writes sequence dictionary to multiline string
    """
    ret = ""
    for header, seq in seq_dict.items():
        ret += f">{header}\n{seq}\n"
    return ret


def write_fasta(seq_dict: dict, out_file: str):
    """
    writes a fasta file from sequence dictionary
    """
    with open(out_file, "w", encoding="utf-8") as outf:
        outf.write(get_fasta_str(seq_dict=seq_dict))


def get_columns(seq_dict: dict) -> dict:
    """
    takes a dictionary of an alignment (key: name, value: sequence) and returns a dictionary of 
    columns (key: position, value: dict{name: state})
    """
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


def get_site_specific_frequencies(seq_dict: dict, smooth: bool=True) -> dict:
    """
    takes a dictionary from get_columns and returns a dictionary of site-specific frequencies
    """
    aa = "ARNDCQEGHILKMFPSTWYV"
    col_dict = get_columns(seq_dict)
    freq_dict = {}  # key is pos, value is np array of freqs
    for pos, column in col_dict.items():
        counts = Counter(char for char in column.values() if char in set(aa))
        tot = sum(counts.values())
        if smooth:  # smoothing
            freqs = np.fromiter((counts[char] + 0.5 for char in aa), dtype=float) / (tot +
                                                                                        10)
            # Jeffrey's prior
        else:  # no smoothing
            freqs = np.fromiter((counts[char] for char in aa), dtype=float) / tot
        freq_dict[pos] = freqs
    return freq_dict


def calc_gene_frequencies(seq_dict: dict) -> np.array:
    """
    calculates the frequencies of amino acid character states in sequence dictionary, ignoring gaps
    """
    aa = "ARNDCQEGHILKMFPSTWYV"
    counts = Counter(char for seq in seq_dict.values() for char in seq if
                     char in set(aa))
    tot = sum(counts.values())
    frequencies = np.fromiter((counts[char] for char in aa), dtype=float) / tot
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


def parse_paml_rates(path: str) -> dict:
    """takes path to a paml 'rates' file containing estimated posterior mean 
    rate per site and returns a dictionary of {site: rate}"""
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
                out[int(line[0]) - 1] = float(line[-2])  # for 0-indexing
    return out


def insert_gaps_by_seq(ref:dict, query:dict):
    """clones gaps from a reference sequence into a query sequence, for example
    to insert gaps in simulated sequences to match an empirical alignment"""
    out = {}
    for header, seq in query.items():
        out[header] = ""
        try:
            for n, _ in enumerate(seq):
                if ref[header][n] == "-":
                    out[header] += "-"
                else:
                    out[header] += query[header][n]
        except KeyError as e:
            raise KeyError(header + " not in reference data, skipping\n") from e
    return out


if __name__ == "__main__":
    seqs = dict(parse_fasta("DODAa_combined_no_og_strict_for_synth.cds.fa.nostop.name.noF.best.fas.trans"))
    site_freqs = get_site_specific_frequencies(seqs)
    print(site_freqs)
