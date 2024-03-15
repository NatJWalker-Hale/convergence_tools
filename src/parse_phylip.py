#!/usr/bin/env python3


"""
Parses phylip-formatted sequence data
"""


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


def write_phylip_str(seq_dict) -> str:
    """
    writes a PHYLIP-formatted string from an aligned sequence dictionary {header: sequence}
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


if __name__ == "__main__":
    s = """
         4 1
        1   ATGC
        2   ATGC
        3   ATGC
        4   ATGC
        """
    seqs = dict(parse_phylip_str(s))
    print(seqs)
    print(write_phylip_str(seqs))
