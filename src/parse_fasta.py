def parse_fasta(path):  # courtesy of Jonathan Chang https://gist.github.com/jonchang/6471846
    """Given a path tries to parse a fasta file. Returns an iterator which
    yields a (name, sequence) tuple"""
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
    """Given a str representing the content of a FASTA-formatted file, return
    an iterator yielding a (name, sequence tuple)"""
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
    """writes sequence dictionary to multiline string"""
    ret = ""
    for header, seq in seq_dict.items():
        ret += f">{header}\n{seq}\n"
    return ret


def write_fasta(seq_dict: dict, out_file: str):
    """writes a fasta file from sequence dictionary"""
    with open(out_file, "w", encoding="utf-8") as outf:
        outf.write(get_fasta_str(seq_dict=seq_dict))


if __name__ == "__main__":
    s = """
        >1
        A
        >2
        A
        >3
        A
        >4
        T
        """
    seqs = dict(parse_fasta_str(s))
    print(get_fasta_str(seqs))