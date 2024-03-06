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
    for line in fa_str.readlines():
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
        ret += f"{header}\n{seq}\n"
    return ret
