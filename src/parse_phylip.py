#!/usr/bin/env python3


"""
Parses phylip-formatted sequence data
"""


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


if __name__ == "__main__":
    s = """
         4 1
        1   A
        2   T
        3   G
        4   C
        """
    print(dict(parse_phylip_str(s)))
