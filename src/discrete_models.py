#!/usr/bin/env python3


from collections import Counter
import numpy as np
from scipy.linalg import expm


class DiscreteModel():
    """
    Discrete state continuous time Markov chain model
    """
    def __init__(self):
        self.nstates = 0
        self.freqs: np.array
        self.mfreqs: np.array
        self.emfreqs: np.array
        self.eqfreqs: np.array
        self.alphabet = ""
        self.R: np.array
        self.Q: np.array

    def set_rate_JTT(self):
        """
        set rate matrix to JTT exchangeabilities
        """
        self.alphabet = "aa"
        self.nstates = 20
        modstr = """58
                54 45
                81 16 528
                56 113 34 10
                57 310 86 49 9
                105 29 58 767 5 323
                179 137 81 130 59 26 119
                27 328 391 112 69 597 26 23
                36 22 47 11 17 9 12 6 16
                30 38 12 7 23 72 9 6 56 229
                35 646 263 26 7 292 181 27 45 21 14
                54 44 30 15 31 43 18 14 33 479 388 65
                15 5 10 4 78 4 5 5 40 89 248 4 43
                194 74 15 15 14 164 18 24 115 10 102 21 16 17
                378 101 503 59 223 53 30 201 73 40 59 47 29 92 285
                475 64 232 38 42 51 32 33 46 245 25 103 226 12 118 477
                9 126 8 4 115 18 10 55 8 9 52 10 24 53 6 35 12
                11 20 70 46 209 24 7 8 573 32 24 8 18 536 10 63 21 71
                298 17 16 31 62 20 45 47 11 961 180 14 323 62 23 38 112 25 16
                0.076748 0.051691 0.042645 0.051544 0.019803 0.040752 0.061830
                0.073152 0.022944 0.053761 0.091904 0.058676 0.023826 0.040126
                0.050901 0.068765 0.058565 0.014261 0.032102 0.066005"""

        modstr = modstr.split("\n")

        ex = np.zeros(shape=(20, 20))

        for i in range(20):
            for j in range(20):
                if j == i:
                    continue
                if j > i:
                    ex[i, j] = float(modstr[j-1].strip().split(" ")[i])
                else:
                    ex[i, j] = float(modstr[i-1].strip().split(" ")[j])

        self.R = ex

        mfreqs = np.array(
                [float(j) for i in [x.strip().split(" ") for x in modstr[19:]]
                for j in i]
        )
        mfreqs = mfreqs / sum(mfreqs)
        self.mfreqs = mfreqs

    def set_model_freqs(self):
        """
        set frequencies to model frequencies
        """
        self.freqs = self.mfreqs

    def set_empirical_freqs(self):
        """
        set frequencies to empirical frequencies
        """
        self.freqs = self.emfreqs

    def set_equal_freqs(self):
        """
        set equal frequencies
        """
        self.eqfreqs = np.ones(self.nstates)
        self.eqfreqs = self.eqfreqs / sum(self.eqfreqs)

    def set_frequencies(self, freqs: np.array):
        """
        set frequencies to the values in freqs
        """
        self.freqs = freqs

    def scale_rate_matrix(self):
        """
        scale the rate matrix to an average rate of 1
        """
        big_pi = np.diag(self.freqs)

        q_unscaled = self.R @ big_pi

        for i in range(20):
            q_unscaled[i, i] = -sum(q_unscaled[i])

        q_ii = np.diag(q_unscaled)
        q = q_unscaled / -sum(q_ii * self.freqs)
        self.Q = q

    def get_P(self, brlen=0.1, rate=1.0):
        """
        calculate transition probability matrix by matrix exponentiation
        """
        p = expm(self.Q * brlen * rate)
        return p

    def calc_empirical_freqs(self, aln: dict):
        """
        takes a dictionary representing a multiple sequence alignment
        (key: header, value: sequence), and calculates state freqs
        """
        if len(set(len(v) for v in aln.values())) > 1:
            raise ValueError("sequences are not aligned")
        alph = set(check_alphabet(seq) for seq in aln.values())
        if len(alph) > 1:
            raise ValueError("multiple alphabets detected")
        alph = list(alph)[0]
        self.alphabet = alph
        states = get_states(alph)
        counts = Counter(char for seq in aln.values() for char in seq if
                         char in set(states))
        tot = sum(counts.values())
        frequencies = np.fromiter((counts[s] for s in states), dtype=float) / tot
        self.emfreqs = frequencies


def check_alphabet(sequence: str) -> str:
    """
    checks if an input sequence is nuc, aa, or other (returns mult)
    """
    aas = "ARNDCQEGHILKMFPSTWYV-.X"
    nucs = "ACGTURYSWKMBDHVN-."
    if len(set(sequence).union(set(nucs))) <= len(set(nucs)):
        return "nuc"
    if len(set(sequence).union(set(aas))) <= len(set(aas)):
        return "aa"
    return "mult"


def get_states(alphabet: str) -> str:
    """
    returns a string of states for a given alphabet
    """
    match alphabet:
        case "aa":
            return "ARNDCQEGHILKMFPSTWYV"
        case "nuc":
            return "ACGT"


if __name__ == "__main__":
    mod = DiscreteModel()
    mod.set_rate_JTT()
    mod.set_model_freqs()
    mod.scale_rate_matrix()

    # let's test that we recover the right likelihood with this using the following

    #             |--0.1--| A
    # |----0.2----|
    # |           |--0.1--| A
    # |
    # |--0.1--| N
    # |
    # |--0.1--| D

    p01 = mod.get_P(0.1)
    p02 = mod.get_P(0.2)
    # tip likelihoods are
    t1 = np.array([1.] + [0.]*19)
    t2 = np.array([1.] + [0.]*19)
    t3 = np.array([0., 0., 1.] + [0.]*17)
    t4 = np.array([0., 0., 0., 1.] + [0.]*16)
    node1_P = np.zeros(20)
    for i in range(20):
        node1_P[i] = sum(p01[i,] * t1) * sum(p01[i ,] * t1)

    root_P = np.zeros(20)
    for i in range(20):
        root_P[i] = sum(p02[i,] * node1_P) * sum(p01[i ,] * t3) * sum(p01[i ,] * t4)

    l = sum(mod.freqs * root_P)
    logl = np.log(l)
    assert round(logl, 6) == -11.045645

    # print(mod.get_P(0.1))

    # experiment with exactly as written J. Zhang paper
    # mod = Discrete_model()
    # mod.set_rate_JTT()
    # print(expm(logm(mod.R)*0.1))

    # seqs = dict(parse_fasta("DODAa_combined_no_og_strict_for_synth.cds.fa.nostop.name.noF.best.fas.trans"))
    # mod.calc_empirical_freqs(seqs)
    # print(mod.emfreqs)
