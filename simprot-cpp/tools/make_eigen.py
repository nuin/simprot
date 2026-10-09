#!/usr/bin/env python3
"""Build the eigen data for an empirical amino acid model (matrix_data.hpp).

    python3 tools/make_eigen.py pam            # print C++ arrays for PAM
    python3 tools/make_eigen.py pam --check    # validate only
    python3 tools/make_eigen.py jtt --check    # same pipeline on JTT, as a control

Inputs are the exchangeabilities S_ij (lower triangle) and the equilibrium
frequencies pi, in the order ARNDCQEGHILKMFPSTWYV, from PAML's dat/ files
(prepared by Z. Yang):
  dayhoff.dat  Dayhoff, Schwartz & Orcutt (1978), Atlas of Protein Sequence
               and Structure 5 suppl. 3, 345-352
  jones.dat    Jones, Taylor & Thornton (1992), CABIOS 8:275-282

The rate matrix is Q_ij = S_ij * pi_j (i != j), scaled so that the expected
number of substitutions per site per unit of t is 0.01 (one PAM). The code
evaluates 100 * distance * rate for PAM and JTT, so a branch length is in
expected substitutions per site, as for PMB.

Layout, as for the JTT and PMB data: eigenvalues[k], and eigenvectors[k][i] =
U[i][k] * sqrt(pi_i) where U holds the orthonormal eigenvectors of the
symmetric matrix pi^1/2 Q pi^-1/2. Row 0 is the zero eigenvalue, so
eigenvectors[0][i] = pi_i, and
    P_ij(t) = sum_k eigenvectors[k][i] * eigenvectors[k][j] * exp(eigenvalues[k] * t) / pi_i
"""
import sys

import numpy as np

ORDER = "ARNDCQEGHILKMFPSTWYV"

DATA = {
    "pam": """
27
98 32
120 0 905
36 23 0 0
89 246 103 134 0
198 1 148 1153 0 716
240 9 139 125 11 28 81
23 240 535 86 28 606 43 10
65 64 77 24 44 18 61 0 7
41 15 34 0 0 73 11 7 44 257
26 464 318 71 0 153 83 27 26 46 18
72 90 1 0 0 114 30 17 0 336 527 243
18 14 14 0 0 0 0 15 48 196 157 0 92
250 103 42 13 19 153 51 34 94 12 32 33 17 11
409 154 495 95 161 56 79 234 35 24 17 96 62 46 245
371 26 229 66 16 53 34 30 22 192 33 136 104 13 78 550
0 201 23 0 0 0 0 0 27 0 46 0 0 76 0 75 0
24 8 95 0 96 0 22 0 127 37 28 13 0 698 0 34 42 61
208 24 15 18 49 35 37 54 44 889 175 10 258 12 48 30 157 0 28
0.087127 0.040904 0.040432 0.046872 0.033474 0.038255 0.049530
0.088612 0.033618 0.036886 0.085357 0.080482 0.014753 0.039772
0.050680 0.069577 0.058542 0.010494 0.029916 0.064718
""",
    "jtt": """
58
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
0.050901 0.068765 0.058565 0.014261 0.032102 0.066005
""",
}


def load(model):
    nums = [float(x) for x in DATA[model].split()]
    tri, pi = nums[:190], np.array(nums[190:210])
    S = np.zeros((20, 20))
    k = 0
    for i in range(1, 20):
        for j in range(i):
            S[i, j] = S[j, i] = tri[k]
            k += 1
    return S, pi / pi.sum()


def decompose(S, pi, rate=0.01):
    Q = S * pi[None, :]
    np.fill_diagonal(Q, 0.0)
    np.fill_diagonal(Q, -Q.sum(axis=1))
    Q *= rate / -np.sum(pi * np.diag(Q))
    sq = np.sqrt(pi)
    A = (sq[:, None] * Q) / sq[None, :]
    A = (A + A.T) / 2.0
    lam, U = np.linalg.eigh(A)
    order = np.argsort(-lam)                 # zero first, then slowest decay
    lam, U = lam[order], U[:, order]
    lam[0] = 0.0
    V = (U * sq[:, None]).T                  # V[k][i] = U[i][k] * sqrt(pi_i)
    if V[0].sum() < 0:
        V[0] = -V[0]
    return Q, lam, V


def transition(lam, V, pi, t):
    return np.einsum("ki,kj,k->ij", V, V, np.exp(lam * t)) / pi[:, None]


def check(model):
    S, pi = load(model)
    Q, lam, V = decompose(S, pi)
    ok = True

    def report(name, value, limit):
        nonlocal ok
        good = value <= limit
        ok &= good
        print(f"  {name:44s} {value:.2e}  {'ok' if good else 'FAIL'}")

    print(f"{model.upper()}: eigenvalues {lam[0]:.3g} .. {lam[-1]:.4g}, "
          f"{int(np.sum(lam > 1e-12))} positive")
    report("row 0 vs pi", np.abs(V[0] - pi).max(), 1e-12)
    report("metric orthonormality", np.abs((V / pi[None, :]) @ V.T - np.eye(20)).max(), 1e-12)
    rebuilt = np.einsum("ki,kj,k->ij", V, V, lam) / pi[:, None]
    report("Q rebuilt from eigen data", np.abs(rebuilt - Q).max(), 1e-14)
    report("substitutions per unit t - 0.01", abs(-np.sum(pi * np.diag(Q)) - 0.01), 1e-15)
    for t in (1.0, 10.0, 100.0, 1000.0):
        P = transition(lam, V, pi, t)
        report(f"P(t={t:g}) row sums - 1", np.abs(P.sum(axis=1) - 1).max(), 1e-12)
        report(f"P(t={t:g}) negative entries", max(0.0, -P.min()), 1e-15)
        report(f"P(t={t:g}) detailed balance", np.abs(pi[:, None] * P - (pi[:, None] * P).T).max(), 1e-15)
    P = transition(lam, V, pi, 1e6)
    report("P(t -> inf) rows vs pi", np.abs(P - pi[None, :]).max(), 1e-12)
    return ok, lam, V


def emit(model, lam, V):
    name = {"pam": "pam", "jtt": "jtt"}[model]
    out = [f"inline constexpr std::array<double, 20> {name}_eigenvalues = {{"]
    for k in range(0, 20, 4):
        out.append("    " + ", ".join(f"{x:.17g}" for x in lam[k:k + 4]) + ",")
    out.append("};")
    out.append("")
    out.append(f"inline constexpr std::array<std::array<double, 20>, 20> {name}_eigenvectors = {{{{")
    for k in range(20):
        out.append("    {")
        for i in range(0, 20, 4):
            out.append("        " + ", ".join(f"{x:.17g}" for x in V[k, i:i + 4]) + ",")
        out.append("    },")
    out.append("}};")
    return "\n".join(out)


def main(argv):
    if len(argv) < 2 or argv[1] not in DATA:
        sys.exit(__doc__)
    ok, lam, V = check(argv[1])
    if not ok:
        sys.exit(1)
    if "--check" not in argv:
        print(emit(argv[1], lam, V))


if __name__ == "__main__":
    main(sys.argv)
