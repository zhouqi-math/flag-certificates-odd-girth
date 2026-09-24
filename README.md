# Exact SageMath certificate verifier

Supplement to **The least signless Laplacian eigenvalue of {C3,C5}-free graphs**, by Qi Zhou and Wei Xiang, revision 4.

This small program independently expands the displayed flags, checks every entry of Table 1, and verifies the exact rational identities at the end of the spectral proof. It is a checking aid for the certificates already given in the paper. It does not reproduce, or require reproducing, the numerical search that originally suggested them.

Repository: https://github.com/zhouqi-math/flag-certificates-odd-girth

The manuscript contains the mathematical proofs, flag definitions, full coefficient table, and rational polynomial certificate. This repository supplies a short independent check of the finite calculations.

## Requirements and use

Use **SageMath** with the `sage` command available. Plain Python is not sufficient. The verifier uses SageMath's native `Graph` class, integer ring `ZZ`, rational field `QQ`, and polynomial ring `PolynomialRing(QQ, "u")`, together with the Python standard-library modules `argparse`, `csv`, `json`, `itertools`, and `pathlib`. No optional Sage package, optimization solver, network access, or numerical-search software is required. All mathematical calculations are exact; there are no floating-point calculations or numerical tolerances. Canonical graph labelling explicitly selects SageMath's built-in `"sage"` algorithm, so the optional `bliss` package is not needed.

From the directory containing this README and `verify_certificates.sage`, run:

```text
sage verify_certificates.sage
```

This checks the certificates without writing files. To write the coefficient table and verification results, run:

```text
sage verify_certificates.sage --output-dir .
```

A successful run prints:

```text
PASS: 1024 labelled graphs; 388 triangle-free graphs; 14 isomorphism classes.
PASS: Table 1 and the weighted identity on all 388 graphs (120 orders each).
PASS: finite-error constant = 175/8.
PASS: clearing-denominators and square-decomposition residuals are zero.
PASS: interval, denominator, radicand and positive-remainder sign checks.
```

The second command also prints `WROTE: verified_coefficients.csv and verification_results.json.` A failed mathematical check raises `ValueError` and terminates with nonzero exit status. Checks use explicit conditions rather than `assert`, so Python optimization does not disable them. Output files are written only after all checks pass; the JSON also records the actual SageMath version. The verifier itself writes no result files unless `--output-dir` is supplied. SageMath may create its normal preparsed `.sage.py` file when launching a `.sage` script.

In a SageMath notebook or console, the same file can instead be executed with an explicit empty argument list, avoiding the notebook's own command-line arguments:

```python
from pathlib import Path
ns = {"__name__": "certificate_verifier"}
exec(Path("verify_certificates.sage").read_text(encoding="utf-8"), ns)
ns["main"]([])
```

The file uses Python-compatible syntax with explicit SageMath imports, so this form needs no preparsing. To save the results, replace the last line with `ns["main"](["--output-dir", "."])`.

## Validation status and reference data

This SageMath conversion has passed static review but has not been executed in SageMath in the editing environment, where SageMath is unavailable. No Python execution result is presented as a SageMath result. `reference_coefficients.csv` and `reference_values.json` contain the previously checked exact values for comparison. They are reference data, not files read by the verifier. A successful SageMath run with `--output-dir` creates fresh `verified_coefficients.csv` and `verification_results.json`, including the actual SageMath version.

## Graphs, labels, and averaging

The vertex set is `{0,1,2,3,4}`. The ten possible edges are listed in lexicographic order by `combinations(range(5), 2)`, and the program tests all `2^10 = 1024` edge masks. It discards graphs containing a triangle. Each edge mask is instantiated as a SageMath `Graph` with all five vertices explicitly retained. The remaining 388 labelled graphs are grouped into 14 isomorphism classes using `graph.canonical_label(algorithm="sage").graph6_string()`. The canonical graph string is used only as an isomorphism-class key; it is not a precomputed table of coefficients. SageMath integer adjacency matrices are used to calculate the flags and induced-subset densities.

For each of the 388 graphs, every permutation `(r0,r1,r2,v,w)` is examined. Its first three vertices carry the ordered root labels `0,1,2`; its last two are distinct extension vertices outside the roots. Type `E` has precisely the root edge `01`; type `P` has precisely the root edges `01,02`. The extension flag index is the bit mask

```text
j = adjacency(v,r0) + 2*adjacency(v,r1) + 4*adjacency(v,r2).
```

The code multiplies the coefficients of the two extension flags in the corresponding displayed linear combination `L_E`, `L_R`, `L_V`, `L_W`, or `L_Z`. A wrong root type contributes **zero**. The averaging is unconditional: there is no division by the number of roots of a particular type. The edge between the two extension vertices is unrestricted. Thus each square column is the sum over all 120 permutations, or equivalently **120 times the unlabelled averaged square coefficient**.

The target is calculated separately, by counting induced vertex subsets:

```text
120 * (p(P3,J)/2 + 4*p(C5,J) - p(P4,J)).
```

Here `p(P3,J)` uses the ten 3-subsets, and `p(P4,J)` uses the five 4-subsets. Sorted induced degree sequences `[1,1,2]` and `[1,1,2,2]` identify these two paths. A five-vertex simple graph in which all degrees are 2 is a pentagon. These induced-subset computations do not use the flag expansion or the expected table entries.

The expected Table 1 entries are an explicit test oracle. The program independently computes all six columns on **every labelled graph**, compares them with its class's expected row, and checks the weighted identity using weights `3/2, 3/8, 1/6, 11/24, 1/2`. The coefficient maxima `3,3,4,2,1` are obtained from the flag definitions, and their weighted squares give the finite-error constant `175/8`.

## Table representatives

The ASCII graph names below correspond to the graph names in Table 1. All representatives use vertices `0,1,2,3,4`, including isolated vertices. With `--output-dir`, the same edge lists and each class's labelled multiplicity are written to `verified_coefficients.csv`.

| Graph | Representative edges |
|---|---|
| `5K1` | none |
| `K2+3K1` | 01 |
| `P3+2K1` | 01, 12 |
| `2K2+K1` | 01, 23 |
| `K1,3+K1` | 01, 02, 03 |
| `P4+K1` | 01, 12, 23 |
| `P3+K2` | 01, 12, 34 |
| `C4+K1` | 01, 12, 23, 03 |
| `K1,4` | 01, 02, 03, 04 |
| `T` | 01, 02, 03, 14 |
| `P5` | 01, 12, 23, 34 |
| `U` | 01, 12, 23, 03, 04 |
| `K2,3` | 02, 03, 04, 12, 13, 14 |
| `C5` | 01, 12, 23, 34, 04 |

Here `T` is the tree with degree sequence `(3,2,1,1,1)`, and `U` is a four-cycle with one pendant edge.

## Final rational identities

Polynomial identities are calculated directly in SageMath's native exact polynomial ring `PolynomialRing(QQ, "u")`. The program constructs the displayed expressions and requires both residual polynomials to equal the zero polynomial. Coefficients are exported in ascending powers of `u` only for readable output; the arithmetic does not use a hand-written polynomial implementation. Put

```text
r0 = 8/105, eta = 89/97, t = 7/5, c = 21/25, s = 18/25,
N(u) = (u-c)(eta^2*u-1),
D(u) = (u+s)^2*(1+eta)*(eta*(u+1)-2).
```

It checks coefficient by coefficient that

```text
205821875 * (r0*D(u) - N(u))
  = 27590000*u^3 - 166092275*u^2 + 318800906*u - 189764295
  = 27590000*(u-19/16)*(u-97/40)^2
    + 482350*(u-1173347/482350)^2 + 758113405787/15435200.
```

It also checks the rational facts supporting positivity: `eta=(1-2*r0)/(1-r0)`, `c=2*t-t^2`, `1+t^2-(1+s)^2=1/625>0`, `1/eta^2=9409/7921>19/16`, positivity of the denominator factors on `u>1/eta^2`, and positivity of the constant remainder. The program requires both polynomial residuals to be exactly zero. With `--output-dir`, these checked facts are recorded in `verification_results.json` as exact integer or rational strings.

## Scope

The enumeration verifies the finite coefficient identity and the arithmetic constant used for the finite-error bound. The polynomial calculation verifies the displayed algebra and its sign conditions. The program is not a formal proof assistant: the manuscript supplies the blow-up limit, the trace and Rayleigh inequalities, and the deductions for graphs of arbitrary order. No sampling of large graphs, floating-point optimization, or reconstruction of the initial numerical search is part of this verification.
