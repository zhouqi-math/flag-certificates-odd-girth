#!/usr/bin/env sage
"""Exact SageMath checks for the two displayed certificates.

Run with sage, not plain Python. Without arguments the check writes no files.
Graphs use SageMath Graph; arithmetic uses ZZ, QQ, and QQ[u]. See README.md.
"""

import argparse
import csv
import json
from itertools import combinations, permutations
from pathlib import Path

from sage.all import Graph, PolynomialRing, QQ, ZZ
from sage.version import version as SAGE_VERSION

VERTICES = tuple(range(5))
PAIRS = tuple(combinations(VERTICES, 2))
ORDERS = tuple(permutations(VERTICES))
WEIGHTS = (QQ(3)/2, QQ(3)/8, QQ(1)/6, QQ(11)/24, QQ(1)/2)
# Entries not listed have coefficient zero. The types are E, P, P, P, P.
FLAGS = ({1: 1, 2: 1, 5: -3, 6: -3}, {0: 2, 2: -3, 4: -3},
         {1: 3, 2: -4, 4: -4, 6: 2}, {2: 1, 4: 1, 6: -2},
         {2: 1, 4: -1})

# Table 1 is an independent test oracle, not input to coefficient generation.
# Each row is (name, representative edge set, E/R/V/W/Z/Target).
TABLE = (
    ("5K1", (), (0, 0, 0, 0, 0, 0)),
    ("K2+3K1", ((0, 1),), (0, 0, 0, 0, 0, 0)),
    ("P3+2K1", ((0, 1), (1, 2)), (0, 16, 0, 0, 0, 6)),
    ("2K2+K1", ((0, 1), (2, 3)), (0, 0, 0, 0, 0, 0)),
    ("K1,3+K1", ((0, 1), (0, 2), (0, 3)), (12, 0, 0, 0, 0, 18)),
    ("P4+K1", ((0, 1), (1, 2), (2, 3)), (4, -48, 0, 0, 0, -12)),
    ("P3+K2", ((0, 1), (1, 2), (3, 4)), (0, 16, 0, 0, 0, 6)),
    ("C4+K1", ((0, 1), (1, 2), (2, 3), (0, 3)), (16, 0, 0, 0, 0, 24)),
    ("K1,4", ((0, 1), (0, 2), (0, 3), (0, 4)), (0, 0, 216, 0, 0, 36)),
    ("T", ((0, 1), (0, 2), (0, 3), (1, 4)), (-24, 36, -32, 4, 4, -24)),
    ("P5", ((0, 1), (1, 2), (2, 3), (3, 4)), (-24, -12, 64, 4, -4, -30)),
    ("U", ((0, 1), (1, 2), (2, 3), (0, 3), (0, 4)),
     (12, 0, -136, -16, 0, -12)),
    ("K2,3", ((0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4)),
     (0, 0, 192, 48, 0, 54)),
    ("C5", ((0, 1), (1, 2), (2, 3), (3, 4), (0, 4)),
     (180, 180, 320, 20, -20, 390)),
)


def require(condition, message):
    """Keep all checks active even when Python is run with -O."""
    if not condition:
        raise ValueError(message)


def graph_from_edges(edges):
    """The isolated vertices are retained in every five-vertex graph."""
    graph = Graph(5, loops=False, multiedges=False)
    graph.add_edges(edges)
    return graph


def triangle_free(graph):
    return not any(graph.has_edge(i, j) and graph.has_edge(i, k)
                   and graph.has_edge(j, k)
                   for i, j, k in combinations(VERTICES, 3))


def canonical(graph):
    """Sage's built-in algorithm; no optional bliss package is required."""
    return graph.canonical_label(algorithm="sage").graph6_string()


def degrees(a, subset):
    return sorted(sum(a[i, j] for j in subset) for i in subset)


def coefficients(graph):
    # Graph(5) has the vertices 0,...,4; fix this order explicitly.
    a = graph.adjacency_matrix(vertices=list(VERTICES))
    require(a.base_ring() is ZZ, "Adjacency matrix is not over ZZ")
    squares = [ZZ(0)] * 5
    for r0, r1, r2, v, w in ORDERS:
        root_edges = (a[r0, r1], a[r0, r2], a[r1, r2])
        if root_edges == (1, 0, 0):
            columns = (0,)  # E has only the edge 01.
        elif root_edges == (1, 1, 0):
            columns = (1, 2, 3, 4)  # P has precisely the edges 01 and 02.
        else:
            continue  # Wrong root type contributes zero, not renormalization.
        roots = (r0, r1, r2)
        j = sum(a[v, r] << i for i, r in enumerate(roots))
        k = sum(a[w, r] << i for i, r in enumerate(roots))
        for column in columns:
            squares[column] += FLAGS[column].get(j, 0) * FLAGS[column].get(k, 0)
    # Induced subset counts, computed separately from the flag expansion.
    p3 = QQ(sum(degrees(a, s) == [1, 1, 2]
                for s in combinations(VERTICES, 3))) / 10
    p4 = QQ(sum(degrees(a, s) == [1, 1, 2, 2]
                for s in combinations(VERTICES, 4))) / 5
    c5 = int(degrees(a, VERTICES) == [2, 2, 2, 2, 2])
    target = 120 * (p3 / 2 + 4 * c5 - p4)
    require(target.denominator() == 1, "Nonintegral scaled target")
    return tuple(int(value) for value in squares) + (int(target),)


def verify_flags():
    oracle = {}
    for name, edges, expected in TABLE:
        graph = graph_from_edges(edges)
        require(triangle_free(graph), "Invalid representative: " + name)
        key = canonical(graph)
        require(key not in oracle, "Repeated isomorphism class: " + name)
        oracle[key] = (name, edges, expected)
    groups = {}
    checked = 0
    for mask in range(1 << len(PAIRS)):
        graph = graph_from_edges(edge for bit, edge in enumerate(PAIRS)
                                 if mask >> bit & 1)
        if not triangle_free(graph):
            continue
        checked += 1
        key = canonical(graph)
        require(key in oracle, "An isomorphism class is missing from Table 1")
        name, _, expected = oracle[key]
        actual = coefficients(graph)
        require(actual == expected, "Table mismatch: " + name + ", mask=" + str(mask))
        weighted = sum(w * value for w, value in zip(WEIGHTS, actual[:5]))
        require(weighted == actual[5], "Flag identity failed: mask=" + str(mask))
        groups[key] = groups.get(key, 0) + 1
    require(checked == 388, "Wrong number of labelled triangle-free graphs")
    require(len(groups) == len(oracle) == 14, "Wrong number of isomorphism classes")
    maxima = tuple(max(map(abs, f.values())) for f in FLAGS)
    require(maxima == (3, 3, 4, 2, 1), "Wrong flag coefficient maxima")
    error_constant = sum(w * m * m for w, m in zip(WEIGHTS, maxima))
    require(error_constant == QQ(175)/8, "Wrong finite-error constant")
    rows = []
    for name, edges, expected in TABLE:
        actual = coefficients(graph_from_edges(edges))
        rows.append([name, ";".join(str(i) + "-" + str(j) for i, j in edges),
                     int(groups[canonical(graph_from_edges(edges))]), *actual,
                     str(sum(w * v for w, v in zip(WEIGHTS, actual[:5])))])
    return rows, {"all_labelled_graphs": int(2 ** len(PAIRS)),
                  "triangle_free_labelled_graphs": int(checked),
                  "isomorphism_classes": int(len(groups)),
                  "orders_per_graph": int(len(ORDERS)),
                  "table_checks_passed": int(checked),
                  "weighted_identity_checks_passed": int(checked),
                  "finite_error_constant": str(error_constant)}


def verify_cubic():
    ring = PolynomialRing(QQ, "u")
    u = ring.gen()
    r, eta, t, s = QQ(8)/105, QQ(89)/97, QQ(7)/5, QQ(18)/25
    c = 2 * t - t ** 2
    require(eta == (1 - 2 * r) / (1 - r), "Wrong spectral endpoint")
    require(c == QQ(21)/25, "Wrong test-vector coefficient")
    gap = 1 + t ** 2 - (1 + s) ** 2
    require(gap == QQ(1)/625 and gap > 0, "Invalid rational square-root bound")
    numerator = (u - c) * (eta ** 2 * u - 1)
    denominator = (u + s) ** 2 * (1 + eta) * (eta * (u + 1) - 2)
    cubic = 27590000*u**3 - 166092275*u**2 + 318800906*u - 189764295
    multiplier = ZZ(205821875)
    residual_clear = multiplier * (r * denominator - numerator) - cubic
    require(residual_clear == ring.zero(), "Clearing-denominators identity failed")
    lower, center, center2 = QQ(19)/16, QQ(97)/40, QQ(1173347)/482350
    constant = QQ(758113405787)/15435200
    decomposition = (27590000 * (u - lower) * (u - center) ** 2
                     + 482350 * (u - center2) ** 2 + constant)
    residual_square = decomposition - cubic
    require(residual_square == ring.zero(), "Square-decomposition identity failed")
    threshold = 1 / eta ** 2
    require(threshold == QQ(9409)/7921, "Wrong admissible-interval endpoint")
    require(0 < eta < 1 and threshold > lower and constant > 0,
            "The cubic positivity certificate has an invalid sign")
    require(threshold > 1 > c and threshold + s > 0
            and eta * (threshold + 1) - 2 > 0, "Invalid denominator sign")
    # Use a fixed coefficient range, including the zero cubic coefficient of
    # each residual; all rational values are exported as exact strings.
    return {"polynomial_ring": str(ring),
            "coefficient_order": "ascending powers of u",
            "cubic": [str(cubic[i]) for i in range(4)],
            "clear_denominators_multiplier": int(multiplier),
            "clear_denominators_residual": [str(residual_clear[i]) for i in range(4)],
            "square_decomposition_residual": [str(residual_square[i]) for i in range(4)],
            "admissible_endpoint": str(threshold), "interval_gap": str(threshold - lower),
            "radicand_gap": str(gap), "positive_remainder": str(constant)}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path,
                        help="write CSV and JSON results here after all checks pass")
    args = parser.parse_args(argv)
    rows, flag = verify_flags()
    cubic = verify_cubic()
    print("PASS: 1024 labelled graphs; 388 triangle-free graphs; 14 isomorphism classes.")
    print("PASS: Table 1 and the weighted identity on all 388 graphs (120 orders each).")
    print("PASS: finite-error constant = 175/8.")
    print("PASS: clearing-denominators and square-decomposition residuals are zero.")
    print("PASS: interval, denominator, radicand and positive-remainder sign checks.")
    if args.output_dir is not None:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        with (args.output_dir / "verified_coefficients.csv").open("w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(["graph", "representative_edges", "labelled_count",
                             "E", "R", "V", "W", "Z", "Target", "weighted_total"])
            writer.writerows(rows)
        result = {"status": "PASS", "software": "SageMath " + str(SAGE_VERSION),
                  "arithmetic": "SageMath ZZ, QQ, and PolynomialRing(QQ, u)",
                  "canonical_label_algorithm": "sage",
                  "flag_certificate": flag, "spectral_certificate": cubic}
        (args.output_dir / "verification_results.json").write_text(
            json.dumps(result, indent=2) + "\n", encoding="utf-8")
        print("WROTE: verified_coefficients.csv and verification_results.json.")


if __name__ == "__main__":
    main()
