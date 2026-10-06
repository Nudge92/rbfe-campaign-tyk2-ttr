#!/usr/bin/env python3
"""Recompute the prevalidation metrics from the raw tables in ../data/.

Nothing is hard-coded: every number printed here comes out of the CSVs. Run it
to check the figures in the README against their sources.

ROC-AUC  (Mann-Whitney U with ties counted as 0.5)
  docking            data/docking_scores.csv     53 actives / 2,640 decoys
                     E_best and the ligand-efficiency column are better when
                     more negative, so their sign is flipped before ranking.
                     The tafamidis row is labelled -1 and excluded.
  fingerprints       data/similarity_all_fp_v2.csv  53 actives / 2,650 decoys
                     (10 decoys have no 3D conformer and so are missing from
                     the docking table; that is the only difference in N.)

Pair accuracy -- "given two molecules on the same scaffold, does the axis put
the more potent one first?" Reported for four separate tests that are often
quoted as one range; they do not share a compound set or an averaging rule.
  benzophenone  8 pairs   data/benzophenone_pairs.csv  + data/bp_summary.csv
  diphenylether 21 pairs  data/diphenylether_pairs.csv
  scaffold-wide 43 pairs  data/scaffold_pairs.csv (micro and macro both shown;
                          filtered to dpchembl >= 1.0 and |dheavy| <= 5)

Writes: ../data/prevalidation_summary.csv  (used by ../../figures/make_figures.py)
"""
import collections
import csv
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(os.path.dirname(HERE), "data")


def roc_auc(pos, neg):
    """Rank-sum AUC; tied scores contribute 0.5."""
    merged = sorted([(v, 1) for v in pos] + [(v, 0) for v in neg])
    rank_sum, i, n = 0.0, 0, len(merged)
    while i < n:
        j = i
        while j < n and merged[j][0] == merged[i][0]:
            j += 1
        rank = (i + 1 + j) / 2.0
        rank_sum += rank * sum(1 for k in range(i, j) if merged[k][1] == 1)
        i = j
    return (rank_sum - len(pos) * (len(pos) + 1) / 2.0) / (len(pos) * len(neg))


def auc_from(path, score_col, sign=1.0, label_col="label"):
    pos, neg = [], []
    with open(os.path.join(DATA, path), encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r[label_col] not in ("0", "1"):
                continue                       # -1 = reference drug, not scored
            (pos if r[label_col] == "1" else neg).append(sign * float(r[score_col]))
    return roc_auc(pos, neg), len(pos), len(neg)


def pair_accuracy(path, axis_col=None, correct_col="correct"):
    with open(os.path.join(DATA, path), encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    if axis_col is None:
        return {"(single axis)": (sum(float(r[correct_col]) for r in rows), len(rows))}
    out = {}
    for r in rows:
        hit, tot = out.get(r[axis_col], (0.0, 0))
        out[r[axis_col]] = (hit + float(r[correct_col]), tot + 1)
    return out


def scaffold_pairs(axis, dpchembl_min=1.0, dheavy_max=5.0):
    with open(os.path.join(DATA, "scaffold_pairs.csv"), encoding="utf-8") as fh:
        rows = [r for r in csv.DictReader(fh)
                if r["axis"] == axis
                and float(r["dpchembl"]) >= dpchembl_min
                and abs(float(r["dheavy"])) <= dheavy_max]
    by_scaffold = collections.defaultdict(list)
    for r in rows:
        by_scaffold[r["scaffold_smiles"]].append(float(r["correct"]))
    micro = sum(float(r["correct"]) for r in rows) / len(rows)
    macro = sum(sum(v) / len(v) for v in by_scaffold.values()) / len(by_scaffold)
    return micro, macro, len(rows), len(by_scaffold)


def main():
    out_rows = []

    print("== ROC-AUC (random = 0.5) ==")
    axes = [
        ("AtomPair fingerprint", "similarity_all_fp_v2.csv", "sim_AtomPair", 1.0),
        ("TopoTorsion fingerprint", "similarity_all_fp_v2.csv", "sim_TopoTorsion", 1.0),
        ("MACCS fingerprint", "similarity_all_fp_v2.csv", "sim_MACCS", 1.0),
        ("ECFP4 fingerprint", "similarity_all_fp_v2.csv", "sim_ECFP4", 1.0),
        ("Docking ligand efficiency", "docking_scores.csv", "LE", -1.0),
        ("Docking score (E_best)", "docking_scores.csv", "E_best", -1.0),
    ]
    for name, path, col, sign in axes:
        auc, npos, nneg = auc_from(path, col, sign)
        print(f"  {name:28s} {auc:.4f}   actives {npos} / decoys {nneg}")
        out_rows.append(dict(metric="roc_auc", name=name, value=round(auc, 4),
                             n_actives=npos, n_decoys=nneg, source=path))

    print("\n== pair accuracy (same scaffold, more potent first) ==")
    bp = pair_accuracy("benzophenone_pairs.csv")["(single axis)"]
    print(f"  benzophenone  MM/GBSA        {bp[0] / bp[1]:.3f}  ({bp[0]:.0f}/{bp[1]})")
    out_rows.append(dict(metric="pair_accuracy", name="benzophenone MM/GBSA (8 pairs)",
                         value=round(bp[0] / bp[1], 4), n_actives="", n_decoys=bp[1],
                         source="benzophenone_pairs.csv"))

    for axis, (hit, tot) in sorted(pair_accuracy("diphenylether_pairs.csv", "axis").items()):
        print(f"  diphenylether {axis:22s} {hit / tot:.3f}  ({hit:.1f}/{tot})")
        out_rows.append(dict(metric="pair_accuracy", name=f"diphenylether {axis} (21 pairs)",
                             value=round(hit / tot, 4), n_actives="", n_decoys=tot,
                             source="diphenylether_pairs.csv"))

    for axis in ("Borda(1:1) 53내순위 [LOO]", "AtomPair 단독 [LOO]",
                 "AtomPair 단독 [LSO]", "중원자 단독 [LOO]", "MW 단독 [LOO]"):
        micro, macro, npair, nscaf = scaffold_pairs(axis)
        print(f"  scaffold-wide {axis:26s} micro {micro:.3f}  macro {macro:.3f}"
              f"  ({npair} pairs / {nscaf} scaffolds)")
        for kind, val in (("micro", micro), ("macro", macro)):
            out_rows.append(dict(metric="pair_accuracy",
                                 name=f"scaffold-wide {axis} {kind} ({npair} pairs)",
                                 value=round(val, 4), n_actives="", n_decoys=npair,
                                 source="scaffold_pairs.csv"))

    dest = os.path.join(DATA, "prevalidation_summary.csv")
    with open(dest, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out_rows[0]))
        w.writeheader()
        w.writerows(out_rows)
    print(f"\nwrote {dest}")


if __name__ == "__main__":
    main()
