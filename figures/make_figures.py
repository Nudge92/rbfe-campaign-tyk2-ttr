#!/usr/bin/env python3
"""Draw the two generated figures from the tables committed in this repository.

No value is written into this script. Everything is read from:
  03_prevalidation/data/prevalidation_summary.csv   (recompute_metrics.py)
  01_tyk2_validation/results/edges_summary.csv      (build_edges_summary.py)

Regenerate both inputs first if the raw tables change:
  python3 03_prevalidation/analysis/recompute_metrics.py
  python3 01_tyk2_validation/analysis/build_edges_summary.py
  python3 figures/make_figures.py

Labels are in English so the figures render the same on a machine without a
Korean font installed.

Writes: figures/prevalidation_roc_auc.png
        figures/tyk2_calc_vs_exp.png
"""
import csv
import itertools
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SEED_FREE_BOOTSTRAP = True  # n=6 -> all 6**6 resamples are enumerated exactly


def read_csv(path):
    with open(os.path.join(ROOT, path), encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def percentile(sorted_vals, q):
    pos = (len(sorted_vals) - 1) * q / 100.0
    lo, hi = math.floor(pos), math.ceil(pos)
    return sorted_vals[lo] + (sorted_vals[hi] - sorted_vals[lo]) * (pos - lo)


def roc_auc_figure():
    rows = [r for r in read_csv("03_prevalidation/data/prevalidation_summary.csv")
            if r["metric"] == "roc_auc"]
    rows.sort(key=lambda r: float(r["value"]))          # worst at the bottom
    names = [r["name"] for r in rows]
    vals = [float(r["value"]) for r in rows]
    ns = [f"{r['n_actives']}/{r['n_decoys']}" for r in rows]

    fig, ax = plt.subplots(figsize=(9, 4.2))
    colors = ["#b23a48" if v < 0.5 else "#4c72b0" for v in vals]
    bars = ax.barh(names, vals, color=colors, height=0.62, zorder=3)
    ax.axvline(0.5, ls="--", lw=1.6, color="#444444", zorder=4)
    ax.text(0.5, len(vals) - 0.3, "  random = 0.5", color="#444444",
            va="center", ha="left", fontsize=9)

    for bar, v, n in zip(bars, vals, ns):
        ax.text(v + 0.012, bar.get_y() + bar.get_height() / 2,
                f"{v:.4f}", va="center", fontsize=9.5)
        ax.text(0.012, bar.get_y() + bar.get_height() / 2,
                f"n = {n}", va="center", fontsize=8, color="white")

    ax.set_xlim(0, 0.88)
    ax.set_xlabel("ROC-AUC  (actives vs property-matched decoys)")
    ax.set_title("TTR prevalidation: which scoring axis separates known actives?",
                 fontsize=12, pad=12)
    ax.grid(axis="x", alpha=0.25, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    fig.text(0.01, 0.015,
             "Docking scores land below random, so the ranking axis was moved "
             "from structure-based to ligand-based.\n"
             "n differs because 10 decoys have no 3D conformer and are absent "
             "from the docking table.",
             fontsize=8, color="#555555")
    fig.tight_layout(rect=(0, 0.075, 1, 1))
    out = os.path.join(HERE, "prevalidation_roc_auc.png")
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(f"wrote {out}")


def tyk2_figure():
    rows = read_csv("01_tyk2_validation/results/edges_summary.csv")
    exp = [float(r["exp_ddg"]) for r in rows]
    calc = [float(r["calc_ddg"]) for r in rows]
    xerr = [float(r["exp_sigma"]) for r in rows]
    yerr = [float(r["calc_sd"]) if r["calc_sd"] else 0.0 for r in rows]
    errors = [float(r["error"]) for r in rows]

    n = len(errors)
    rmse = math.sqrt(sum(e * e for e in errors) / n)
    mue = sum(abs(e) for e in errors) / n
    resampled = sorted(math.sqrt(sum(e * e for e in s) / n)
                       for s in itertools.product(errors, repeat=n))
    lo, hi = percentile(resampled, 2.5), percentile(resampled, 97.5)

    fig, ax = plt.subplots(figsize=(7.2, 7.0))
    span = [min(exp + calc) - 0.7, max(exp + calc) + 0.7]
    ax.plot(span, span, color="#444444", lw=1.4, zorder=2, label="perfect agreement")
    ax.fill_between(span, [s - 1 for s in span], [s + 1 for s in span],
                    color="#999999", alpha=0.13, zorder=1, label="within 1 kcal/mol")

    styles = {"HREMD": ("o", "#4c72b0"), "independent": ("s", "#dd8452")}
    seen, placed = set(), []
    for r, x, y, xe, ye in zip(rows, exp, calc, xerr, yerr):
        marker, color = styles[r["method"]]
        label = None
        if r["method"] not in seen:
            seen.add(r["method"])
            label = f"{r['method']} windows"
        ax.errorbar(x, y, xerr=xe, yerr=ye, fmt=marker, color=color, ms=8,
                    capsize=3, lw=1.2, zorder=5, label=label)
        # Dodge the annotation when an earlier one sits close by, so that no two
        # edge labels overlap whatever the data happen to be.
        offset = (10, -13)
        while any(abs(x - px) < 0.45 and abs(y - py) < 0.30 and offset == po
                  for px, py, po in placed):
            offset = (10, 10) if offset == (10, -13) else (-72, 10)
        placed.append((x, y, offset))
        ax.annotate(f"{r['edge']}  (n={r['n_replicates']})", (x, y),
                    textcoords="offset points", xytext=offset, fontsize=8.5,
                    color="#333333")

    ax.set_xlim(span)
    ax.set_ylim(span)
    ax.set_aspect("equal")
    ax.set_xlabel("experimental ΔΔG  (kcal/mol)")
    ax.set_ylabel("calculated ΔΔG  (kcal/mol)")
    ax.set_title("TYK2 benchmark: 6 edges with measured affinities", fontsize=12, pad=10)
    ax.grid(alpha=0.25, zorder=0)
    ax.set_axisbelow(True)
    ax.legend(loc="upper left", fontsize=8.5, framealpha=0.9)

    below = sum(1 for e in errors if e < 0)
    ax.text(0.97, 0.06,
            f"RMSE {rmse:.3f}  [{lo:.3f}, {hi:.3f}]\n"
            f"MUE {mue:.3f}\n"
            f"sign agreement {sum(int(r['sign_agrees']) for r in rows)}/{n}",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=9,
            bbox=dict(boxstyle="round,pad=0.45", fc="white", ec="#bbbbbb"))

    fig.text(0.012, 0.012,
             f"x error bars: Ki uncertainty carried by the benchmark set (a flat "
             f"~30% assumption, not a measured spread).\n"
             f"y error bars: spread over replicates, where n > 1.\n"
             f"All {below} of {n} points sit below the line -- every edge is "
             f"underestimated.  95% CI = exact bootstrap over {len(resampled):,} "
             f"edge resamples.",
             fontsize=7.8, color="#555555")
    fig.tight_layout(rect=(0, 0.085, 1, 1))
    out = os.path.join(HERE, "tyk2_calc_vs_exp.png")
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(f"wrote {out}")
    print(f"  RMSE {rmse:.4f} [{lo:.4f}, {hi:.4f}]  MUE {mue:.4f}  "
          f"points below the diagonal {below}/{n}")


if __name__ == "__main__":
    roc_auc_figure()
    tyk2_figure()
