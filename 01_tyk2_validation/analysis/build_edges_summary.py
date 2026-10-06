#!/usr/bin/env python3
"""Build results/edges_summary.csv for the 6 scored TYK2 edges.

Every number is read from a file committed in this repository -- nothing is
hard-coded except the mapping of each edge to its ligand pair and to the run
that was scored. That mapping is documented in ../TYK2_계산값_대장.md.

Calculated ddG (MBAR, kcal/mol), per replicate:
  edge1  HREMD, n=1          results/hr2_ddg.json      ["edge1"]["ddg"]["MBAR"]
  edge2  independent, n=3    results/step4_raw.json    ["edge2"][rep]["ddg"]["MBAR"]
  edge3  independent, n=1    results/step4_raw.json    ["edge3"][rep]
  edge4  independent, n=1    results/step4_raw.json    ["edge4"][rep]
  edge5  HREMD, n=3          results/hrB_ddg.json      [rep]["ddg"]["MBAR"]
  edge6  HREMD, n=3          results/hrA_ddg.json      [rep]["MBAR"]

Experimental ddG is derived from the Ki values in results/ligands_ki.yml
(doi 10.1016/j.ejmech.2013.03.070, Table 4) as RT*ln(Ki_B / Ki_A) at 300 K,
the simulation temperature. Its uncertainty is propagated from the error field
of the same file, sigma(ddG) = RT*sqrt((sA/KiA)^2 + (sB/KiB)^2), which gives
~0.25 kcal/mol per edge. Note that this error field is a flat 29-31% of the
value for all 13 ligands, so it is an assumption carried by the benchmark set
rather than a measured spread.

Writes: results/edges_summary.csv
"""
import csv
import json
import math
import os
import re
import statistics

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(os.path.dirname(HERE), "results")
RT = 0.001987 * 300.0  # kcal/mol

# edge -> (ligand A, ligand B, sampling method, source file, how to reach the MBAR values)
EDGES = [
    ("edge1", "ejm_31", "ejm_42", "HREMD",       "hr2_ddg.json",   "edge"),
    ("edge2", "ejm_46", "ejm_42", "independent", "step4_raw.json", "edge"),
    ("edge3", "ejm_46", "ejm_50", "independent", "step4_raw.json", "edge"),
    ("edge4", "ejm_46", "ejm_47", "independent", "step4_raw.json", "edge"),
    ("edge5", "ejm_31", "ejm_50", "HREMD",       "hrB_ddg.json",   "reps_nested"),
    ("edge6", "ejm_42", "ejm_50", "HREMD",       "hrA_ddg.json",   "reps_flat"),
]


def load(name):
    with open(os.path.join(RESULTS, name), encoding="utf-8") as fh:
        return json.load(fh)


def ki_nanomolar():
    """Read each ligand's Ki and its error; avoids a PyYAML dependency.

    The file mixes units (the ejm_* series is uM, the jmc_* series nM), so every
    value is converted to nM. Only ratios are used, so the common unit chosen
    does not affect the result. Returns {ligand: (Ki_nM, sigma_nM)}.
    """
    scale = {"uM": 1000.0, "nM": 1.0, "mM": 1e6}
    path = os.path.join(RESULTS, "ligands_ki.yml")
    out, lig, unit, value, err = {}, None, None, None, None

    def flush():
        if lig and value is not None:
            if unit not in scale:
                raise SystemExit(f"unexpected Ki unit for {lig}: {unit}")
            out[lig] = (value * scale[unit], (err or 0.0) * scale[unit])

    with open(path, encoding="utf-8") as fh:
        for line in fh:
            m = re.match(r"^lig_(\S+):\s*$", line)
            if m:
                flush()
                lig, unit, value, err = m.group(1), None, None, None
                continue
            m = re.match(r"^\s+unit:\s*(\S+)\s*$", line)
            if m:
                unit = m.group(1)
            m = re.match(r"^\s+value:\s*(\S+)\s*$", line)
            if m:
                value = float(m.group(1))
            m = re.match(r"^\s+error:\s*(\S+)\s*$", line)
            if m:
                err = float(m.group(1))
    flush()
    return out


def mbar_values(edge, how, fname):
    d = load(fname)
    if how == "edge":
        node = d[edge]
        if "ddg" in node:                       # hr2_ddg.json
            return [node["ddg"]["MBAR"][0]]
        return [node[r]["ddg"]["MBAR"][0]       # step4_raw.json
                for r in sorted(node, key=int)]
    if how == "reps_nested":                    # hrB_ddg.json
        return [d[r]["ddg"]["MBAR"][0] for r in sorted(d)]
    if how == "reps_flat":                      # hrA_ddg.json
        return [d[r]["MBAR"][0] for r in sorted(d)]
    raise ValueError(how)


def main():
    ki = ki_nanomolar()
    rows, errors = [], []
    for edge, a, b, method, fname, how in EDGES:
        reps = mbar_values(edge, how, fname)
        calc = sum(reps) / len(reps)
        (ki_a, s_a), (ki_b, s_b) = ki[a], ki[b]
        exp = RT * math.log(ki_b / ki_a)
        sigma_exp = RT * math.hypot(s_a / ki_a, s_b / ki_b)
        err = calc - exp
        errors.append(err)
        rows.append(dict(
            edge=edge, ligand_a=a, ligand_b=b, method=method, n_replicates=len(reps),
            calc_ddg=round(calc, 4),
            calc_sd=round(statistics.stdev(reps), 4) if len(reps) > 1 else "",
            exp_ddg=round(exp, 4), exp_sigma=round(sigma_exp, 4),
            error=round(err, 4),
            sign_agrees=int(calc * exp > 0),
            replicates=" ".join(f"{v:+.4f}" for v in reps),
            source=fname,
        ))

    out = os.path.join(RESULTS, "edges_summary.csv")
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    n = len(errors)
    rmse = math.sqrt(sum(e * e for e in errors) / n)
    mue = sum(abs(e) for e in errors) / n
    print(f"wrote {out}")
    print(f"n={n}  RMSE={rmse:.4f}  MUE={mue:.4f}  "
          f"sign agreement={sum(r['sign_agrees'] for r in rows)}/{n}  "
          f"errors all negative={all(e < 0 for e in errors)}")


if __name__ == "__main__":
    main()
