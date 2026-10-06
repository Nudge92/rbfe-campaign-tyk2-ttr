#!/usr/bin/env python3
"""TYK2 6-edge RMSE/MUE with a bootstrap 95% CI (edge-level resampling).

Reproduces the interval reported for the TYK2 pilot (RMSE 0.793, 95% CI 0.434-1.053).
Standard library only. Two intervals are printed:

  1. Monte Carlo bootstrap: N_BOOT resamples of the 6 edges with replacement, SEED fixed.
  2. Exact bootstrap: all 6**6 = 46,656 equally likely resamples enumerated, so the
     result does not depend on a seed. With n = 6 this is cheap and is the reference.

Percentiles use linear interpolation between order statistics (numpy's default).

Inputs are per-replicate MBAR ddG values (kcal/mol) and experimental Ki (nM).
Sources in the TYK2 pilot tree (07_step3_lambda/):
  edge1            hr2_summary.txt   (HREMD, 1 run)
  edge2, 3, 4      step4_summary.txt (independent windows; 3, 1, 1 runs)
  edge5            hrB_ddg.txt       (HREMD, 3 runs)
  edge6            hrA_ddg.txt       (HREMD, 3 runs)
  Ki               01_inputs/00_data/ligands.yml (doi 10.1016/j.ejmech.2013.03.070)
"""
import itertools
import math
import random

SEED = 20261007
N_BOOT = 10_000
RT = 0.001987 * 300.0  # kcal/mol at 300 K, the simulation temperature

KI_NM = {"ejm_31": 96.0, "ejm_42": 64.0, "ejm_46": 4.8, "ejm_47": 74.0, "ejm_50": 250.0}

# edge: (ligand A, ligand B, per-replicate calculated ddG A->B)
EDGES = {
    "edge1": ("ejm_31", "ejm_42", [-0.4705]),
    "edge2": ("ejm_46", "ejm_42", [0.6724, 0.8163, 0.6175]),
    "edge3": ("ejm_46", "ejm_50", [1.1038]),
    "edge4": ("ejm_46", "ejm_47", [0.5437]),
    "edge5": ("ejm_31", "ejm_50", [0.4035, 0.5117, 0.1057]),
    "edge6": ("ejm_42", "ejm_50", [0.1714, -0.4939, 1.3842]),
}


def rmse(errs):
    return math.sqrt(sum(e * e for e in errs) / len(errs))


def mue(errs):
    return sum(abs(e) for e in errs) / len(errs)


def percentile(sorted_vals, q):
    """q in [0, 100]; linear interpolation between order statistics."""
    pos = (len(sorted_vals) - 1) * q / 100.0
    lo = math.floor(pos)
    hi = math.ceil(pos)
    return sorted_vals[lo] + (sorted_vals[hi] - sorted_vals[lo]) * (pos - lo)


def ci95(samples):
    s = sorted(samples)
    return percentile(s, 2.5), percentile(s, 97.5)


def main():
    errors = []
    print(f"{'edge':6s} {'n':>2s} {'calc':>8s} {'exp':>8s} {'error':>8s}")
    for name, (a, b, reps) in EDGES.items():
        calc = sum(reps) / len(reps)
        exp = RT * math.log(KI_NM[b] / KI_NM[a])
        errors.append(calc - exp)
        print(f"{name:6s} {len(reps):2d} {calc:+8.4f} {exp:+8.4f} {calc - exp:+8.4f}")

    n = len(errors)
    agree = sum(1 for (a, b, reps) in EDGES.values()
                if (sum(reps) / len(reps)) * math.log(KI_NM[b] / KI_NM[a]) > 0)
    print(f"\nRMSE {rmse(errors):.4f}   MUE {mue(errors):.4f}   sign agreement {agree}/{n}")

    rng = random.Random(SEED)
    mc = [[rng.choice(errors) for _ in range(n)] for _ in range(N_BOOT)]
    lo, hi = ci95([rmse(s) for s in mc])
    mlo, mhi = ci95([mue(s) for s in mc])
    print(f"\nMonte Carlo bootstrap  (seed {SEED}, {N_BOOT} resamples)")
    print(f"  RMSE 95% CI [{lo:.3f}, {hi:.3f}]   MUE 95% CI [{mlo:.3f}, {mhi:.3f}]")

    exact = list(itertools.product(errors, repeat=n))
    lo, hi = ci95([rmse(s) for s in exact])
    mlo, mhi = ci95([mue(s) for s in exact])
    print(f"\nExact bootstrap        (all {len(exact)} resamples)")
    print(f"  RMSE 95% CI [{lo:.3f}, {hi:.3f}]   MUE 95% CI [{mlo:.3f}, {mhi:.3f}]")

    uppers = []
    for seed in range(200):
        r = random.Random(seed)
        uppers.append(ci95([rmse([r.choice(errors) for _ in range(n)]) for _ in range(N_BOOT)])[1])
    print(f"\nSeed sensitivity of the RMSE upper bound over 200 seeds: "
          f"{min(uppers):.3f} to {max(uppers):.3f}")


if __name__ == "__main__":
    main()
