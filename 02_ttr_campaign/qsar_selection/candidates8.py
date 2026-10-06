#!/usr/bin/env python
"""작업 1 — 디페닐에터 후보 전수 기술. ★저장된 예측만 사용(새 학습 없음)."""
import csv, json, os, sys
import numpy as np
T = "/home/nudge/Project/CADD/ttr"
sys.path.insert(0, os.path.join(T, "qsar/src"))
os.chdir(T)
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFMCS
RDLogger.DisableLog("rdApp.*")
from prereg_baseline import fp as fp2048
DPE = "c1ccc(Oc2ccccc2)cc1"
CO2 = Chem.MolFromSmarts("[CX3](=O)[OX2H1,OX1-]")
OHc = Chem.MolFromSmarts("[OX2H][c]")
REP = ["A1", "B1", "A3"]
M9 = ["A1", "A2", "A3", "B1", "B2", "B3", "C1", "C2", "C3"]


def mcs_n(a, b):
    r = rdFMCS.FindMCS([a, b], timeout=5, ringMatchesRingOnly=True, completeRingsOnly=True,
                       atomCompare=rdFMCS.AtomCompare.CompareElements,
                       bondCompare=rdFMCS.BondCompare.CompareOrderExact)
    return 0 if (r.canceled or r.numAtoms == 0) else r.numAtoms


def main():
    G = [dict(smiles=r["smiles"], scaffold=r.get("scaffold", ""))
         for r in csv.DictReader(open("generation/filtered/clustered_18167.csv", encoding="utf-8"))
         if Chem.MolFromSmiles(r["smiles"]) is not None]
    for g in G:
        g["_m"] = Chem.MolFromSmiles(g["smiles"]); g["_fp"] = fp2048(g["smiles"])
    PB = np.load("qsar/reports/drugB_gen_predictions.npz", allow_pickle=True)
    P0 = np.load("qsar/reports/D_gen_predictions.npz", allow_pickle=True)
    RB = np.vstack([np.argsort(np.argsort(-PB[m])) for m in REP])
    R03 = np.vstack([np.argsort(np.argsort(-P0[m])) for m in REP])
    R09 = np.vstack([np.argsort(np.argsort(-P0[m])) for m in M9])
    cB = np.argsort(np.argsort(RB.sum(axis=0))) + 1        # drugB 합의 순위(1=최상)
    c03 = np.argsort(np.argsort(R03.sum(axis=0))) + 1      # 원본 3종
    c09 = np.argsort(np.argsort(R09.sum(axis=0))) + 1      # 원본 9종
    iqrB = np.percentile(RB, 75, axis=0) - np.percentile(RB, 25, axis=0)
    # 앵커
    props = {r["chembl_id"]: r for r in csv.DictReader(open("data/ttr_compound_properties.csv",
                                                            encoding="utf-8"))}
    A = []
    for r in csv.DictReader(open("data/anchor_assay.csv")):
        p = props.get(r["node"])
        if p:
            A.append(dict(id=r["node"], pKd=float(r["pchembl"]), smiles=p["smiles"],
                          m=Chem.MolFromSmiles(p["smiles"]), fp=fp2048(p["smiles"])))
    print(f"앵커 {len(A)} · 생성분자 {len(G)}")

    def describe(idx):
        rows = []
        for i in idx:
            g = G[i]; m = g["_m"]
            tan = [DataStructs.TanimotoSimilarity(g["_fp"], a["fp"]) for a in A]
            mcs = [mcs_n(m, a["m"]) for a in A]
            j = int(np.argmax(tan))
            rows.append(dict(
                drugB_rank=int(cB[i]), orig3_rank=int(c03[i]), orig9_rank=int(c09[i]),
                rank_IQR=int(iqrB[i]),
                carboxyl=int(m.HasSubstructMatch(CO2)),
                n_phenol=len(m.GetSubstructMatches(OHc)),
                n_halogen=sum(1 for a_ in m.GetAtoms() if a_.GetSymbol() in ("Cl", "Br", "F", "I")),
                heavy=m.GetNumHeavyAtoms(),
                best_anchor=A[j]["id"], best_anchor_pKd=A[j]["pKd"],
                best_tan=round(float(tan[j]), 4),
                mcs_best=int(mcs[j]),
                changed_vs_best=int(max(m.GetNumHeavyAtoms(), A[j]["m"].GetNumHeavyAtoms()) - mcs[j]),
                mcs_max=int(max(mcs)), tan_median=round(float(np.median(tan)), 4),
                smiles=g["smiles"]))
        return rows

    out = {}
    for N, tag in ((300, "top300_dpe_carboxyl"), (100, "top100_dpe")):
        idx = [i for i in np.argsort(cB) if cB[i] <= N and G[i]["scaffold"] == DPE]
        if tag.endswith("carboxyl"):
            idx = [i for i in idx if G[i]["_m"].HasSubstructMatch(CO2)]
        rows = describe(idx)
        out[tag] = rows
        print(f"\n{'='*100}\n★★ {tag} — {len(rows)}개")
        print(f"{'drugB':>6s} {'원본3':>6s} {'원본9':>6s} {'IQR':>5s} {'CO2H':>5s} {'OH':>3s} "
              f"{'X':>3s} {'중원자':>5s} {'최근접앵커':>14s} {'pKd':>5s} {'Tan':>6s} "
              f"{'MCS':>4s} {'변환':>5s}")
        for r in rows:
            print(f"{r['drugB_rank']:6d} {r['orig3_rank']:6d} {r['orig9_rank']:6d} "
                  f"{r['rank_IQR']:5d} {r['carboxyl']:5d} {r['n_phenol']:3d} {r['n_halogen']:3d} "
                  f"{r['heavy']:5d} {r['best_anchor']:>14s} {r['best_anchor_pKd']:5.2f} "
                  f"{r['best_tan']:6.3f} {r['mcs_best']:4d} {r['changed_vs_best']:5d}")
        print("  SMILES:")
        for r in rows:
            print(f"    #{r['drugB_rank']:<4d} {r['smiles']}")
    # 학습셋 의존성
    print(f"\n{'='*100}\n★★ 학습셋 의존성 — 두 학습셋 모두에서 상위권인가")
    for tag in out:
        rows = out[tag]
        for th in (300, 1000):
            ok = [r for r in rows if r["orig9_rank"] <= th or r["orig3_rank"] <= th]
            print(f"  {tag}: 원본(3종 또는 9종) 상위 {th} 안에도 드는 것 "
                  f"★{len(ok)}/{len(rows)} — drugB 순위 {[r['drugB_rank'] for r in ok]}")
    json.dump(out, open("qsar/reports/candidates8.json", "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    with open("qsar/reports/candidates_dpe.csv", "w", newline="", encoding="utf-8") as fh:
        allr = [dict(set=k, **r) for k, v in out.items() for r in v]
        wtr = csv.DictWriter(fh, fieldnames=list(allr[0].keys())); wtr.writeheader(); wtr.writerows(allr)
    print("\n저장 qsar/reports/candidates8.json · candidates_dpe.csv")


if __name__ == "__main__":
    main()
