#!/usr/bin/env python
"""표적 도킹 대상 목록 생성 — ①앵커 9 ②후보 8 ③drugB 상위 100 ④원본 상위 100."""
import csv, json, os, sys
import numpy as np
T = "/home/nudge/Project/CADD/ttr"
sys.path.insert(0, os.path.join(T, "qsar/src"))
os.chdir(T)
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")
DPE = "c1ccc(Oc2ccccc2)cc1"
CO2 = Chem.MolFromSmarts("[CX3](=O)[OX2H1,OX1-]")
REP = ["A1", "B1", "A3"]
M9 = ["A1", "A2", "A3", "B1", "B2", "B3", "C1", "C2", "C3"]

G = [dict(smiles=r["smiles"], scaffold=r.get("scaffold", ""))
     for r in csv.DictReader(open("generation/filtered/clustered_18167.csv", encoding="utf-8"))
     if Chem.MolFromSmiles(r["smiles"]) is not None]
PB = np.load("qsar/reports/drugB_gen_predictions.npz", allow_pickle=True)
P0 = np.load("qsar/reports/D_gen_predictions.npz", allow_pickle=True)
cB = np.argsort(np.argsort(np.vstack([np.argsort(np.argsort(-PB[m])) for m in REP]).sum(0))) + 1
c9 = np.argsort(np.argsort(np.vstack([np.argsort(np.argsort(-P0[m])) for m in M9]).sum(0))) + 1

rows, seen = [], {}


def add(cid, smi, role, extra=""):
    if smi in seen:
        seen[smi]["role"] += "|" + role
        return
    r = dict(cid=cid, smiles=smi, role=role, extra=extra)
    rows.append(r); seen[smi] = r


props = {r["chembl_id"]: r for r in csv.DictReader(open("data/ttr_compound_properties.csv",
                                                        encoding="utf-8"))}
for r in csv.DictReader(open("data/anchor_assay.csv")):
    p = props.get(r["node"])
    if p:
        add(r["node"], p["smiles"], "anchor", f"pKd={r['pchembl']}")
for i in np.argsort(cB):
    if cB[i] <= 300 and G[i]["scaffold"] == DPE and Chem.MolFromSmiles(G[i]["smiles"]).HasSubstructMatch(CO2):
        add(f"cand{cB[i]:04d}", G[i]["smiles"], "cand8", f"drugB={cB[i]}")
for i in np.argsort(cB)[:100]:
    add(f"dB{cB[i]:04d}", G[i]["smiles"], "drugB100", f"drugB={cB[i]}")
for i in np.argsort(c9)[:100]:
    add(f"or{c9[i]:04d}", G[i]["smiles"], "orig100", f"orig9={c9[i]}")
from collections import Counter
print(f"★ 표적 총 {len(rows)} (중복 제거 후)")
print("  역할 분포:", dict(Counter(r["role"] for r in rows)))
with open("qsar/data/docking/targets.csv", "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=["cid", "smiles", "role", "extra"])
    w.writeheader(); w.writerows(rows)
print("저장 qsar/data/docking/targets.csv")
