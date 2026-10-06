#!/usr/bin/env python
"""작업 5-4 항목 9 — 외부 검증.

★Tox24 로 학습한 모델이 ★우리 F 축 tan040 의 역가 순위를 맞히는가.
  학습 = tox24_aromatic 전체(873) 또는 tox24_full(1,500) · 하이퍼파라미터는
  ★그 데이터 안의 골격 그룹 5-겹으로만 고른다(외부셋을 보지 않는다).
  평가 = tan040 중 ★실측 pIC50 보유 화합물.
★ 축이 다르다(결합 치환 % vs 피브릴 IC50) — prereg v4 external_validation_v4 에 명시.
"""
import csv, json, os, sys, warnings
import numpy as np
warnings.filterwarnings("ignore")
T = "/home/nudge/Project/CADD/ttr"
sys.path.insert(0, os.path.join(T, "qsar/src")); sys.path.insert(0, os.path.join(T, "ef"))
os.chdir(T)
from rdkit import Chem, DataStructs, RDLogger
RDLogger.DisableLog("rdApp.*")
from features import props20, ion_scalars, make_fp
from prereg_baseline import fp as fp2048
import model_loso as M
from model_tox24 import inner5
MIDS = ["A1", "A2", "A3", "B1", "B2", "B3", "C1", "C2", "C3"]
SEED = 42


def feat(rows, smi_key="smiles", ph_key="smiles_ph74"):
    for r in rows:
        m = Chem.MolFromSmiles(r[smi_key])
        r["_p20"] = props20(m); r["_ion"] = ion_scalars(m)
        r["_fp_neu"] = make_fp(r[smi_key]); r["_fp_ion"] = make_fp(r[ph_key])
        r["_fp2048"] = fp2048(r[smi_key])
    return rows


def boot_rho(x, y, B=2000, seed=SEED):
    rng = np.random.RandomState(seed); n = len(x); out = []
    for _ in range(B):
        i = rng.randint(0, n, n)
        if len(set(x[i])) > 1 and len(set(y[i])) > 1:
            r = M.spearman(x[i], y[i])
            if not np.isnan(r):
                out.append(r)
    out = np.sort(np.array(out))
    return [round(float(out[int(.025*len(out))]), 4), round(float(out[int(.975*len(out))]), 4)]


def main():
    # ── 외부셋: tan040 중 pIC50 보유 ──
    A = {a["chembl_id"]: a for a in csv.DictReader(open("qsar/data/interim/actives_tan040.csv"))}
    tgt = {t["chembl_id"]: t for t in csv.DictReader(open("qsar/data/interim/target_tan040.csv"))}
    ext = []
    for cid, a in A.items():
        t = tgt.get(cid, {})
        if t.get("target_pIC50"):
            r = dict(a); r["_pic50"] = float(t["target_pIC50"]); ext.append(r)
    feat(ext)
    y = np.array([r["_pic50"] for r in ext])
    print(f"★ 외부셋 = tan040 중 실측 pIC50 보유 {len(ext)} · pIC50 {y.min():.2f}~{y.max():.2f}")

    # ── 겹침 확인 (누수) ──
    res = {"n_external": len(ext), "pIC50_range": [float(y.min()), float(y.max())]}
    out = {}
    for SET in ("tox24_aromatic", "tox24_full"):
        TR = list(csv.DictReader(open(f"qsar/data/interim/{SET}.csv")))
        for r in TR:
            r["_y"] = float(r["value"])
        feat(TR)
        ik_tr = {r["inchikey"] for r in TR}
        ik_ext = {Chem.MolToInchiKey(Chem.MolFromSmiles(r["smiles"])) for r in ext}
        leak = ik_tr & ik_ext
        print(f"\n=== 학습 {SET} ({len(TR)}) · ★외부셋과 InChIKey 겹침 {len(leak)}")
        keep = [r for r in ext if Chem.MolToInchiKey(Chem.MolFromSmiles(r["smiles"])) not in ik_tr]
        yy = np.array([r["_pic50"] for r in keep])
        print(f"  겹침 제거 후 외부셋 {len(keep)}")
        # 2D 기준선: Tox24 actives 에 대한 maxTan
        act = [r for r in TR if r["label"] == "1"]
        bs = np.array([max(DataStructs.BulkTanimotoSimilarity(r["_fp2048"],
                                                              [a["_fp2048"] for a in act]))
                       for r in keep])
        brho = M.spearman(yy, bs)
        print(f"  ★2D 기준선(Tox24 actives maxTan) rho {brho:+.4f} {boot_rho(yy, bs)}")
        d = {"n_train": len(TR), "n_leak": len(leak), "n_external_used": len(keep),
             "baseline_rho": round(brho, 4), "baseline_ci": boot_rho(yy, bs), "models": {}}
        for mid in MIDS:
            hp, ins = inner5(mid, TR)
            p = np.array(M.fit_predict(mid, hp, TR, keep), float)
            rho = M.spearman(yy, p)
            ci = boot_rho(yy, p)
            d["models"][mid] = dict(label=M.LABEL[mid], hp=hp, inner_spearman=round(float(ins), 4),
                                    rho=round(float(rho), 4), ci95=ci,
                                    pred_median=round(float(np.median(p)), 3),
                                    pred_sd=round(float(p.std()), 3))
            flag = "★" if ci[0] > 0 else " "
            print(f"  {flag}{mid} {M.LABEL[mid][:36]:38s} rho {rho:+.4f} {ci} "
                  f"· 예측 중앙 {np.median(p):7.2f} sd {p.std():6.2f}")
        out[SET] = d
    res["sets"] = out
    json.dump(res, open("qsar/reports/D_external.json", "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    print("\n저장 qsar/reports/D_external.json")


if __name__ == "__main__":
    main()
