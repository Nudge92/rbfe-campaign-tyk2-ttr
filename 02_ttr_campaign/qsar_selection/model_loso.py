#!/usr/bin/env python
"""4단계 — 기술자 3안 × 모델 9종, LOSO 평가. prereg v2 지표.

★ 모델 수 9 (prereg 상한 12) — 사전 선언:
   A1 A-Ridge · A2 A-GP(RBF) · A3 A-RF(depth<=4)
   B1 B-KRR(RBF(props) ⊕ Tanimoto(이온화형 FP)) · B2 B-kNN(Tanimoto) · B3 B-RF(props+유사도요약)
   C1 C-KRR(Tanimoto(중성 FP) ⊕ RBF(이온화 스칼라)) · C2 C-kNN(Tanimoto 중성) · C3 C-RF(props+이온+유사도요약)
★ 중첩 CV — 하이퍼파라미터는 train 안의 LOSO 로만 고른다. 스케일링·유사도요약도 폴드 안에서만.
★ 딥러닝·GNN 없음.
"""
import csv, json, os, sys, warnings
import numpy as np
warnings.filterwarnings("ignore")
T = "/home/nudge/Project/CADD/ttr"
sys.path.insert(0, os.path.join(T, "qsar/src"))
os.chdir(T)
from sklearn.linear_model import Ridge
from sklearn.kernel_ridge import KernelRidge
from sklearn.ensemble import RandomForestRegressor
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel, WhiteKernel
from sklearn.preprocessing import StandardScaler
from scipy.stats import rankdata
from rdkit import Chem

from features import (props20, ion_scalars, make_fp, tanimoto_matrix, sim_summary,
                      PROP_NAMES, ION_NAMES)
from metrics import roc_auc, ef, ef_band, bedroc, bedroc_band, tie_diagnostics, rank_percentile

SET = os.environ.get("QSAR_SET", "S0_covfree")
SEED = 42
N_BOOT = 2000


# ───────────────────────── 데이터 ─────────────────────────
def load():
    A = list(csv.DictReader(open(f"qsar/data/interim/actives_{SET}.csv")))
    D = list(csv.DictReader(open(f"qsar/data/interim/decoys_{SET}.csv")))
    print(f"[{SET}] actives {len(A)} · decoys {len(D)}")
    # ★F 축 회귀 타깃 (prereg v3 regression_target_v3) — pIC50 = 9 − log10(IC50_nM)
    tgt = {}
    tp = f"qsar/data/interim/target_{SET}.csv"
    if os.path.exists(tp):
        for r in csv.DictReader(open(tp)):
            if r["target_pIC50"]:
                tgt[r["chembl_id"]] = float(r["target_pIC50"])
    for a in A:
        m = Chem.MolFromSmiles(a["smiles"])
        mi = Chem.MolFromSmiles(a["smiles_ph74"])
        a["_p20"] = props20(m); a["_ion"] = ion_scalars(m)
        a["_fp_neu"] = make_fp(a["smiles"]); a["_fp_ion"] = make_fp(a["smiles_ph74"])
        # ★F 축은 pchembl 이 없다 → target_{SET}.csv 의 pIC50 을 쓴다.
        # ★타깃이 없는 actives 는 학습에서만 빠지고 순위 평가에는 남는다(%FF→IC50 변환 금지).
        if a.get("pchembl"):
            a["_y"] = float(a["pchembl"]); a["_has_y"] = True
        elif a["chembl_id"] in tgt:
            a["_y"] = tgt[a["chembl_id"]]; a["_has_y"] = True
        else:
            a["_y"] = float("nan"); a["_has_y"] = False
    for d in D:
        m = Chem.MolFromSmiles(d["smiles"])
        d["_p20"] = props20(m); d["_ion"] = ion_scalars(m)
        d["_fp_neu"] = make_fp(d["smiles"])
        pr, _ = __import__("protonate").protonate(m)
        d["_fp_ion"] = make_fp(Chem.MolToSmiles(pr))
    ny = sum(1 for a in A if a["_has_y"])
    print(f"  ★회귀 타깃 보유 {ny}/{len(A)} · 타깃 없는 {len(A) - ny} 는 학습 제외(순위 평가는 포함)")
    return A, D


# ───────────────────── 모델 (fit/predict) ─────────────────────
def _gp_rbf(alpha):
    k = ConstantKernel(1.0, (1e-3, 1e3)) * RBF(1.0, (1e-2, 1e3)) + WhiteKernel(alpha, (1e-5, 1e1))
    return GaussianProcessRegressor(kernel=k, normalize_y=True, random_state=SEED, n_restarts_optimizer=2)


GRIDS = {
    "A1": [{"alpha": a} for a in (0.1, 1.0, 10.0, 100.0)],
    "A2": [{"alpha": a} for a in (1e-2, 1e-1, 1.0)],
    "A3": [{"max_depth": d, "max_features": f} for d in (2, 3, 4) for f in (0.3, 0.6)],
    "B1": [{"alpha": a, "w": w} for a in (1e-2, 1e-1, 1.0) for w in (0.3, 0.5, 0.7)],
    "B2": [{"k": k} for k in (1, 3, 5, 7)],
    "B3": [{"max_depth": d, "max_features": f} for d in (2, 3, 4) for f in (0.3, 0.6)],
    "C1": [{"alpha": a, "w": w} for a in (1e-2, 1e-1, 1.0) for w in (0.3, 0.5, 0.7)],
    "C2": [{"k": k} for k in (1, 3, 5, 7)],
    "C3": [{"max_depth": d, "max_features": f} for d in (2, 3, 4) for f in (0.3, 0.6)],
}
DESC = {"A1": "A", "A2": "A", "A3": "A", "B1": "B", "B2": "B", "B3": "B",
        "C1": "C", "C2": "C", "C3": "C"}
LABEL = {
    "A1": "A-Ridge (물성20)", "A2": "A-GP(RBF) (물성20)", "A3": "A-RF d<=4 (물성20)",
    "B1": "B-KRR (RBF물성 ⊕ Tanimoto 이온화FP)", "B2": "B-kNN (Tanimoto 이온화FP)",
    "B3": "B-RF d<=4 (물성 + 유사도요약3)",
    "C1": "C-KRR (Tanimoto 중성FP ⊕ RBF 이온스칼라)", "C2": "C-kNN (Tanimoto 중성FP)",
    "C3": "C-RF d<=4 (물성 + 이온스칼라 + 유사도요약3)",
}


def fit_predict(mid, hp, tr_rows, te_rows):
    """★train 통계(스케일러·유사도 기준)는 tr_rows 에서만 만든다."""
    ytr = np.array([r["_y"] for r in tr_rows])
    if mid in ("A1", "A2", "A3"):
        Xtr = np.array([r["_p20"] for r in tr_rows]); Xte = np.array([r["_p20"] for r in te_rows])
        sc = StandardScaler().fit(Xtr); Xtr, Xte = sc.transform(Xtr), sc.transform(Xte)
        if mid == "A1":
            m = Ridge(alpha=hp["alpha"]).fit(Xtr, ytr)
        elif mid == "A2":
            m = _gp_rbf(hp["alpha"]).fit(Xtr, ytr)
        else:
            m = RandomForestRegressor(n_estimators=500, max_depth=hp["max_depth"],
                                      max_features=hp["max_features"], random_state=SEED,
                                      n_jobs=1).fit(Xtr, ytr)
        return m.predict(Xte)

    fpkey = "_fp_ion" if DESC[mid] == "B" else "_fp_neu"
    ftr = [r[fpkey] for r in tr_rows]; fte = [r[fpkey] for r in te_rows]

    if mid in ("B2", "C2"):
        Tte = tanimoto_matrix(fte, ftr)
        k = min(hp["k"], len(ftr))
        idx = np.argsort(-Tte, axis=1)[:, :k]
        w = np.take_along_axis(Tte, idx, axis=1) + 1e-9
        return (ytr[idx] * w).sum(axis=1) / w.sum(axis=1)

    if mid in ("B1", "C1"):
        Ktr_fp = tanimoto_matrix(ftr, ftr); Kte_fp = tanimoto_matrix(fte, ftr)
        side = "_p20" if mid == "B1" else "_ion"
        Str = np.array([r[side] for r in tr_rows]); Ste = np.array([r[side] for r in te_rows])
        sc = StandardScaler().fit(Str); Str, Ste = sc.transform(Str), sc.transform(Ste)
        g = 1.0 / max(Str.shape[1], 1)
        d_tr = ((Str[:, None] - Str[None]) ** 2).sum(-1); d_te = ((Ste[:, None] - Str[None]) ** 2).sum(-1)
        Ktr_s, Kte_s = np.exp(-g * d_tr), np.exp(-g * d_te)
        w = hp["w"]
        Ktr = w * Ktr_fp + (1 - w) * Ktr_s
        Kte = w * Kte_fp + (1 - w) * Kte_s
        m = KernelRidge(alpha=hp["alpha"], kernel="precomputed").fit(Ktr, ytr)
        return m.predict(Kte)

    # B3 · C3 — 트리 + 유사도 요약 (★train 에 대한 유사도만)
    Ttr = tanimoto_matrix(ftr, ftr).copy()
    np.fill_diagonal(Ttr, -1.0)                      # ★자기제외
    Str3 = sim_summary(Ttr)
    Ste3 = sim_summary(tanimoto_matrix(fte, ftr))
    if mid == "B3":
        Xtr = np.hstack([np.array([r["_p20"] for r in tr_rows]), Str3])
        Xte = np.hstack([np.array([r["_p20"] for r in te_rows]), Ste3])
    else:
        Xtr = np.hstack([np.array([r["_p20"] for r in tr_rows]),
                         np.array([r["_ion"] for r in tr_rows]), Str3])
        Xte = np.hstack([np.array([r["_p20"] for r in te_rows]),
                         np.array([r["_ion"] for r in te_rows]), Ste3])
    sc = StandardScaler().fit(Xtr)
    m = RandomForestRegressor(n_estimators=500, max_depth=hp["max_depth"],
                              max_features=hp["max_features"], random_state=SEED,
                              n_jobs=1).fit(sc.transform(Xtr), ytr)
    return m.predict(sc.transform(Xte))


def spearman(x, y):
    if len(x) < 3 or len(set(x)) < 2 or len(set(y)) < 2:
        return float("nan")
    rx, ry = rankdata(x), rankdata(y)
    mx, my = rx.mean(), ry.mean()
    num = ((rx - mx) * (ry - my)).sum()
    den = np.sqrt(((rx - mx) ** 2).sum() * ((ry - my) ** 2).sum())
    return float(num / den) if den else float("nan")


def inner_select(mid, tr_rows):
    """★중첩 CV — train 안의 LOSO 로만 하이퍼파라미터를 고른다."""
    folds = sorted({r["murcko"] for r in tr_rows})
    best, best_s = None, -np.inf
    for hp in GRIDS[mid]:
        preds, ys = [], []
        for f in folds:
            itr = [r for r in tr_rows if r["murcko"] != f]
            ite = [r for r in tr_rows if r["murcko"] == f]
            if len(itr) < 5 or not ite:
                continue
            try:
                p = fit_predict(mid, hp, itr, ite)
            except Exception:
                p = np.full(len(ite), np.nan)
            preds.extend(p); ys.extend([r["_y"] for r in ite])
        preds, ys = np.array(preds, float), np.array(ys, float)
        ok = ~np.isnan(preds)
        if ok.sum() < 5:
            continue
        s = spearman(ys[ok], preds[ok])
        if np.isnan(s):
            s = -np.mean((ys[ok] - preds[ok]) ** 2)      # 퇴화 시 −MSE 로 대체
        if s > best_s:
            best_s, best = s, hp
    return best if best is not None else GRIDS[mid][0], best_s


# ───────────────────── 평가 ─────────────────────
def cluster_bootstrap(folds, y, s_model, s_base, B=N_BOOT, seed=SEED):
    """★골격 클러스터 단위 부트스트랩 → Δ(모델 − 기준선) 95% CI."""
    rng = np.random.RandomState(seed)
    uf = sorted(set(folds))
    out = {"EF5": [], "BEDROC": [], "ROC_AUC": []}
    for _ in range(B):
        pick = rng.choice(len(uf), len(uf), replace=True)
        idx = np.concatenate([np.where(folds == uf[i])[0] for i in pick])
        yy = y[idx]
        if yy.sum() == 0 or yy.sum() == len(yy):
            continue
        sm, sb = s_model[idx], s_base[idx]
        out["EF5"].append(ef(sm, yy, 0.05) - ef(sb, yy, 0.05))
        out["BEDROC"].append(bedroc(sm, yy) - bedroc(sb, yy))
        out["ROC_AUC"].append(roc_auc(sm, yy) - roc_auc(sb, yy))
    res = {}
    for k, v in out.items():
        v = np.sort(np.array(v, float))
        v = v[~np.isnan(v)]
        res[k] = dict(n=len(v), lo=round(float(v[int(.025 * len(v))]), 4),
                      hi=round(float(v[int(.975 * len(v))]), 4),
                      med=round(float(np.median(v)), 4)) if len(v) > 50 else None
    return res


def boot_metric_ci(folds, vals, B=N_BOOT, seed=SEED, fn=np.mean):
    rng = np.random.RandomState(seed)
    uf = sorted(set(folds))
    out = []
    for _ in range(B):
        pick = rng.choice(len(uf), len(uf), replace=True)
        idx = np.concatenate([np.where(folds == uf[i])[0] for i in pick])
        out.append(fn(vals[idx]))
    out = np.sort(np.array(out, float))
    return round(float(out[int(.025 * len(out))]), 4), round(float(out[int(.975 * len(out))]), 4)


def boot_spearman_ci(folds, x, y, B=N_BOOT, seed=SEED):
    rng = np.random.RandomState(seed)
    uf = sorted(set(folds)); out = []
    for _ in range(B):
        pick = rng.choice(len(uf), len(uf), replace=True)
        idx = np.concatenate([np.where(folds == uf[i])[0] for i in pick])
        r = spearman(x[idx], y[idx])
        if not np.isnan(r):
            out.append(r)
    if len(out) < 50:
        return None, None
    out = np.sort(np.array(out))
    return round(float(out[int(.025 * len(out))]), 4), round(float(out[int(.975 * len(out))]), 4)


def bh_fdr(pvals):
    p = np.asarray(pvals, float); n = len(p)
    o = np.argsort(p); q = np.empty(n)
    prev = 1.0
    for rank, i in enumerate(o[::-1]):
        k = n - rank
        prev = min(prev, p[i] * n / k)
        q[i] = prev
    return q


def boot_pvalue(deltas):
    """부트스트랩 분포에서 Δ<=0 인 비율 (단측)."""
    d = np.asarray(deltas, float)
    d = d[~np.isnan(d)]
    if len(d) == 0:
        return 1.0
    return float(max(1.0 / len(d), (d <= 0).mean()))


# ───────────────────── 폴드 단위 병렬 실행 ─────────────────────
# ★폴드끼리 완전히 독립이고 random_state 가 고정이므로 워커 수와 무관하게 결과가 같다.
#   fork 로 전역(_G)을 물려주므로 큰 배열을 태스크마다 피클하지 않는다.
#   ★RF 는 n_jobs=1 — 46행 학습에서 스레드 4개는 이득이 없고(498 vs 534 ms) 예측은 완전 동일하다.
_G = {}
NW = int(os.environ.get("QSAR_NW", "1"))


def _fold_job(args):
    """한 (모델, 폴드) 조합. 순수 함수 — 전역은 읽기만 한다."""
    mid, f, singleton = args
    A, clus, byact = _G["A"], _G["clus"], _G["byact"]
    te_act = clus[f]
    tr_act = [a for a in A if a["murcko"] != f and a["_has_y"]]
    hp, inner_s = inner_select(mid, tr_act)
    te_dec = [d for a in te_act for d in byact.get(a["chembl_id"], [])]
    te_rows = te_act + te_dec
    p = fit_predict(mid, hp, tr_act, te_rows)
    ids = [r.get("chembl_id", r.get("zinc_id")) for r in te_rows]
    return (mid, f, singleton, hp, float(inner_s), ids, [float(v) for v in p],
            len(te_act), len(te_dec), len(tr_act))


def run_jobs(tasks, tag=""):
    if NW <= 1:
        return [_fold_job(t) for t in tasks]
    import multiprocessing as mp
    ctx = mp.get_context("fork")
    done = []
    with ctx.Pool(NW) as pool:
        for i, r in enumerate(pool.imap_unordered(_fold_job, tasks, chunksize=1), 1):
            done.append(r)
            if i % 25 == 0 or i == len(tasks):
                print(f"    [{tag}] 진행 {i}/{len(tasks)}", flush=True)
    return done


def main():
    A, D = load()
    byact = {}
    for d in D:
        byact.setdefault(d["active_id"], []).append(d)
    clus = {}
    for a in A:
        clus.setdefault(a["murcko"], []).append(a)
    eval_folds = [k for k, v in clus.items() if len(v) >= 2]
    sing_folds = [k for k, v in clus.items() if len(v) == 1]
    print(f"골격 {len(clus)} · 평가 폴드(크기>=2) {len(eval_folds)} · 단일 폴드 {len(sing_folds)}")

    # 2D 기준선 점수 (같은 규칙: LSO max Tanimoto · ECFP4 2048 — prereg 와 동일)
    from prereg_baseline import score_lso
    ev_ids = {a["chembl_id"] for k in eval_folds for a in clus[k]}
    bid, by, bs, bfold = score_lso(A, D, ev_ids)
    base_score = {i: s for i, s in zip(bid, bs)}

    # ★전역을 채운 뒤 fork — 워커는 읽기만 한다
    _G["A"], _G["clus"], _G["byact"] = A, clus, byact
    MIDS = ["A1", "A2", "A3", "B1", "B2", "B3", "C1", "C2", "C3"]
    print(f"★ 병렬 워커 {NW} · 평가 폴드 작업 {len(MIDS) * len(eval_folds)}개", flush=True)
    RES = {}
    for r in run_jobs([(m, f, False) for m in MIDS for f in eval_folds], "평가폴드"):
        RES[(r[0], r[1])] = r

    rows_matrix = []
    summary = {}
    for mid in MIDS:
        oof = {}
        per_fold = []
        hp_used = {}
        for f in eval_folds:
            te_act = clus[f]
            _, _, _, hp, inner_s, ids_f, p_f, na, n_dec, n_tr = RES[(mid, f)]
            hp_used[f] = hp
            for i_, v in zip(ids_f, p_f):
                oof[i_] = v
            p = np.array(p_f, float)
            te_dec = [None] * n_dec
            tr_act = [None] * n_tr
            # ★역가 지표는 타깃 보유 홀드아웃 actives 에서만 — n 을 반드시 병기한다
            hm = np.array([r["_has_y"] for r in te_act], bool)
            ya = np.array([r["_y"] for r in te_act], float)[hm]
            pa = np.array(p[:na], float)[hm]
            yy = np.array([1] * na + [0] * len(te_dec))
            ss = np.array(p, float)
            per_fold.append(dict(
                model=mid, fold=f, n_actives=na, n_potency=int(hm.sum()),
                n_decoys=len(te_dec), n_train=len(tr_act),
                hp=json.dumps(hp, ensure_ascii=False), inner_spearman=round(float(inner_s), 4),
                ROC_AUC=round(roc_auc(ss, yy), 4), EF5=round(ef(ss, yy, 0.05), 3),
                BEDROC=round(bedroc(ss, yy), 4),
                spearman_actives=None if len(ya) < 3 else round(spearman(ya, pa), 4),
                RMSE=None if not len(ya) else round(float(np.sqrt(((ya - pa) ** 2).mean())), 4),
                MAE=None if not len(ya) else round(float(np.abs(ya - pa).mean()), 4),
                pred_std=None if not len(ya) else round(float(pa.std()), 4),
                exp_std=None if not len(ya) else round(float(ya.std()), 4),
            ))
        # pooled
        ids = list(bid)
        sm = np.array([oof[i] for i in ids], float)
        y = by; folds = bfold
        act_mask = y == 1
        id2y = {a["chembl_id"]: (a["_y"] if a["_has_y"] else np.nan)
                for k in eval_folds for a in clus[k]}
        exp_a = np.array([id2y[i] for i in np.array(ids)[act_mask]], float)
        pred_a = sm[act_mask]
        fold_a = folds[act_mask]
        pm = ~np.isnan(exp_a)                      # ★역가 평가 가능한 actives
        exp_a, pred_a, fold_a = exp_a[pm], pred_a[pm], fold_a[pm]
        n_pot = int(pm.sum())
        rho = spearman(exp_a, pred_a)
        rlo, rhi = boot_spearman_ci(fold_a, exp_a, pred_a)
        err = np.abs(exp_a - pred_a)
        rmse = float(np.sqrt((err ** 2).mean())); mae = float(err.mean())
        rlo_r, rhi_r = boot_metric_ci(fold_a, (exp_a - pred_a) ** 2, fn=lambda v: np.sqrt(v.mean()))
        mlo, mhi = boot_metric_ci(fold_a, err, fn=np.mean)
        vr = float(pred_a.std() / exp_a.std())
        bs_res = cluster_bootstrap(folds, y, sm, bs)
        summary[mid] = dict(
            label=LABEL[mid], descriptor=DESC[mid],
            ROC_AUC=round(roc_auc(sm, y), 4), EF5=round(ef(sm, y, 0.05), 3),
            EF1=round(ef(sm, y, 0.01), 3), BEDROC=round(bedroc(sm, y), 4),
            delta_vs_baseline=bs_res,
            held_out_actives_spearman=dict(pooled=round(rho, 4), n=n_pot,
                                           ci95=[rlo, rhi] if rlo is not None else None),
            RMSE=round(rmse, 4), RMSE_ci=[rlo_r, rhi_r], MAE=round(mae, 4), MAE_ci=[mlo, mhi],
            variance_ratio=round(vr, 4),
            tie_bands=dict(EF5=ef_band(sm, y, 0.05), BEDROC=bedroc_band(sm, y)),
            tie_diagnostics=tie_diagnostics(sm, y),
            hp_per_fold={k: v for k, v in hp_used.items()},
        )
        rows_matrix.extend(per_fold)
        print(f"  {mid} {LABEL[mid][:40]:42s} AUC {summary[mid]['ROC_AUC']:.4f} "
              f"EF5 {summary[mid]['EF5']:.2f} BED {summary[mid]['BEDROC']:.4f} "
              f"rho {summary[mid]['held_out_actives_spearman']['pooled']:+.3f} vr {vr:.2f}", flush=True)

    # 단일 폴드 (★분리 보고)
    print(f"★ 단일 폴드 작업 {len(MIDS) * len(sing_folds)}개", flush=True)
    SRES = {}
    for r in run_jobs([(m, f, True) for m in MIDS for f in sing_folds], "단일폴드"):
        SRES[(r[0], r[1])] = r
    sing = {}
    for mid in MIDS:
        rec = []
        for f in sing_folds:
            a = clus[f][0]
            _, _, _, hp, _is, _ids, p_s, _na, n_dec, _nt = SRES[(mid, f)]
            p = np.array(p_s, float)
            te_dec = [None] * n_dec
            rec.append(dict(chembl_id=a["chembl_id"],
                            pchembl=(a["_y"] if a["_has_y"] else None),
                            murcko=a["murcko"], n_decoys=len(te_dec),
                            rank_percentile=round(rank_percentile(float(p[0]), np.array(p[1:], float)), 4)))
        sing[mid] = rec
        rp = [r["rank_percentile"] for r in rec]
        taf = [r for r in rec if r["chembl_id"] == "CHEMBL2103837"]
        print(f"  [단일폴드] {mid} 중앙 rp {np.median(rp):.3f} · 타파미디스 "
              f"{taf[0]['rank_percentile'] if taf else '-'}")

    with open(f"qsar/reports/A3_loso_matrix_{SET}.csv" if SET != "S0_covfree"
              else "qsar/reports/04_loso_matrix.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows_matrix[0].keys()))
        w.writeheader(); w.writerows(rows_matrix)
    json.dump({"summary": summary, "singleton_folds": sing,
               "n_models_declared": 9, "seed": SEED, "set": SET},
              open(f"qsar/reports/A3_results_{SET}.json" if SET != "S0_covfree"
                   else "qsar/reports/04_results.json", "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    tag = "04" if SET == "S0_covfree" else f"A3_{SET}"
    print(f"\n저장 qsar/reports/ {tag} matrix·results")


if __name__ == "__main__":
    main()
