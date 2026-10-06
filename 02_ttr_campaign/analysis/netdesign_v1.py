#!/usr/bin/env python3
"""그물 재설계(설계만) — ★시뮬레이션 0. 기존에 돌린 엣지를 살리고, 구조 기준을 채우는 데 필요한 ★최소 추가 엣지를 찾는다.

입력: qsar/reports/cooh45_v9.csv (10 노드 45쌍의 변환 난이도 — cooh_fail·dummies_heavy) · edges_v7.txt · edges_v12_add.txt
      final_report_20260921/roundtrip.json (수렴 플래그가 붙은 엣지)
기준 (prereg v37 에 계산 전 등재):
  K1 앵커 삼각형(N01·N04·N07) 포함
  K2 이중연결 — 어느 노드 하나를 빼도 나머지가 연결돼 있다
  K3 모든 노드가 길이 ≤4 인 단순 고리 2개 이상에 속한다 (network_design.md 의 기존 기준)
  K4 모든 엣지가 길이 ≤5 인 단순 고리 2개 이상에 속한다
  K5 앵커가 없는 치환 패턴 부류(ortho·meta)마다 앵커 부류(para)로 가는 직접 엣지가 2개 이상이고, 그중 두 개는 양 끝 노드를 공유하지 않는다
  K6 최대 차수 ≤5
목적 (사전순): ① 추가하는 어려운 엣지(cooh_fail=1) 수 최소 → ② 추가 엣지 수 최소 → ③ 추가 엣지의 dummies_heavy 합 최소 → ④ 엣지 이름 사전순
두 판: [유지17] 돌린 17엣지를 전부 유지 · [유지14] 수렴 플래그 엣지를 빼고 14엣지만 유지(뺀 엣지도 다시 추가할 후보에 든다)
탐색 상한: 추가 7엣지. 그 안에 해가 없으면 '없음'으로 적고 기준을 바꾸지 않는다.

사용: netdesign_v1.py --pretest [--json <경로>]     고리 세기가 기존 기록(RULER_AND_NETWORK_20260929.md)을 재현하는지 + 음성 대조
      netdesign_v1.py --run <출력 디렉터리>         prereg 마지막 항목이 v37 등재이고 스크립트 sha256 이 같을 때만 돈다
"""
import csv, hashlib, itertools, json, re, subprocess, sys
from pathlib import Path

T = Path("/home/nudge/Project/CADD/ttr"); Q = T / "qsar"; B = Q / "data/fep_m1"
REGKEY = "network_redesign_paper_only_v37"
ANCH = ["CHEMBL240808", "CHEMBL241454", "CHEMBL438498"]
MAXDEG, MAXADD, L_NODE, L_EDGE = 5, 7, 4, 5


def sha(p):
    h = hashlib.sha256(Path(p).read_bytes()).hexdigest(); assert len(h) == 64; return h


def load():
    pairs = {}
    for r in csv.DictReader(open(Q / "reports/cooh45_v9.csv", encoding="utf-8")):
        k = frozenset((r["a"], r["b"]))
        pairs[k] = dict(a=r["a"], b=r["b"], pat={r["a"]: r["pat_a"], r["b"]: r["pat_b"]}, hard=int(r["cooh_fail"]), dum=int(r["dummies_heavy"]))
    assert len(pairs) == 45
    nodes = sorted({x for k in pairs for x in k}); assert len(nodes) == 10
    pat = {}
    for k, v in pairs.items(): pat.update(v["pat"])
    run = {}
    for fn in ("edges_v7.txt", "edges_v12_add.txt"):
        for ln in open(B / fn):
            q = ln.split()
            if len(q) == 3: run[q[0]] = frozenset((q[1], q[2]))
    assert len(run) == 17 and all(k in pairs for k in run.values())
    rt = json.load(open(Q / "docs/stage1_out/final_report_20260921/roundtrip.json", encoding="utf-8"))
    flagged = sorted(e for e in run if e in rt and any(isinstance(v, dict) and v.get("flag") for v in rt[e].values()))
    return pairs, nodes, pat, run, flagged


def cycles(nodes, edges, L):
    """길이 ≤L 인 단순 고리 전부 (노드 튜플). 시작 = 고리의 최소 노드 · 방향 중복 제거."""
    idx = {n: i for i, n in enumerate(nodes)}; adj = {i: set() for i in range(len(nodes))}
    for e in edges:
        a, b = tuple(e); adj[idx[a]].add(idx[b]); adj[idx[b]].add(idx[a])
    out = []
    def dfs(s, path):
        u = path[-1]
        for v in adj[u]:
            if v == s and len(path) >= 3 and path[1] < path[-1]: out.append(tuple(path))
            elif v > s and v not in path and len(path) < L: dfs(s, path + [v])
    for s in range(len(nodes)): dfs(s, [s])
    return [tuple(nodes[i] for i in c) for c in out]


def metrics(nodes, edges):
    cy5 = cycles(nodes, edges, L_EDGE); cy4 = [c for c in cy5 if len(c) <= L_NODE]
    npart = {n: 0 for n in nodes}; epart = {e: 0 for e in edges}
    for c in cy4:
        for n in c: npart[n] += 1
    for c in cy5:
        for i in range(len(c)):
            e = frozenset((c[i], c[(i + 1) % len(c)])); epart[e] += 1
    deg = {n: sum(1 for e in edges if n in e) for n in nodes}
    return dict(ncyc4=len(cy4), ncyc5=len(cy5), npart=npart, epart=epart, deg=deg)


def connected(nodes, edges, skip=None):
    ns = [n for n in nodes if n != skip]; seen = {ns[0]}; st = [ns[0]]
    while st:
        u = st.pop()
        for e in edges:
            if u in e and skip not in e:
                v = next(x for x in e if x != u)
                if v not in seen: seen.add(v); st.append(v)
    return len(seen) == len(ns)


def k5(edges, pat):
    res = {}
    for cls in ("ortho", "meta"):
        br = [e for e in edges if sorted(pat[x] for x in e) == sorted(("para", cls))]
        disj = any(not (a & b) for a, b in itertools.combinations(br, 2))
        res[cls] = (len(br), disj)
    return res


def check(nodes, edges, pat, tri):
    """기준 위반 목록(빈 목록 = 전부 충족)."""
    bad = []
    if not all(t in edges for t in tri): bad.append("K1")
    if not connected(nodes, edges) or any(not connected(nodes, edges, skip=n) for n in nodes): bad.append("K2")
    m = metrics(nodes, edges)
    if any(v < 2 for v in m["npart"].values()): bad.append("K3")
    if any(v < 2 for v in m["epart"].values()): bad.append("K4")
    if any(n < 2 or not d for n, d in k5(edges, pat).values()): bad.append("K5")
    if max(m["deg"].values()) > MAXDEG: bad.append("K6")
    return bad, m


def search(nodes, pairs, pat, base, tri):
    cand = sorted((k for k in pairs if k not in base), key=lambda k: tuple(sorted(k)))
    easy = [k for k in cand if not pairs[k]["hard"]]; hard = [k for k in cand if pairs[k]["hard"]]
    degb = {n: sum(1 for e in base if n in e) for n in nodes}
    best, feas = None, []
    for kh in range(0, MAXADD + 1):
        for hs in itertools.combinations(hard, kh):
            for es in itertools.chain.from_iterable(itertools.combinations(easy, j) for j in range(len(easy) + 1)):
                add = hs + es
                if len(add) > MAXADD: continue
                d = dict(degb)
                for e in add:
                    for n in e: d[n] += 1
                if max(d.values()) > MAXDEG: continue
                edges = set(base) | set(add)
                bad, m = check(nodes, edges, pat, tri)
                if bad: continue
                key = (kh, len(add), sum(pairs[e]["dum"] for e in add), tuple(sorted(tuple(sorted(e)) for e in add)))
                feas.append((key, add))
        if feas: break          # ① 어려운 엣지 수가 최소인 층에서 멈춘다
    feas.sort(key=lambda x: x[0])
    return feas, len(cand), len(easy), len(hard)


def pretest():
    pairs, nodes, pat, run, flagged = load(); res = {}; ok = True
    ruler = (Q / "docs/RULER_AND_NETWORK_20260929.md").read_text(encoding="utf-8")
    m16 = re.search(r"동결 실행 그물 \(N01~N16\)\*\*\s*\|\s*(\d+)\s*\|", ruler); m17 = re.search(r"설계 전체 \(\+N17\)\s*\|\s*(\d+)\s*\|", ruler)
    e16 = {v for k, v in run.items() if k != "N17"}; e17 = set(run.values())
    a16, a17 = metrics(nodes, e16), metrics(nodes, e17)
    res["executed16_cycles_len4"] = [a16["ncyc4"], int(m16.group(1))]; res["plus_N17_cycles_len4"] = [a17["ncyc4"], int(m17.group(1))]
    lt16 = {n: v for n, v in a16["npart"].items() if v < 2}; lt17 = {n: v for n, v in a17["npart"].items() if v < 2}
    exp16 = {"cand0028": 0, "cand0118": 1, "cand0042": 1, "cand0132": 1, "cand0271": 1}; exp17 = {"cand0028": 0, "cand0118": 1, "cand0132": 1}
    res["executed16_nodes_lt2"] = [lt16, exp16]; res["plus_N17_nodes_lt2"] = [lt17, exp17]
    ok &= a16["ncyc4"] == int(m16.group(1)) and a17["ncyc4"] == int(m17.group(1)) and lt16 == exp16 and lt17 == exp17
    # 음성 대조: 엣지 하나(N12)를 빼면 고리 수가 줄어야 한다
    neg = metrics(nodes, e16 - {run["N12"]}); res["negctl_drop_N12_cycles_len4"] = neg["ncyc4"]; ok &= neg["ncyc4"] < a16["ncyc4"]
    # 기준 판정기: 돌린 17엣지는 기준을 어겨야 한다(어떤 기준인지 기록) · 완전 그래프는 K6 만 어겨야 한다
    bad17, _ = check(nodes, e17, pat, [run["N01"], run["N04"], run["N07"]]); res["executed17_violations"] = bad17; ok &= len(bad17) > 0
    badK, _ = check(nodes, set(pairs), pat, [run["N01"], run["N04"], run["N07"]]); res["complete_graph_violations"] = badK; ok &= badK == ["K6"]
    res["flagged"] = flagged; res["script_sha256"] = sha(__file__); res["all_as_expected"] = bool(ok)
    res["at"] = subprocess.run(["date", "-Is"], capture_output=True, text=True).stdout.strip()
    if "--json" in sys.argv: Path(sys.argv[sys.argv.index("--json") + 1]).write_text(json.dumps(res, ensure_ascii=False, indent=1), encoding="utf-8")
    print(json.dumps(res, ensure_ascii=False, indent=1)); sys.exit(0 if ok else 1)


def main(out):
    pr = json.load(open(Q / "models/prereg.json", encoding="utf-8")); last = pr["changelog"][-1]
    if last.get("키") != REGKEY or last["(h)★스크립트"]["sha256"] != sha(__file__):
        print("★거부: prereg 마지막 항목이 이 설계의 등재가 아니거나 스크립트가 등재 뒤 바뀌었다"); sys.exit(2)
    pairs, nodes, pat, run, flagged = load(); name = {v: k for k, v in run.items()}
    tri = [run["N01"], run["N04"], run["N07"]]
    out.mkdir(parents=True, exist_ok=False)
    result = {}; rows = []
    for tag, keep in (("유지17", sorted(run)), ("유지14", [e for e in sorted(run) if e not in flagged])):
        base = {run[e] for e in keep}
        bad0, m0 = check(nodes, base, pat, tri)
        feas, nc, ne, nh = search(nodes, pairs, pat, base, tri)
        r = dict(keep=keep, base_violations=bad0, n_candidates=nc, n_easy=ne, n_hard=nh, n_optimal_layer=len(feas))
        if feas:
            key, add = feas[0]; edges = base | set(add); _, m = check(nodes, edges, pat, tri)
            ties = [f for f in feas if f[0][:3] == key[:3]]
            r.update(add=[dict(a=sorted(e)[0], b=sorted(e)[1], hard=pairs[e]["hard"], dum=pairs[e]["dum"], was=name.get(e, "")) for e in add],
                     n_add=len(add), n_add_hard=key[0], dum_sum=key[2], n_ties=len(ties),
                     ties=[[sorted(e) for e in f[1]] for f in ties[:10]],
                     n_edges=len(edges), ncyc4=m["ncyc4"], ncyc5=m["ncyc5"], npart=m["npart"], deg=m["deg"],
                     epart_min=min(m["epart"].values()), k5={k: list(v) for k, v in k5(edges, pat).items()},
                     edges=[dict(a=sorted(e)[0], b=sorted(e)[1], name=name.get(e, "신규"), hard=pairs[e]["hard"], dum=pairs[e]["dum"],
                                 cyc5=m["epart"][e], status=("유지" if e in base else "추가")) for e in sorted(edges, key=lambda e: (name.get(e, "Z"), tuple(sorted(e))))])
            for d in r["edges"]: rows.append([tag, d["name"], d["a"], d["b"], d["status"], "어려움" if d["hard"] else "쉬움", d["dum"], d["cyc5"]])
        else:
            r.update(add=None)
        r["before"] = dict(n_edges=len(base), ncyc4=m0["ncyc4"], npart=m0["npart"], deg=m0["deg"], epart_min=min(m0["epart"].values()) if m0["epart"] else None,
                           epart_lt2=sum(1 for v in m0["epart"].values() if v < 2), k5={k: list(v) for k, v in k5(base, pat).items()})
        result[tag] = r
    result["meta"] = dict(at=subprocess.run(["date", "-Is"], capture_output=True, text=True).stdout.strip(), script_sha256=sha(__file__), regkey=REGKEY,
                          nodes=nodes, pat=pat, flagged=flagged, criteria=dict(MAXDEG=MAXDEG, MAXADD=MAXADD, L_NODE=L_NODE, L_EDGE=L_EDGE))
    (out / "design.json").write_text(json.dumps(result, ensure_ascii=False, indent=1), encoding="utf-8")
    with open(out / "design_edges.tsv", "x", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t"); w.writerow(["판", "엣지", "A", "B", "상태", "난이도(cooh_fail)", "dummies_heavy", "길이≤5 고리 참여 수"]); w.writerows(rows)
    for tag in ("유지17", "유지14"):
        r = result[tag]
        print(f"== {tag}: 유지 {len(r['keep'])} · 현재 위반 {r['base_violations']} · 후보 {r['n_candidates']} (쉬움 {r['n_easy']} · 어려움 {r['n_hard']})")
        if r["add"] is None: print("   ★탐색 상한 안에 해 없음"); continue
        print(f"   추가 {r['n_add']} (어려움 {r['n_add_hard']}) · dummies 합 {r['dum_sum']} · 같은 점수의 해 {r['n_ties']}개 · 총 엣지 {r['n_edges']} · 길이≤4 고리 {r['ncyc4']}")
        for a in r["add"]: print(f"     + {a['a']} – {a['b']}  ({'어려움' if a['hard'] else '쉬움'} · 더미 {a['dum']}{' · 기존 ' + a['was'] if a['was'] else ''})")
        print("   노드별 길이≤4 고리 참여:", r["npart"]); print("   차수:", r["deg"]); print("   엣지 고리 참여 최소:", r["epart_min"], "· K5:", r["k5"])


if __name__ == "__main__":
    if len(sys.argv) >= 2 and sys.argv[1] == "--pretest": pretest()
    elif len(sys.argv) == 3 and sys.argv[1] == "--run": main(Path(sys.argv[2]).resolve())
    else: print(__doc__); sys.exit(2)
