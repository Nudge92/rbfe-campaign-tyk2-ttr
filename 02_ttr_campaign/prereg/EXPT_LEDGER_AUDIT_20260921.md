# 실험값 원장 감사 + 계 정의 후속 (2026-09-21)

★새 시뮬레이션 0 · 파드 기동 0 · 기존 파일만 읽음. ★못 찾은 것은 "미상".

---

## (1) Kd1 / Kd2 원장 — ★앵커 3종

원장 파일: **`/home/nudge/Project/CADD/ttr/qsar/data/interim/lit_validated_10_1021_jm0700159.csv`**
(65 레코드 · Kd **26건** = 13 화합물 × Kd1/Kd2 · `extracted_by` = ★**`claude_pdf_read`** — ★ChEMBL 덤프가 아니라 ★논문 PDF 판독본)

| 화합물 | 실제로 쓴 실험값 | ★Kd1 인가 Kd2 인가 | 출처 레코드 ID | 출처 파일 절대경로 |
|---|---:|---|---|---|
| **CHEMBL438498** | **0.2 nM** (pKd 9.70) | ★**Kd1** (`notes` = `Kd1 = 1차 결합자리`) | `Table 1` / `S.N. 6` / CSV 줄 **27** | `/home/nudge/Project/CADD/ttr/qsar/data/interim/lit_validated_10_1021_jm0700159.csv` |
| **CHEMBL241454** | **27 nM** (pKd 7.57) | ★**Kd1** (`notes` = `Kd1 = 1차 결합자리`) | `Table 1` / `S.N. 10` / CSV 줄 **47** | 〃 |
| **CHEMBL240808** | **27.0 nM** (pKd 7.57) | ★**Kd1** (`notes` = `Kd1 = 1차 결합자리`) | `Table 1` / `S.N. 7` / CSV 줄 **32** | 〃 |

### 같은 문서의 다른 쪽 값(Kd2)과 그 값으로 바꿨을 때의 오차

| 화합물 | Kd2 | 출처 레코드 ID | CSV 줄 |
|---|---:|---|---|
| CHEMBL438498 | **8300 nM** | `S.N. 6 (Kd2)` | 28 |
| CHEMBL241454 | **292.0 nM** | `S.N. 10 (Kd2)` | 48 |
| CHEMBL240808 | **3.0 nM** ★역전 | `S.N. 7 (Kd2)` | 33 |

ΔΔG(A→B) = RT·ln(Kd_B/Kd_A) · RT = 0.61603 kcal/mol @ 310 K

| 간선 | 계산 ΔΔG | 실험(Kd1) | ★오차(Kd1) | 실험(Kd2) | ★오차(Kd2) |
|---|---:|---:|---:|---:|---:|
| **N01** (241454→240808) | −3.867 | +0.000 | **−3.867** | −2.820 | **−1.047** |
| **N04** (240808→438498) | −0.490 | −3.022 | **+2.532** | +4.882 | **−5.372** |
| **N07** (241454→438498) | −5.109 | −3.022 | **−2.087** | +2.062 | **−7.171** |

### ★★한 줄 명시 — 종류 일치 여부

★★**셋 다 Kd1 이다. 섞여 있지 않다.** (원장 `notes` 열에 세 건 모두 `Kd1 = 1차 결합자리` 로 적혀 있다.)

### ★fep_nodes_smiles.csv ↔ 최종 보고 대조

| 경로 | 값 |
|---|---|
| `/home/nudge/Project/CADD/ttr/qsar/reports/fep_nodes_smiles.csv` | 438498 pKd **9.7** · 241454 **7.57** · 240808 **7.57** · `assay` 전부 `CHEMBL900184` |
| `/home/nudge/Project/CADD/ttr/qsar/src/final_report_v11.py` (24행) | `TRUE={"N01":(0.000,0.206),"N04":(-3.022,0.307),"N07":(-3.022,0.264)}` ★**하드코딩** |
| `/home/nudge/Project/CADD/ttr/qsar/src/edge_table_v11.py` (26행) | `TRUE={"N01":0.000,"N04":-3.022,"N07":-3.022}` ★하드코딩 |
| `/home/nudge/Project/CADD/ttr/qsar/src/node_potential_v11.py` (24행) | `EXP={"N01":0.000,"N04":-3.022,"N07":-3.022}` ★하드코딩 |

★**코드 추적 결과**: `fep_nodes_smiles.csv` 를 읽는 코드는 ★**`/home/nudge/Project/CADD/ttr/qsar/src/net_minus1.py` (22행)** 과 ★**`.../src/net_minus1_konn.py` (19행)** ★둘뿐이고, ★**최종 보고 3종은 이 파일을 읽지 않는다.**
★**값은 일치한다** (pKd 9.7 − 7.57 = 2.13 log × 1.41847 = 3.021 ≈ 하드코딩 3.022 · 차이는 pKd 반올림 대 원 Kd 사용). ★★**단 일치는 ★데이터 흐름이 아니라 ★손 입력의 결과다 — 파일을 고쳐도 보고는 안 바뀐다.**

---

## (2) 결합 자리 개수와 계산이 대응하는 양

| 항목 | 내용 | 출처 |
|---|---|---|
| 대표 계의 리간드 수 | ★**1 개** (`MOL` 블록 1 · 원자 28 · resid 118) | `/home/nudge/Project/CADD/ttr/qsar/data/fep_m1/pass1_full/N01_cplx_L3_r3/w00/prod.gro`<br>측정 `/home/nudge/Project/CADD/ttr/qsar/src/sysdef_audit_v11.py` |
| 어느 자리인가 | ★**site A** — Lys15 사본 2개가 2.02 / 3.45 Å, 나머지 2개는 16.41 / 16.57 Å → ★**두 T4 자리 중 1개만 점유** | 〃 |
| 설계 근거 | `자리 = site A (사슬 A+C가 형성)` | `/home/nudge/Project/CADD/ttr/pocket_design.md` (§0 표) |
| ★대응하는 양 | ★★**이 배치(2자리 중 1자리 점유)는 Kd1(1차 결합자리)에 대응한다.** | 위 3줄 |
| 통계 인자 | ★`−RT·ln2` = **0.411 kcal/mol** (298.15 K) / **0.427** (310 K · 우리 MD 온도 · 병기 규약 `ANCHOR_TRIANGLE_PREREG.md:45-48`). ★★**전 간선에 동일하게 들어가므로 ΔΔG 에서 상쇄된다.** | — |
| Adair 식 | ★**2자리 협동 모형** — `fluorescence titration · Adair equation (two binding sites)` (`assay_protocol` 열 · 26/26 레코드 동일) | `/home/nudge/Project/CADD/ttr/qsar/data/interim/lit_validated_10_1021_jm0700159.csv` |
| 실험 문서가 Kd1/Kd2 를 분리 보고하는가 | ★★**한다.** 한 화합물당 2 레코드(`S.N. n` / `S.N. n (Kd2)`) · 13 화합물 × 2 = 26 건 | 〃 |

★협동성 크기는 추정하지 않았다. ★문서에 있는 숫자만 옮겼다.

---

## (3) 6E72 홀로 여부 — ★파일로

### ㄱ. `6E72_tetramer_1model.pdb` 의 레코드 — ★전문 인용

```
$ grep -E "^(HEADER|TITLE|COMPND|SOURCE|EXPDTA|HET |HETNAM|HETSYN|FORMUL)" \
    /home/nudge/Project/CADD/ttr/data/reference/6E72_tetramer_1model.pdb
(출력 없음)

$ awk '{print substr($0,1,6)}' ... | sort | uniq -c
   3632 ATOM
     36 HETATM
      1 END

$ awk '/^HETATM/{print substr($0,18,3)}' ... | sort | uniq -c
     36 OCS
```

★★**TITLE · COMPND · HET · HETNAM · FORMUL 레코드가 ★하나도 없다.** ★HETATM 은 ★`OCS`(시스테인설폰산 · 번역후 변형 잔기) ★36 원자뿐이고 ★**리간드는 들어 있지 않다.**

### ㄴ. 결합해 있던 리간드 — ★원본 CIF 원문

파일 `/home/nudge/Project/CADD/ttr/pairtest/pdb/6E72.cif`

```
_struct.title    'Structure of Human Transthyretin Val30Met Mutant in Complex with Tafamidis'

_pdbx_entity_nonpoly.name                                     _pdbx_entity_nonpoly.comp_id
2 '2-(3,5-dichlorophenyl)-1,3-benzoxazole-6-carboxylic acid'   3MI
3 water                                                       HOH

_chem_comp:
3MI non-polymer . '2-(3,5-dichlorophenyl)-1,3-benzoxazole-6-carboxylic acid' Tafamidis 'C14 H7 Cl2 N O3' 308.116
OCS 'L-peptide linking' n 'CYSTEINESULFONIC ACID' ? 'C3 H7 N O5 S' 169.156
```

★★**결합해 있던 리간드 = `3MI` = ★타파미디스.** ★제목 자체가 ★`Val30Met Mutant in Complex with Tafamidis` 다.

### ㄷ. 어느 단계에서 제거됐는가 — ★스크립트 경로

| 단계 | 상태 |
|---|---|
| 제작 내용 기술 | `/home/nudge/Project/CADD/ttr/pocket_design.md` (43행) — ★`사슬 재명명 A1→A · B1→B · A2→C · B2→D, 물 제거, 원자번호 재부여` |
| ★**타파미디스 제거 언급** | ★★**그 문장에 없다.** 물 제거만 적혀 있다. 그런데 파일에는 3MI 원자가 0개다 |
| ★**제작 스크립트 경로** | ★★**미상** — `tetramer_1model` 문자열을 ★쓰는(write) 코드를 `.py`·`.sh` 전수 검색했으나 ★0건. 읽는 코드만 있다 |
| 추출된 리간드 파일 | `/home/nudge/Project/CADD/ttr/data/reference/6E72_3MI_siteA_orient1.pdb` (20원자 · 2026-09-06 14:14)<br>`/home/nudge/Project/CADD/ttr/docking_validation/ligands/3MI_crystal_orient2.sdf` (2026-09-06 14:48) |

### ㄹ. 우리 리간드 배치가 무엇을 기준으로 정해졌는가 — ★코드 원문

파일 `/home/nudge/Project/CADD/ttr/pairtest/constrained.py`

```python
13: REF_DEPTH=5.89; SHIFT=-2.05
16: # ── 채널축 (6E72 사량체 site A) ──
17: st=PDBParser(...).get_structure("r", str(T/"data/reference/6E72_tetramer_1model.pdb"))
18: K=np.array([r["NZ"].coord ... if r.id[1]==15 ...])        # Lys15 NZ
19: S=np.array([r["OG"].coord ... if r.id[1]==117 ...])       # Ser117 OG
20: PROT=np.array([a.coord for ch in st[0] for r in ch for a in r if a.element!="H"])
23: REF=loadsdf(LIG/"3MI_crystal_orient2.sdf")   # ★기준 배향 2
64: FREF=frame(REF)
72:     R,t=kabsch(F,FREF); P=(R@P0.T).T+t
76:     tgt = REF_DEPTH if mode=="a" else (depth(P0)+SHIFT)
```

★배치의 ★**모든** 기준이 홀로 구조에서 나온다:
- **프레임** `FREF` = ★타파미디스 결정 포즈의 고리계 프레임
- **채널축** `AX` = ★홀로 사량체 site A 의 Lys15 NZ → Ser117 OG
- **깊이 목표** `REF_DEPTH = 5.89` = ★타파미디스의 깊이 (★유도 과정은 파일에 없음 — 하드코딩)
- **충돌 채점** `PROT` = ★홀로 단백질 좌표

### ★★판정

★★**순환논법이다.** ★배치 기준(프레임·축·깊이·충돌)이 ★전부 ★타파미디스가 결합한 홀로 포켓 ★하나에서 나왔고, ★그 배치를 다시 ★같은 홀로 포켓 안에서 채점한다.
★★**따라서 ★포켓 모양이 우리 디페닐에터 리간드에 대해 틀렸더라도 ★이 절차는 그것을 검출할 수 없다.** ★"디페닐에터가 타파미디스처럼 결합한다"는 것은 ★검증된 것이 아니라 ★**가정으로 들어가 있다.**
★변호하지 않는다.

---

## (4) TRUE σ 0.206 / 0.307 / 0.264 의 출처

| 항목 | 결과 |
|---|---|
| 전수 검색 범위 | `docs/` · `src/` · `analysis/` · `data/interim/` · `reports/` |
| 결과 | ★★**출처 없음.** 걸린 문자열은 전부 무관한 맥락(`TABLE5_A977.md:20` 0.2067 · `WINDOW_CANDIDATE_RESULT.md:103` 0.2063 · `COLLECTION_PLAN.md:390` CI 상한 · `FINDINGS.md:463` p값 0.307 · `TOX24_ASSESSMENT.md` Tanimoto 0.2642) |
| ★등재 | ★**RISK-6 · "출처 미상 하드코딩"** (`/home/nudge/Project/CADD/ttr/qsar/docs/RISK_REGISTER.md`) |

### ★앵커 삼각형 0.65σ 판정이 이 세 숫자에 의존하는가 — ★코드·산술로 확인

```
문서에 적힌 값 : |닫힘| / σ_prop = 0.751434 / 1.149987 = 0.65 σ
   → /home/nudge/Project/CADD/ttr/qsar/docs/stage1_out/pass1/EDGE_TABLE_20260920.md:73

계산 σ̂ 전파   = sqrt(0.30805² + 1.08134² + 0.24140²) = ★1.149987   ← ★정확히 일치
하드코딩 σ 전파 = sqrt(0.206²  + 0.307²   + 0.264²)   =  0.454292   ← 불일치
```

★코드 확인: `/home/nudge/Project/CADD/ttr/qsar/src/final_report_v11.py` 에서 `TRUE` 가 나오는 줄은 ★**24행(정의)과 76행(`if e in TRUE:`) 뿐**이다 — ★**튜플의 둘째 원소(σ)는 한 번도 읽히지 않는다.**

★★**판정: 0.65σ 판정은 ★이 세 숫자에 ★의존하지 않는다.** 계산된 회차간 σ̂ 를 쓴다. ★세 숫자는 ★**읽히지 않는 죽은 값**이다.
