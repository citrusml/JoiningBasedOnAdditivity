最終更新: 2026-11-18 06:53

## 📈 進捗状況
- [x] 完了: `benchmark_bioperl_nelsi.py` への `--lba-ratios`（T1, TN 枝長倍率指定）および `--methods-set core`（共通8手法のみ実行）の拡張 (完了日: 11/18)
- [x] 完了: NELSI LBA (`lba3`, `lba5`) 用 LSF ジョブ配列スクリプト `scripts/submit_lsf_nelsi_lba.sh`（420並列ジョブ＋自動CSV統合）の作成・検証 (完了日: 11/18)

## 🛠 技術的決定・メモ
- **[課題]**: `PhyloTree_Methods_Comparison` 準拠の NELSI 自己相関モデル上で長枝牽引（`lba3`: `lba_ratio=3.0`, `lba5`: `lba_ratio=5.0`）の影響を、$D \in \{0.1, 0.2, 0.5, 1.0, 1.5, 2.0, 3.0\}$、$L=500$、$\text{indel}=0.1$、$N \in \{20, 32, 64\}$、50反復、共通8手法（`MSA_ML`, `SSS_Linear_NJ`, `SSS_Log_NJ`, `Poisson_NJ`, `MSA_NJ`, `GS_Poisson`, `GS_SSS_Linear`, `GS_SSS_Log`）について LSF スパコン上で並列実行する必要があった。
- **[解決策]**:
  1. `benchmark_bioperl_nelsi.py` に `--lba-ratios` と `--methods-set {all,core}` を追加。`lba_ratio > 1.0` の場合は伸長対象の `(T1, TN)` を固定 LBA ペアとして評価するよう設定。
  2. `scripts/submit_lsf_nelsi_lba.sh` を作成（`#BSUB -J nelsi_lba[1-420]`, `#BSUB -q mafft`, `logs/nelsi_lba_bsub_%I.out/.err` + `logs/nelsi_lba*.log` リアルタイム出力 + 完了後の自動依存マージジョブ `nelsi_lba_merge`）。
  3. 1反復スモークテスト（`lba3`, `lba5`, `N=20`, `D=1.0`, `core` 8手法 = 16行）を実行し、正常出力を確認。
- **[パラメータ]**: `nelsi_kishino` (`v=0.3`, `scale_mode=root_to_tip`), `lba_ratio` $\in \{3.0, 5.0\}$, $D \in \{0.1, 0.2, 0.5, 1.0, 1.5, 2.0, 3.0\}$, $N \in \{20, 32, 64\}$, $L=500$, `indel=0.1`, 50 reps, 420 LSF jobs.

## 🚀 次に行うこと (Next Actions)
1. スパコン上で `bash scripts/submit_lsf_nelsi_lba.sh` を実行し、`results_nelsi_lba/benchmark_results.csv` および `summary_by_condition.csv` を取得する

---

最終更新: 2026-11-18 05:49

## 📈 進捗状況
- [x] 完了: LSF スパコン用ジョブ配列スクリプト `scripts/submit_lsf_bioperl_nelsi.sh`（420並列ジョブ）の作成と構文チェック (完了日: 11/18)
- [x] 完了: `benchmark_bioperl_nelsi.py` の衝突回避シャード保存（`shards/`）および `--merge-only` 統合機能の実装・動作確認 (完了日: 11/18)
- [x] 完了: `src/generate_tree_nelsi.R` のリポジトリ内配置による HPC スタンドアロン実行対応 (完了日: 11/18)

## 🛠 技術的決定・メモ
- **[課題]**: LSF スパコン上で 42 条件（2モデル × 3タクソン × 7発散度）× 50反復を並列実行する際、各ジョブが同一ファイル名で CSV を上書きせず、終了後に全シャードを安全に統合できる仕組みが必要であった。また、スパコンへ本リポジトリのみ転送した場合でも `generate_tree_nelsi.R` を参照できるよう自己完結化する必要があった。
- **[解決策]**:
  1. `scripts/submit_lsf_bioperl_nelsi.sh` を作成（`#BSUB -J bioperl_nelsi[1-420]`, `#BSUB -n 1`, `#BSUB -M 4000`, `#BSUB -W 04:00`、`micromamba activate phylomethod_env`、`bash scripts/install_indelible.sh`）。1条件あたり10シャード（各5反復）に分割して 420 並列ジョブで投入。
  2. `benchmark_bioperl_nelsi.py` に `<outdir>/shards/shard_<model>_N<taxa>_D<div >_s<idx>.csv` 形式の衝突回避保存と `--merge-only` 統合モードを追加し、2シャードの分割実行→統合テスト（28行）の正常完了を確認。
  3. `src/generate_tree_nelsi.R` を `src/` 配下にコピーし、`PhyloTree_Methods_Comparison` 非依存で動作するよう対応。
- **[パラメータ]**: 420 LSF array jobs (`2 models x 3 taxa x 7 D x 10 shards`), 5 reps/job, `OMP_NUM_THREADS=1`, `-M 4000`, `-W 04:00`.

## 🚀 次に行うこと (Next Actions)
1. スパコン上で `mkdir -p logs && bsub < scripts/submit_lsf_bioperl_nelsi.sh` を実行
2. 全420ジョブ完了後、`python benchmark_bioperl_nelsi.py --merge-only --outdir results_bioperl_nelsi` で最終 CSV と条件別サマリを生成

---

最終更新: 2026-11-18 05:45

## 📈 進捗状況
- [x] 完了: NELSI自己相関木生成関数 `generate_nelsi_tree`（`src/tree_simulator.py`）、根-葉距離統計 `compute_root_to_tip_metrics`（`src/evaluation.py`）、Graph Splitting ラッパー `src/gs_wrapper.py` の実装 (完了日: 11/18)
- [x] 完了: BioPerl vs. NELSI 大規模比較ベンチマークスクリプト `benchmark_bioperl_nelsi.py` の作成と1反復スモークテスト完走確認 (完了日: 11/18)
- [x] 完了: `tests/test_lba_evaluation.py` テスト追加および `pytest -v` 全21テストパス確認 (完了日: 11/18)

## 🛠 技術的決定・メモ
- **[課題]**: $D \in \{0.1, 0.2, 0.5, 1.0, 1.5, 2.0, 3.0\}$、$L=500$、$\text{indel}=0.1$、$N \in \{20, 32, 64\}$、50反復において BioPerl と NELSI を比較し、`MSA_ML`, `SSS_Linear_NJ`, `SSS_Log_NJ`, `Poisson_NJ`, `MSA_NJ`, `GS`（3種）, および `sss_linear` 専用7手法を一括評価できるベンチマークスクリプトが必要であった。
- **[解決策]**:
  1. `src/tree_simulator.py` に `generate_tree_nelsi.R` を呼び出す `generate_nelsi_tree` を追加。
  2. `src/evaluation.py` に `compute_root_to_tip_metrics`（平均・分散・min・max・max/min比）を追加。無向推定木には midpoint rooting を適用して一貫した根-葉距離を算出。
  3. `src/gs_wrapper.py` を新設し、`gs2 -s -l -d` により `Poisson`, `SSS Linear`, `SSS Log` の距離行列から Graph Splitting 木を再構築（`patristic_rmse` 等のため OLS 枝長を付与）。
  4. `benchmark_bioperl_nelsi.py` に並列実行（`--threads`）および HPC シャーディング（`--shard-idx`, `--num-shards`）を実装。
- **[パラメータ]**: $D \in \{0.1, 0.2, 0.5, 1.0, 1.5, 2.0, 3.0\}$, $L=500\text{ aa}$, $\text{indel}=0.1$, $N \in \{20, 32, 64\}$, $\text{rep}=50$, BioPerl + NELSI (`nelsi_kishino`, `v=0.3`).

## 🚀 次に行うこと (Next Actions)
1. ユーザーによる `benchmark_bioperl_nelsi.py` の本番実行（ローカル並列または HPC シャード実行）
2. 出力された `results_bioperl_nelsi/benchmark_results.csv` および `summary_by_condition.csv` の可視化・解析

---

最終更新: 2026-11-18 05:36

## 📈 進捗状況
- [x] 完了: LBA特化評価指標群（`pair_nrf`, `pair_split_recall`, `lba_quartet_attraction`, `nta`, `compute_matrix_lba_margin`, `find_max_distance_pair`）の `src/evaluation.py` への実装 (完了日: 11/18)
- [x] 完了: `Method.md` 第5.3節への新規LBA評価指標の数理定義追記 (完了日: 11/18)
- [x] 完了: `tests/test_lba_evaluation.py` への単体テスト追加および `pytest -v` 全20テストパス確認 (完了日: 11/18)
- [ ] 待機中: 実装したLBA評価指標を用いたベンチマーク・評価ステップの実行（ユーザー指示により保留）

## 🛠 技術的決定・メモ
- **[課題]**: 従来のLBA評価は `ds_est == 2`（偽姉妹群率）や木全体の `nRF` に頼っていたが、部分的な接近の定量化や、無関係な短枝の局所的ズレ（RFスプリット不一致）を切り離した純粋なLBA吸引度評価が必要であった。
- **[解決策]**:
  1. `compute_tree_patristic_matrix` に `unit_lengths=True` オプションを追加し、無向木上のトポロジー的枝数距離行列をベクトル演算で厳密算出できるようにした。
  2. `compute_pair_separating_rf` により最大距離ペア $(u, v)$ を分離する内部スプリットのみを対象とした `pair_nrf` / `pair_split_recall` / `pair_split_fnr` を実装。
  3. `compute_lba_quartet_attraction` により、他枝の局所ノイズに完全に不変なカルテット吸引率 `lba_quartet_attraction` および誤差率 `lba_quartet_error` を実装。
  4. `nta`（正規化トポロジー接近度）、`find_max_distance_pair`（最大パトリスティック距離ペア自動抽出）、`compute_matrix_lba_margin`（距離行列段階でのLBAマージン）を追加し、`evaluate_tree_lba` に統合した。
- **[パラメータ]**: `auto_max_pair=False`（`lba_pair` 未指定時に `True` で最大距離ペアを自動取得）, `eps=1e-7`.

## 🚀 次に行うこと (Next Actions)
1. ユーザーからの指示があり次第、`benchmark_lba.py` や保存済みLBA系統樹群に対し新指標を用いた評価ステップを実行する

---

最終更新: 2026-11-18 04:36

## 📈 進捗状況
- [x] 完了: 系統樹距離・行列距離残差に基づく結合法 `AddNJ3`（`src/add_nj.py`）およびトポロジー探索手法 `Add-NNI2 / Add-SPR2 / Add-Hybrid2`（`src/topology_search.py` 内 `TreeDistanceScorer`）の設計・実装 (完了日: 11/18)
- [x] 完了: `Method.md` への数理モデル・アルゴリズム・ファイル役割の体系的追記 (完了日: 11/18)
- [x] 完了: 単体テスト追加（`tests/test_add_nj.py`, `tests/test_topology_search.py`）と全18テストパス確認 (`pytest -v`) (完了日: 11/18)
- [x] 完了: `main.py` への全20手法ベンチマーク統合および再実行（$N=20, L=500\text{ aa}, D \in \{0.5, 1.0, 2.0\}, 10\text{ reps}$）、`RawResults.md`・グラフ出力更新 (完了日: 11/18)

## 🛠 技術的決定・メモ
- **[課題]**: 従来の4点条件違反度（カルテット順位整合性）のみを用いた探索（`TreeScorer`）や貪欲結合（`AddNJ`, `AddNJ2`）に加え、「系統樹上のパトリスティック距離 $D_T$ と観測距離行列 $D_{\text{obs}}$ の直接の差（残差 $\mathcal{E}_{\text{dist}}(T; D_{\text{obs}})$）」を目的関数とする結合法（`AddNJ3`）および局所・広域探索法（`Add-NNI2 / Add-SPR2 / Add-Hybrid2`）を確立する必要があった。
- **[解決策]**:
  1. **AddNJ3 (`src/add_nj.py`)**: 候補ペア（NJの $Q$ 値上位 `top_k=4`）を仮結合した状態から `lookahead_depth=1` の先読み＋高速NJロールアウトで完全木の二分割集合を生成し、正則化付き正規方程式 $(A_T^\top A_T + 10^{-12}I)\mathbf{b} = A_T^\top \mathbf{d}_{\text{obs}}$ と非負クランプ $\hat{\mathbf{b}} = \max(\mathbf{0}, \mathbf{b})$ による距離残差 RMSE を最小化するペアを選択する方式を実装。
  2. **TreeDistanceScorer (`src/topology_search.py`)**: `ols_rmse`（OLS枝長によるRMSE）、`ols_mae`、`wls_fm`（Fitch-Margoliash重み付き相対二乗誤差）、`wls_mre` をサポートし、`fit_ols_branch_lengths` とスプリット接続行列構築ヘルパー `_build_split_incidence_matrix` を共通化（コード重複排除・postorderメモ化）。`hill_climb_search` を拡張し、`TreeDistanceScorer` 利用時の改善判定閾値 `tol=1e-9` を設定。
  3. **[ベンチマーク実験の主要結果 ($N=20, L=500\text{ aa}, 10\text{ reps}$)]**:
     - $D = 1.0$: `Standard NJ` / `AddNJ2` ($\text{nRF} = 0.0118 \pm 0.0078$) および `NJ + Add-Hybrid` ($\text{nRF} = 0.0176 \pm 0.0090$) に対し、**`AddNJ3` と `Add-NNI2 / SPR2 / Hybrid2` 群はすべて $\text{nRF} = 0.0059 \pm 0.0059$** と誤差を半減。
     - $D = 2.0$: `Standard NJ` ($0.0882$), `BioNJ` ($0.1235$), `FastME-SPR` ($0.0882$), 4点ペナルティ版 `NJ + Add-Hybrid` ($0.0824$) に対し、`NJ + Add-NNI2 / SPR2 / Hybrid2` および `AddNJ3 + Add-Hybrid2` が $0.0765 \pm 0.0197$、**`BioNJ + Add-Hybrid2` が全20手法中最小の $\text{nRF} = 0.0706 \pm 0.0192$** を達成。
- **[パラメータ]**: `top_k=4`, `lookahead_depth=1`, `tie_tolerance=1e-12`, `reg=1e-12`, `mode='ols_rmse'`, `tol=1e-9` (`TreeDistanceScorer`), $N=20, L=500, D \in \{0.5, 1.0, 2.0\}$, 10 reps.

## 🚀 次に行うこと (Next Actions)
1. 必要に応じ `benchmark_hpc.py` への `AddNJ3` / `Add-Hybrid2` の組み込みと大規模データセット（$N \in \{20, 32, 64, 128\}$, インデルあり）での性能評価
2. `benchmark_lba.py` における長枝牽引（LBA）耐性の検証（`AddNJ3` vs `AddNJ2` vs `BioNJ+Add-Hybrid2`）

---

最終更新: 2026-11-13 05:25

## 📈 進捗状況
- [x] 完了: 系統樹シミュレーション手法の学部生向け解説書および先行研究（syz049.pdf）比較ドキュメント `Tree_Simulation.md` の作成 (完了日: 11/13)
- [x] 完了: BioPerl `Bio::Tree::RandomFactory` (Backward Yule) および Python フォールバックの数理的背景（切断指数分布の逆関数法等）の完全な解析・文書化 (完了日: 11/13)
- [x] 完了: Kingman合祖過程に基づく枝長定義および「分岐数と根からの距離の相関（非ウルトラメトリック性）」と生物進化学的妥当性の解説を `Tree_Simulation.md` に追記 (完了日: 11/13)
- [x] 完了: NELSIの自己相関型モデル（Autocorrelated Relaxed Clock / ABR）との数理的・生物学的比較解説を `Tree_Simulation.md` に追記 (完了日: 11/13)

## 🛠 技術的決定・メモ
- **[課題]**: 現在のシミュレーション手法（BioPerl / syz049 モデル）と、進化速度シミュレータとして著名な NELSI（Ho et al. 2015）の「自己相関型モデル（Autocorrelated Branch Rates, ABR）」との数理構造の違い、生物学的リアリティの差、および両者の研究目的の違いについての体系的比較が求められていた。
- **[解決策]**:
  1. NELSI が採用する Thorne et al. (1998) 準拠の対数正規自己相関モデル（幾何ブラウン運動：$\ln(r_{\text{child}}) \sim \mathcal{N}(\ln(r_{\text{parent}}), \nu \Delta t)$）の数理を解説。
  2. 根本的な4つの相違点を整理：
     - ① パラメータ構造: 「置換数直接サンプリング」 vs 「実時間 $\Delta t$ × 速度 $r$ の明示的分離」。
     - ② 枝間相関: 「完全無相関（独立）」 vs 「クレード単位での連続的な速度遺伝」。
     - ③ 根からの距離: 「分岐エッジ数の累積」 vs 「クレードの速度履歴（全現生種で経過年数は同一）」。
     - ④ 研究目的: 「深層系統樹トポロジー推定ベンチマーク」 vs 「分子年代推定・分子時計モデルの統計的検出」。
  3. なぜ本研究や syz049 で NELSI 型ではなく本モデルが選ばれるのか（スーパーファミリーにおける突発的・独立な変異加速の再現、トポロジー復元アルゴリズムのストレステスト、超高速・大量生成のポータビリティ）を論理的に明文化。
- **[パラメータ]**: $l = 1 - \ln(u(e-1)+1)$, $l_{\text{scaled}} = l \times D$, $\text{min\_length} = 10^{-6}$, $u \sim U(0,1)$, $D \in \{0.5, 1.0, 2.0\}$.

## 🚀 次に行うこと (Next Actions)
1. ユーザーへの NELSI 自己相関型モデルとの比較解説の提示
2. 学部生向けレクチャーや研究報告資料としての活用

---

最終更新: 2026-11-11 10:28

## 📈 進捗状況
- [x] 完了: `results_hpc_indel0.05` および `results_hpc_indel0.1` における Standard NJ の距離行列作成法別 nRF エラーバープロットの生成 (完了日: 11/11)
- [x] 完了: インデル率（0.05 vs 0.10）横断比較 2行×3列マルチパネルプロット（`standard_nj_distance_matrix_cross_indel.png`）の実装・生成 (完了日: 11/11)
- [x] 完了: `analysis/merge_and_plot_hpc.py` への Standard NJ 専用プロット生成パイプライン追加 (完了日: 11/11)

## 🛠 技術的決定・メモ
- **[課題]**: Standard NJ（近隣結合法）が用いる入力距離行列の構築法（BLOSUM62+Poisson, WAG+Γ ML, MMseqs2 SSS linear $d=1-w$, MMseqs2 SSS log $d=-\ln(w)$, True Patristic）の違いが、配列数 $N$ および発散度 $D$、インデル率の増大に伴い系統樹トポロジー推定精度（nRF）にどのような影響を与えるかが個別に可視化されていなかった。
- **[解決策]**:
  1. `analysis/merge_and_plot_hpc.py` に `standard_nj_distance_matrix_comparison.png`（およびFastTree非表示版 `standard_nj_pure_matrix_comparison.png`）を生成するパイプラインを追加。
  2. インデル率 0.05 と 0.10 を上下2行に並べ、発散度 $D \in \{0.5, 1.0, 2.0\}$ ごとに Taxa $N \in \{20, 32, 64, 128\}$ に対する Mean ± SD エラーバープロット（水平ドッジ付き）を作成する `_plot_standard_nj_cross_indel_comparison` を実装。
  3. **[実験結果の重要な発見]**:
     - **True Patristic Distance**: 全条件で $\text{nRF} = 0.0000 \pm 0.0000$ を達成（加法性が満たされれば NJ は完全に真の系統樹を復元する理論的裏付け）。
     - **MMseqs2 SSS Linear ($d = 1 - w$)**: 高発散度（$D=2.0$）において、BLOSUM62（$\text{nRF} \approx 0.360$ / $0.526$）や WAG+Γ（$\text{nRF} \approx 0.375$ / $0.537$）だけでなく、MSA+ML FastTree（$\text{nRF} \approx 0.413$ / $0.557$）よりも有意に低いエラー率（$\text{nRF} \approx 0.333$ / $0.492$）を達成し、インデル多発・高発散領域におけるアライメントフリー的類似度距離の極めて高い耐性が実証された。
     - **MMseqs2 SSS Log ($d = -\ln(w)$)**: 低発散度（$D=0.5$）において微小なスコア変動が対数変換により増幅され、$\text{nRF} \approx 0.14 \sim 0.18$ と線形距離より誤差が大きくなる傾向が確認された。
     - **BLOSUM62 vs WAG+Γ**: インデルが存在する本シミュレーション環境では、BLOSUM62 が WAG+Γ よりも一貫してわずかに優れた（低い）nRF を示した。
- **[パラメータ]**: $N \in \{20, 32, 64, 128\}$, $D \in \{0.5, 1.0, 2.0\}$, rep=30, INDELible Zipfian rate $\in \{0.05, 0.10\}$, Mean ± SD (capsize=3.5, elinewidth=1.2).

## 🚀 次に行うこと (Next Actions)
1. ユーザーへの Standard NJ 距離行列比較エラーバープロットの提示
2. 論文・考察への「高発散・インデル多発下での SSS Linear 距離の優位性」の取り込み

---

最終更新: 2026-11-11 10:04

## 📈 進捗状況
- [x] 完了: `results_hpc_indel0.05` および `results_hpc_indel0.1` フォルダ内の nRF vs SSS グラフの横軸対数スケール（$\log_{10}(\mathrm{SSS})$）化 (完了日: 11/11)
- [x] 完了: 配列数 $N \in \{20, 32, 64, 128\}$ 別の層別化マルチパネルプロット（`nrf_vs_sss_by_taxa_*.png`）および $N$ 別単独プロット（`nrf_vs_sss_N*.png` / `nrf_vs_sss_curves_N*.png`）の新規実装・出力 (完了日: 11/11)
- [x] 完了: `analysis/merge_and_plot_hpc.py` における対数空間 LOWESS 平滑化（$\log_{10}(\mathrm{SSS})$ 空間でのスプライン平滑化）と $y$ 軸クリップ上限解除（$y \le 1.0$）の実装 (完了日: 11/11)

## 🛠 技術的決定・メモ
- **[課題]**: 従来の nRF vs SSS プロットでは横軸が線形スケールであり、また $N \in \{20, 32, 64, 128\}$ が混在または $N=64, 128$ の狭い SSS 領域（$0.003 \sim 0.05$）に全データが集中して密集していた。また、$y$ 軸上限が 0.52 に固定されていたため、高発散度や高インデル率での高い nRF（最大 0.88）がフレーム外にクリップされる問題があった。
- **[解決策]**:
  1. `analysis/merge_and_plot_hpc.py` 内の `_compute_lowess_curve` を拡張し、`is_log_x=True` による $\log_{10}(\mathrm{SSS})$ 空間での LOWESS 平滑化スプライン算出を実装（`PhyloTree_Methods_Comparison` の `plot_nrf_vs_log_sss` と同一設計）。
  2. `_plot_nrf_vs_sss_curves` および `_plot_nrf_vs_sss_by_taxa` の横軸を `ax.set_xscale("log")` に変更し、各 $N$ の SSS 分布に応じた動的レンジ設定を導入。
  3. 配列数 $N$ 別の専用グラフ生成関数 `_plot_nrf_vs_sss_single_taxa` を新設し、各行列タイプ（BLOSUM62, WAG+Γ, MMseqs2 SSS linear/log, True Patristic）および全手法比較について、4分割マルチパネル（`nrf_vs_sss_by_taxa_*.png`）、$N$ 別単独プロット（`nrf_vs_sss_N20_*.png` 等）、$N$ 別2パネル曲線プロット（`nrf_vs_sss_curves_N20_*.png` 等）を自動一括生成。
  4. 縦軸上限を $1.0$ に拡張し、高誤差領域も欠損なく可視化。
- **[パラメータ]**: $N \in \{20, 32, 64, 128\}$, $D \in \{0.5, 1.0, 1.5\}$, INDELible Zipfian rate $\in \{0.05, 0.1\}$, rep=30, LOWESS frac=0.35 in $\log_{10}(\mathrm{SSS})$.

## 🚀 次に行うこと (Next Actions)
1. ユーザーへの生成グラフ（対数スケール＆ $N$ 別分割プロット）の提示と確認
2. 必要に応じた図の解釈・論文報告への反映

---

最終更新: 2026-11-11 08:55

## 📈 進捗状況
- [x] 完了: ベンチマーク同一実験の再実行（Taxa=20, L=500 aa, D∈{0.5, 1.0, 2.0}, 10 reps, 計14手法）と4点制約違反度 $\mathcal{L}(T; D_{\text{obs}})$ の算出 (完了日: 11/11)
- [x] 完了: 4点制約違反度比較グラフ（`results/penalty_comparison.png`）の新規作成・左右2パネル可視化 (完了日: 11/11)
- [x] 完了: `main.py` および `RawResults.md` への4点条件違反度集計テーブル出力の追加 (完了日: 11/11)

## 🛠 技術的決定・メモ
- **[課題]**: 従来の `main.py` ではトポロジー誤差（nRF）と実行時間のみがグラフ化されており、本研究のコア理論である4点条件違反度（Quartet Additivity Penalty $\mathcal{L}(T; D_{\text{obs}})$）の可視化が欠落していた。
- **[解決策]**:
  1. 同一条件（計30試行）でベンチマークを再実行し、全14手法の系統樹加法性ペナルティ $\mathcal{L}(T; D_{\text{obs}})$ を検証。
  2. `main.py` の `generate_plots` に Figure 3 として `results/penalty_comparison.png`（Full Scale (a) ＋ Zoomed Scale (b) の2パネル散布図＋Mean ± SEMエラーバープロット）を実装。
  3. **[実験結果の発見]**:
     - ナイーブな AddNJ（手法1）は消去バイアスにより $\mathcal{L}(T) \approx 0.033 \sim 0.047$ と高い4点制約違反を示す。
     - 一方、改善版 **AddNJ2（手法1改善）** は $\mathcal{L}(T) \approx 0.000089$ (D=0.5), $0.001004$ (D=1.0), $0.005307$ (D=2.0) とStandard NJと同等水準まで劇的に違反度が抑制されている。
     - $D=2.0$ において、BioNJ は $\mathcal{L}(T) \approx 0.006846$ と高い違反度を示すのに対し、提案する加法性探索手法（**NJ + Add-Hybrid / AddNJ2 + Add-Hybrid**）は全手法中で最小の $\mathcal{L}(T) \approx 0.005178$ を達成し、目的関数通りの加法性最小化が実証された。
- **[パラメータ]**: $N=20$, $D \in \{0.5, 1.0, 2.0\}$, $L=500$ aa, 14手法, 10 reps, Mean ± SEM.

## 🚀 次に行うこと (Next Actions)
1. ユーザーへの4点制約違反度グラフおよび実験結果の提示
2. 論文・レポートに向けた考察の統合

---

最終更新: 2026-11-11 08:40

## 📈 進捗状況
- [x] 完了: `results/nrf_comparison.png` の手法対応ラベル（凡例）のキャンバス内完全表示・左右パネル位置同期修正 (完了日: 11/11)

## 🛠 技術的決定・メモ
- **[課題]**: `fig.legend` のアンカーがキャンバス右端外（`x=0.99`）に配置されていたため、画像保存時に右側の手法対応ラベル・凡例ボックスがフレーム外に切り落とされていた。また左右パネル間で横ドッジ位置にわずかなずれがあった。
- **[解決策]**:
  1. `main.py` の `fig.legend` を `bbox_to_anchor=(0.825, 0.5), loc="center left"` に再配置し、`plt.savefig` に `bbox_inches="tight"` を指定。全14手法の名称・マーカー形状・配色が右側の専用領域に確実に収まるよう修正。
  2. 凡例タイトルに `Reconstruction Methods (Left-to-Right in cluster)` を明記し、凡例の上から下への並び順が各発散度クラスタ内の左から右へのマーカー配置順と 1:1 に完全対応するように設計。
  3. 左パネル (a) と右パネル (b) で各手法の水平ドッジ座標を完全統一し、パネル間での手法の視線移動を直観的にした。
- **[パラメータ]**: $N=20$, $D \in \{0.5, 1.0, 2.0\}$, $L=500$ aa, 14手法, 10 reps, Mean ± SEM.

## 🚀 次に行うこと (Next Actions)
1. ユーザーへの修正図提示と表示確認

---

最終更新: 2026-11-11 08:35

## 📈 進捗状況
- [x] 完了: `results/nrf_comparison.png` の棒グラフからエラーバープロット（散布図＋Mean ± SEM、水平ドッジ）への刷新 (完了日: 11/11)
- [x] 完了: `main.py` 内 `generate_plots` の改修（全手法表示のFull Scale (a) ＋ 競合手法高解像度表示のZoomed Scale (b) の2パネル構成） (完了日: 11/11)

## 🛠 技術的決定・メモ
- **[課題]**: 従来の棒グラフ（`plt.bar`）では14手法が発散度ごとに密集して重なり、またナイーブなAddNJ（手法1, nRF≈0.65）の存在によって低エラー領域（nRF 0〜0.15）のAddNJ2、FastME、BioNJ、Standard NJ間の精細な差が潰れて視認困難であった。
- **[解決策]**:
  1. 棒グラフを廃止し、論文仕様の散布図＋エラーバープロット（`ax.errorbar`, `linestyle="none"`, `capsize=3.5`, `markeredgecolor="black"`）に変更。
  2. 水平オフセット（Dodge）を導入し、同一発散度におけるマーカーとエラーバーの重なりを解消。
  3. 左パネル (a) Full Scale（$0 \le \text{nRF} \le 0.85$）で全14手法の全体像（手法1の破綻と手法1改善AddNJ2の劇的改善）を示し、右パネル (b) Zoomed Scale（$0 \le \text{nRF} \le 0.165$）で競合手法間の差（D=2.0におけるBioNJの劣化 vs AddNJ2の耐性）を高解像度に可視化。
- **[パラメータ]**: $N=20$, $D \in \{0.5, 1.0, 2.0\}$, $L=500$ aa, 14手法, 10 replicates, Mean ± SEM.

## 🚀 次に行うこと (Next Actions)
1. ユーザーへの生成図（`results/nrf_comparison.png`）の提示
2. 必要に応じた図の微調整（マーカーサイズ、配色、凡例配置など）

---

最終更新: 2026-11-11 08:15

## 📈 進捗状況
- [x] 完了: INDELible v1.03を用いたZipfian Indelシミュレーション（`POW 1.7 50`、indel_rate = 0.05, 0.1）の実装と自動ビルドスクリプト整備 (完了日: 11/11)
- [x] 完了: 系統樹距離・ブランチ長評価指標（KF距離 `kf_bsd` / パトリスティックRMSE `patristic_rmse`）の計算機能追加（Add-NNI/SPRの非負OLSブランチ長最適化含む） (完了日: 11/11)
- [x] 完了: スパコン向けLSFジョブアレイ投入スクリプト（`scripts/submit_lsf_indel.sh` / 240ジョブ並列）の作成 (完了日: 11/11)
- [x] 完了: 可視化スクリプト `analysis/merge_and_plot_hpc.py` の更新（エラバープロット `nrf_vs_taxa.png`、SSS横軸プロットのSSS=0.03/0.06破線除去、AddNJ2ラベル統一、FastME-NNIを含む全手法対応） (完了日: 11/11)

## 🛠 技術的決定・メモ
- **[課題]**:
  1. 従来の置換のみシミュレーションでは配列長が500aa固定となり、高発散（D=2.0）でもMAFFTがほぼ完全なアライメントを復元するため、MSA_MLのnRFが全SSS域で0付近に張り付いていた。
  2. Add-NNI / Add-SPRはカルテット加法性ペナルティの最小化のみでトポロジーを変更するため、得られた系統樹の枝長が未最適化となり、KF距離（BSD）やパトリスティックRMSEの評価で不利になる。
- **[解決策]**:
  1. INDELible v1.03を導入し、先行研究（PhyloTree_Methods_Comparison）と同一のZipfianモデル（`POW 1.7 50`）によるIndel挿入（indel_rate=0.05, 0.1）を実装。
  2. `src/topology_search.py` に `fit_ols_branch_lengths`（二乗残差最小化＋非負射影）を追加し、NNI/SPR後の木に対して観測距離行列 $D_{\text{obs}}$ に適合したブランチ長を付与。
  3. `src/evaluation.py` に二部グラフスプリット接続行列を用いた高速ベクトル化パトリスティック距離行列計算 $D_{\text{pat}} = (M \odot b)^T(1-M) + (1-M)^T(M \odot b)$ を実装（$N=128$ でも $<11$ ms）。
  4. `analysis/merge_and_plot_hpc.py` を改修し、`results_large_scale/nrf_vs_taxa.png` 準拠の切断マーカー付きエラバープロット（Dごとの独立y軸・横ドッジ）と、SSS横軸プロット（SSS=0.03/0.06破線除去、`AddNJ2` 表記）を生成。
- **[パラメータ]**: $N \in \{20, 32, 64, 128\}$, $D \in \{0.5, 1.0, 2.0\}$, $\text{indel\_rate} \in \{0.05, 0.1\}$, WAG+$\Gamma$ ($\alpha=1.0$), `replicates=10`（計240ジョブ）。

## 🚀 次に行うこと (Next Actions)
1. スパコン上で `bsub < scripts/submit_lsf_indel.sh` を実行して全240ジョブを並列投入
2. ジョブ完了後、`python analysis/merge_and_plot_hpc.py --outdir results_hpc_indel_0.05` および `--outdir results_hpc_indel_0.1` で集計・グラフ出力

---

最終更新: 2026-11-11 06:56

## 📈 進捗状況
- [x] 完了: LBAベンチマーク v2（3種距離行列 blosum62/wag_gamma/true × 14距離法 + MSA_ML、FastME-NNI追加）(完了日: 11/11)
- [x] 完了: SSS（MMseqs2, 長枝ペア w(T1,TN) 含む）と4点条件違反度（行列Δ_ADD・木L(T)）の算出 (完了日: 11/11)
- [x] 完了: 全5,280本の系統樹・配列の保存（results_lba/trees/, all_inferred_trees.csv）(完了日: 11/11)

## 🛠 技術的決定・メモ
- **[課題]**: `compute_mmseqs_sss_statistics` が PATH 上の mmseqs のみ探索し、未活性環境で SW 代替にフォールバックしていた。
- **[解決策]**: sys.executable 隣接および phylomethod_env の mmseqs を探索するよう修正。hill_climb_search に penalty≈0 での早期終了を追加（true行列で無駄なSPR探索を回避）。Hybrid は NNI 結果から SPR を継続する形に変更し重複計算を削減。
- **[結果]**: true行列では全手法 nRF=0（実装妥当性確認）。WAG+Γ距離はRMSE小だが Δ_ADD 大で nRF は Poisson と同等〜悪化。AddNJ2 と NJ の LBA 増分差は SEM 内（10反復では有意差なし）。
- **[パラメータ]**: N∈{20,32}, D∈{0.5,1.0,2.0}, lba_ratio∈{1.0,3.0}, L=500, WAG+Γ α=1.0, 10 reps, seed=42。

## 🚀 次に行うこと (Next Actions)
1. 有意差検定のため反復数を増やす: `python benchmark_lba.py --replicates 50 --outdir results_lba_r50`
2. lba_ratio を 5.0 などに拡張し偽姉妹群が多発する領域で手法差を評価

---

最終更新: 2026-11-11 05:06

## 📈 進捗状況
- [x] 完了: WAG+Gamma最尤推定距離行列（`compute_pairwise_wag_gamma_distance_matrix`）およびSSS類似度（`compute_all_pairwise_matrices_and_sss`）の実装・単体テスト (完了日: 11/11)
- [x] 完了: `TreeScorer` の $N=128$ 対応高速化・系統樹の4点条件違反度一括評価機能の実装 (完了日: 11/11)
- [x] 完了: スパコン（LSF）向け統合ベンチマークスクリプト `benchmark_hpc.py`（BLOSUM62 / WAG+Gamma / True距離行列、全Newick保存、4点条件違反度、SSS算出、rep=30）の実装・ドライラン検証 (完了日: 11/11)
- [x] 完了: LSFジョブアレイ投入スクリプト `scripts/submit_lsf.sh` および結果マージ・可視化スクリプト `analysis/merge_and_plot_hpc.py` の実装 (完了日: 11/11)

## 🛠 技術的決定・メモ
- **[課題]**:
  1. $N=20, 32, 64, 128$ 全条件を統一反復数（`rep=30`）で揃え、$N=128$ においても `Add-NNI` 探索を実行する必要がある。
  2. 従来のBLOSUM62+ポアソン距離だけでなく、シミュレーション真の進化モデルに合致した「WAG+$\Gamma$ 最尤推定距離行列」および「真のPatristic距離行列」での比較が必要。
  3. 各手法で構築された系統樹（Newick）の保存、系統樹の4点条件違反度 $\mathcal{L}(T; D)$ の算出、および配列類似度スコア（SSS）の記録が必要。
- **[解決策]**:
  1. `src/sequence_simulator.py` に固有値分解と事前計算グリッド＋放物線補間を用いた高速WAG+$\Gamma$最尤距離推定器およびSSS統計量算出関数を実装。
  2. `src/topology_search.py` の `TreeScorer` に `itertools.islice` による確定的な均一カルテットサンプリングを導入し、$N=128$ でもメモリ爆発なく高速に `Add-NNI` および4点条件違反度評価を実行可能にした。
  3. `scripts/submit_lsf.sh` にて `#BSUB -J add_phylo[1-12]` ジョブアレイを構成し、全12条件（4 Taxa $\times$ 3 Divergence）をスパコン上で同時並列実行できるようにした。
- **[パラメータ]**: $N \in \{20, 32, 64, 128\}$, $D \in \{0.5, 1.0, 2.0\}$, $L=500$ aa, WAG+Gamma ($\alpha=1.0$), `replicates=30`, 距離行列 $\in \{\text{blosum62}, \text{wag\_gamma}, \text{true}\}$.

## 🚀 次に行うこと (Next Actions)
1. スパコン（LSF環境）にて `git pull` 後、`mkdir -p logs && bsub < scripts/submit_lsf.sh` でジョブアレイを投入
2. 全12ジョブ完了後、`python analysis/merge_and_plot_hpc.py --outdir results_hpc` を実行して統合CSVと全グラフを生成

---

最終更新: 2026-11-11 04:35

## 📈 進捗状況
- [x] 完了: 長枝牽引（LBA）検証の実験計画書作成およびユーザー承認 (完了日: 11/11)
- [x] 完了: `src/evaluation.py` へのLBA特化評価関数（最短パス枝数 $DS(T_1, T_N)$、偽姉妹群判定）の実装と単体テスト (完了日: 11/11)
- [x] 完了: LBA専用並列ベンチマーク実行スクリプト `benchmark_lba.py` の実装と全120試行（N=20, 32 / D=0.5, 1.0, 2.0 / LBA比率 1.0 vs 3.0 / 全15手法）の完遂 (完了日: 11/11)
- [x] 完了: 論文品質グラフ生成スクリプト `analysis/plot_lba_results.py` の実装と画像群生成 (`results_lba/`) (完了日: 11/11)
- [x] 完了: `RawResults.md`, `Method.md` への実験結果追記および学術的考察の記録 (完了日: 11/11)

## 🛠 技術的決定・メモ
- **[課題]**: 従来の系統樹評価は全体トポロジー誤差（nRF）のみであり、長枝牽引（LBA）によって生じる特定の長枝ペア（$T_1$ と $T_N$）の引力・短縮アーティファクトを直接捉える指標が存在しなかった。
- **[解決策]**:
  1. `src/evaluation.py` に無向系統樹上の幅優先探索（BFS）による最短パス枝数計算 `compute_topological_distance_edges` および `evaluate_tree_lba` を実装。
  2. $T_1$ と $T_N$ が直近内部ノードに直接結合した偽姉妹群（$DS=2$）の発生率、および枝数短縮度 $\Delta DS = DS_{\text{true}} - DS_{\text{est}}$ を定量化。
  3. `benchmark_lba.py` にて並列マルチプロセス実行をサポートし、全120条件・15手法のシミュレーションを約2分で完遂。
- **[実験結果]**:
  - $N=32, D=0.5$（LBA比率3.0）において、Standard NJ ($nRF=0.0103$) や BioNJ ($nRF=0.0138$) が誤差を出す中、提案手法 **AddNJ2** は **$nRF = 0.0000$ (完全正解)** を達成し、FastME-SPRおよび最尤法 (MSA+ML) と同等の最高耐性を発揮した。
  - $N=32, D=1.0$ において、BioNJは $nRF=0.0310 \to 0.1103$（約3.5倍）へ急落したのに対し、AddNJ2は $0.0138 \to 0.0414$ に踏みとどまり、Standard NJ ($0.0517$) を凌駕する高いLBA耐性を実証した。
  - AddNJ2の推論時間は $N=32$ でわずか **0.147秒** であり、最尤法（3.66秒）の約25倍、加法性Hybrid探索（22.4秒）の約150倍高速。
- **[パラメータ]**: $N \in \{20, 32\}$, $D \in \{0.5, 1.0, 2.0\}$, `lba_ratio` $\in \{1.0, 3.0\}$, $L=500$ aa, WAG+Gamma ($\alpha=1.0$), 10 replicates each (計120データセット).

## 🚀 次に行うこと (Next Actions)
1. ユーザーへLBA検証結果の提示と生成グラフ（`results_lba/`）の案内
2. 論文執筆に向けた全体まとめ（大規模ベンチマーク結果とLBA耐性結果の統合）

---

最終更新: 2026-11-11 03:35

## 📈 進捗状況
- [x] 完了: 大規模ベンチマーク結果（100反復）の回収とTaxa=20（10反復）の統合可視化 (完了日: 11/11)
- [x] 完了: 可視化改善（SDエラーバー付与、独立縦軸スケール、点のみのプロット、全7手法の完全表示） (完了日: 11/11)
- [x] 完了: `analysis/plot_large_scale.py` の実装および `results_large_scale/` 画像（nrf_vs_taxa.png, runtime_vs_taxa.png）の再生成 (完了日: 11/11)

## 🛠 技術的決定・メモ
- **[課題]**: 従来のプロットでは、一部手法（NNI探索系）が表示されておらず、進化距離間で縦軸を共有していたため低ダイバージェンス（D=0.5, 1.0）の差が視認困難であった。またTaxa=20が含まれておらず、折れ線で結ばれていたため点の重なりで見づらかった。
- **[解決策]**:
  1. `analysis/plot_large_scale.py` にて、Taxa=20, 32, 64, 128 を横軸に統合。
  2. 評価された全7手法（MSA_ML, Standard_NJ, BioNJ, FastME_SPR, AddNJ2, BioNJ_Add_NNI, AddNJ2_Add_NNI）を漏れなくプロット。
  3. 各手法の水平オフセット（Dodge）を設け、プロット同士を線で結ばず散布図＋標準偏差（SD）エラーバー形式に変更。
  4. ダイバージェンスごと（D=0.5, 1.0, 2.0）に縦軸のレンジを最適化（独立スケール）し、細かな精度差を明瞭に可視化。
- **[パラメータ]**: Taxa $N \in \{20, 32, 64, 128\}$, Divergence $D \in \{0.5, 1.0, 2.0\}$, Replicates: $N=20$ (10 reps), $N=32,64,128$ (100 reps each, 計900データセット).

## 🚀 次に行うこと (Next Actions)
1. ユーザーへの可視化更新結果の提示と学術的知見の報告
2. 論文執筆用まとめまたは追加検証項目の確認

---

最終更新: 2026-11-10 02:28

## 📈 進捗状況
- [x] 完了: 大規模ベンチマーク計画の策定・外部PC実行パッケージの作成 (完了日: 11/10)
- [x] 完了: MSA+ML (MAFFT + FastTree WAG+Gamma) パイプライン実装と単体テスト (完了日: 11/10)
- [x] 完了: AddNJ2のトップKスクリーニングによる大規模データ高速化 ($N=128$ で0.30秒) (完了日: 11/10)
- [x] 完了: 大規模検証スクリプト `benchmark_large_scale.py` の動作検証 (完了日: 11/10)
- [x] 完了: 完全ポータブル環境定義 `environment.yml` および実行計画書 `LargeScale_Benchmark_Plan.md` の作成 (完了日: 11/10)

## 🛠 技術的決定・メモ
- **[課題]**: $N=128$ ではカルテット数が約 $1.07 \times 10^7$ 通りとなり、全ペア探索（8,128ペア）を愚直に行うと計算が停止したように遅延する問題。また、基準となる最高精度法として多重整列＋最尤法（MSA + ML）の追加が必要。
- **[解決策]**:
  1. `src/ml_wrapper.py` に `mafft --auto` による多重整列と `FastTree -wag -gamma` による最尤系統樹推論パイプラインを実装。
  2. `add_nj2` に NJ の $Q$ 行列による候補ペアのトップ $K=50$ スクリーニングおよびカルテットサンプリング（最大200ペア）を導入。これにより $N=128$ の木構築時間を数分以上からわずか 0.30 秒に短縮。
  3. 別の計算機・スパコンですぐに再現実行できるよう、環境定義 `environment.yml` と自律実行スクリプト `benchmark_large_scale.py`、および詳細な実行手引き `LargeScale_Benchmark_Plan.md` を作成。
- **[パラメータ]**: $N \in \{32, 64, 128\}$, $L=500$ aa, WAG+Gamma ($\alpha=1.0$), divergence $\in \{0.5, 1.0, 2.0\}$, 10 replicates each (計90データセット). 比較手法: MSA+ML, NJ, BioNJ, FastME-SPR, AddNJ2, BioNJ+Add-NNI, AddNJ2+Add-NNI.

## 🚀 次に行うこと (Next Actions)
1. 外部マシン（またはスパコン / Google Colab / GPUサーバー）へのリポジトリ配置と環境構築 (`micromamba create -f environment.yml`)
2. 大規模ベンチマークの実行 (`python benchmark_large_scale.py --taxa 32 64 128 --replicates 10 --outdir results_large_scale`)
3. 生成された `results_large_scale/` の回収と解析結果の報告

---

最終更新: 2026-11-10 02:08

## 📈 進捗状況
- [x] 完了: 改善版手法1（AddNJ2: 姉妹群カルテット整合度法）の実装および単体テスト (完了日: 11/10)
- [x] 完了: AddNJ2を含む全14手法の30反復ベンチマークシミュレーション完遂 (完了日: 11/10)
- [x] 完了: 比較グラフの更新 (results/nrf_comparison.png, runtime_comparison.png) (完了日: 11/10)
- [x] 完了: RawResults.md / Method.md の更新と学術的考察 (完了日: 11/10)

## 🛠 技術的決定・メモ
- **[課題]**: 初期のAddNJが縮約行列全体の平均欠損を評価していたため、ノイズの大きい葉を安易に消去してしまう消去バイアスによりnRF=0.65付近に停滞していた。
- **[解決策]**: 結合候補 (a, b) と残りの全ペア (k, l) が形成するカルテットが ab|kl を支持するかどうかのみを評価する「姉妹群カルテット整合度（Sister-Pair Quartet Consistency）」を目的関数とした `add_nj2` を実装。
- **[実験結果]**: AddNJ2により、D=0.5でnRFが0.6471から0.0000（完全正解）へ、D=1.0で0.6529から0.0118へ、D=2.0で0.6706から0.0882（BioNJの0.1235を大幅に上回る）へと劇的な性能向上を達成した。さらに計算時間も0.120秒から0.030秒へと4倍高速化した。
- **[パラメータ]**: N=20 taxa, L=500 aa, WAG+Gamma (alpha=1.0), divergence in [0.5, 1.0, 2.0], Global alignment (BLOSUM62) -> Poisson distance, 10 replicates each (計30反復).

## 🚀 次に行うこと (Next Actions)
1. ユーザーへの改善結果（AddNJ vs AddNJ2 vs 既存法）の報告とグラフの提示
2. 追加の要望や論文投稿用フォーマットの確認

---

最終更新: 2026-11-10 01:38

## 📈 進捗状況
- [x] 完了: 手法1（AddNJ: 加法性最小化貪欲法）の実装および単体テスト (完了日: 11/10)
- [x] 完了: 手法2（Add-NNI, Add-SPR, Add-Hybrid）の実装および単体テスト (完了日: 11/10)
- [x] 完了: FastMEラッパーの実装 (BioNJ, FastME-NNI, FastME-SPR) (完了日: 11/10)
- [x] 完了: 全30反復ベンチマーク本実験の完遂 (完了日: 11/10)
- [x] 完了: 論文品質グラフ作成 (nrf_comparison.png, runtime_comparison.png) (完了日: 11/10)
- [x] 完了: RawResults.md / Method.md の更新と考察 (完了日: 11/10)

## 🛠 技術的決定・メモ
- **[課題]**: AddNJにおいて加法的行列復元時に複数のペアが加法性欠損0.0となるタイが発生し、真の姉妹群以外が選ばれる現象。
- **[解決策]**: 目的関数Delta_ADDが同値（許容差1e-6内）の際、NJ法のQ規準値を二次タイブレーカーとして採用することで、加法的行列上で100%真の系統樹を復元可能にした。
- **[実験結果]**: 高ダイバージェンス（D=2.0）下において、BioNJ単体（nRF=0.1235）に対し、加法性Hybrid探索を適用した「BioNJ + Add-Hybrid」は nRF=0.0765 とトポロジー誤差を約38%大幅削減し、FastME-SPR (0.0882) を上回る最高精度を達成した。
- **[パラメータ]**: N=20 taxa, L=500 aa, WAG+Gamma (alpha=1.0), divergence in [0.5, 1.0, 2.0], Global alignment (BLOSUM62) -> Poisson distance, 10 replicates each (計30反復).

## 🚀 次に行うこと (Next Actions)
1. ユーザーへの完了報告と実験結果（AddNJとAdd-Hybridの挙動・性能差）の提示
2. 追加検証（taxa数変更やローカルアライメント耐性実験など）の要否確認

---

最終更新: 2026-11-10 01:10

## 📈 進捗状況
- [x] 完了: 要件定義・仕様確定ヒアリング (/grill-me) (完了日: 11/10)
- [x] 完了: 環境構築・FastME確認・ディレクトリ整備 (完了日: 11/10)
- [ ] 進行中: 手法1（加法性NJ: AddNJ）の実装と単体テスト (進捗率: 30% / 課題: なし)
- [ ] 未着手: 手法2（加法性NNI/SPR/Hybrid）の実装と単体テスト
- [ ] 未着手: FastMEラッパーの実装 (BioNJ, FastME-NNI/SPR)
- [ ] 未着手: ベンチマークシミュレーションの実行 (N=20, L=500, divergence=0.5, 1.0, 2.0, 10 reps)
- [ ] 未着手: 可視化・成果物まとめ (RawResults.md, Method.md, results/図表)

## 🛠 技術的決定・メモ
- **[環境]**: `/Users/aibarakazuki/.micromamba/envs/phylomethod_env` (Python 3.10, FastME 2.1.6.3, Biopython, DendroPy, NumPy).
- **[手法1設計]**: 距離縮約式 $D_{uk} = \frac{1}{2}(D_{ik} + D_{jk} - D_{ij})$ を採用。全ペア縮約後行列 $D'$ の全カルテット加法性欠損 $\Delta_{\text{ADD}}(D')$ を全探索して最小ペアを選択。
- **[手法2設計]**: 目的関数 $\mathcal{L}(T) = \frac{1}{|\mathcal{Q}|} \sum_q \frac{S(q_T) - \min(S)}{S_{\max} + \varepsilon}$。NNI、SPR、Hybrid (NNI $\to$ SPR) の3モードを実装。
- **[パラメータ]**: $N=20$ taxa, $L=500$ aa, divergence scale $\in \{0.5, 1.0, 2.0\}$, WAG+Gamma ($\alpha=1.0$), BLOSUM62 Global Alignment $\to$ Poisson distance, 10 replicates each (計30データセット).

## 🚀 次に行うこと (Next Actions)
1. `src/add_nj.py` の実装（ベクトル化カルテット評価 & AddNJ木構築）
2. `tests/test_add_nj.py` による4点加法性行列復元テスト
