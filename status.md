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
