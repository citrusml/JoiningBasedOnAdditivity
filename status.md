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
