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
