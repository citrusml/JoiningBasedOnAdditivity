# 大規模系統樹ベンチマーク計画書 (N=32, 64, 128 & MSA+ML)

本ドキュメントは、Taxa数（葉の数）を **$N = 32, 64, 128$** に拡張し、従来の距離法（NJ, BioNJ, FastME-SPR）や提案手法（AddNJ2, Add-NNI）に加え、**多重配列整列＋最尤法（MSA + ML: MAFFT + FastTree WAG+Gamma）** を統合した大規模ベンチマーク実験を**別のPCや計算クラスタ（スパコン）で即座に実行できるようにまとめた完全実行計画書**です。

---

## 1. 実験の目的と背景

1. **Taxa数増加に伴うスケーラビリティの検証**:
   - 小規模（$N=20$）で高い精度を示した `AddNJ2`（姉妹群カルテット整合度法）および加法性局所探索（`Add-NNI`）が、大規模配列群（$N=32, 64, 128$）においても既存の距離法（NJ, BioNJ, FastME）を凌駕できるかを検証する。
2. **最尤法（MSA + ML）との精度・速度トレードオフの比較**:
   - 進化学研究で最高精度とされる「MAFFT による多重整列 ＋ FastTree による WAG+$\Gamma$ 最尤推論」を基準ベンチマークに加え、距離行列ベースの手法（特に `AddNJ2`）がどこまで最尤法に迫れるか、また計算時間でどれほどの優位性があるかを明らかにする。

---

## 2. 比較対象手法一覧

| 手法表記 | アルゴリズム概要 | 使用ツール / スクリプト |
| :--- | :--- | :--- |
| **`MSA_ML`** | 多重配列整列 (MAFFT --auto) $\to$ 最尤系統樹推論 (FastTree -wag -gamma) | [`src/ml_wrapper.py`](file:///Users/aibarakazuki/JoiningBasedOnAdditivity/src/ml_wrapper.py) |
| **`Standard_NJ`** | 標準 Neighbor-Joining 法（BLOSUM62 ポアソン補正距離） | Biopython [`DistanceTreeConstructor`](file:///Users/aibarakazuki/JoiningBasedOnAdditivity/src/tree_inference.py) |
| **`BioNJ`** | 距離分散最小化重み付き NJ 法 (Gascuel 1997) | FastME `-m B` ([`src/fastme_wrapper.py`](file:///Users/aibarakazuki/JoiningBasedOnAdditivity/src/fastme_wrapper.py)) |
| **`FastME_SPR`** | Balanced Minimum Evolution (BME) 基準による SPR 探索 | FastME `-m B -s` ([`src/fastme_wrapper.py`](file:///Users/aibarakazuki/JoiningBasedOnAdditivity/src/fastme_wrapper.py)) |
| **`AddNJ2`** | **提案手法1改善**: 姉妹群カルテット整合度最小化貪欲結合法 | [`src/add_nj.py`](file:///Users/aibarakazuki/JoiningBasedOnAdditivity/src/add_nj.py) (`add_nj2`) |
| **`BioNJ_Add_NNI`** | **提案手法2**: BioNJ初期木 ＋ 加法性ペナルティ最小化 NNI 探索 | [`src/topology_search.py`](file:///Users/aibarakazuki/JoiningBasedOnAdditivity/src/topology_search.py) |
| **`AddNJ2_Add_NNI`**| **提案手法ハイブリッド**: AddNJ2初期木 ＋ 加法性ペナルティ NNI 探索 | [`src/topology_search.py`](file:///Users/aibarakazuki/JoiningBasedOnAdditivity/src/topology_search.py) |

---

## 3. 実験条件とパラメータ設計

- **Taxa数 ($N$)**: **32, 64, 128**（3水準）
- **進化モデル**: WAG + $\Gamma$ (4カテゴリ離散ガンマ分布, $\alpha=1.0$)
- **配列長 ($L$)**: 500 aa
- **Divergence Scale ($D$)**: **0.5, 1.0, 2.0** substitutions/site（3水準: 近縁〜高変異ノイズ）
- **反復回数**: 各条件 10 反復（合計 $3 \times 3 \times 10 = 90$ データセット）
- **評価指標**:
  1. **トポロジー精度**: 真の木との正規化 Robinson-Foulds 距離（nRF $\in [0, 1]$）
  2. **計算効率**: 1試行あたりの実行時間（秒）

---

## 4. 別のパソコン / スパコンでの環境構築手順

本リポジトリには依存パッケージを定義した [`environment.yml`](file:///Users/aibarakazuki/JoiningBasedOnAdditivity/environment.yml) が同梱されています。

### 方法 A: Conda / Micromamba による構築（推奨）

```bash
# 1. リポジトリをクローンまたはコピー
git clone <URL> JoiningBasedOnAdditivity
cd JoiningBasedOnAdditivity

# 2. environment.yml から一発で環境構築
micromamba create -f environment.yml -y
# （conda の場合: conda env create -f environment.yml）

# 3. 環境を有効化
micromamba activate phylomethod_env
# （conda の場合: conda activate phylomethod_env）

# 4. 単体テストの正常動作を確認
pytest tests/
```

### 方法 B: 外部ツールの個別インストール（Linux / Mac）

Conda を使わない場合、以下の 2 つの外部バイナリを PATH に通してください：
1. **FastME** (2.1.6.3): `conda install -c bioconda fastme` または公式サイトからコンパイル
2. **FastTree** (2.2.0): `conda install -c bioconda fasttree`
3. **MAFFT** (7.526): Ubuntu/Debian なら `sudo apt-get install mafft`、Mac なら `brew install mafft`

---

## 5. ベンチマーク実行コマンド

大規模ベンチマーク専用スクリプト [`benchmark_large_scale.py`](file:///Users/aibarakazuki/JoiningBasedOnAdditivity/benchmark_large_scale.py) を実行します。

### ① クイック動作確認（1反復のみ）
```bash
python benchmark_large_scale.py --taxa 32 64 --divergences 0.5 --replicates 1 --outdir test_results
```

### ② 本番実行（$N=32, 64, 128$, 全条件 10 反復）
```bash
python benchmark_large_scale.py \
  --taxa 32 64 128 \
  --divergences 0.5 1.0 2.0 \
  --replicates 10 \
  --outdir results_large_scale
```

> **[TIP] バックグラウンド実行（SSH 切断対策）:**
> 長時間実行となるため、リモートサーバーでは `nohup` や `tmux` の使用を推奨します：
> ```bash
> nohup python benchmark_large_scale.py --taxa 32 64 128 --replicates 10 > benchmark_large_scale.log 2>&1 &
> ```

---

## 6. 自動生成される出力成果物

実行完了後、`results_large_scale/` フォルダに以下の解析成果物が自動生成されます：

1. **`results_large_scale/large_scale_summary.csv`**:
   各レプリケートの全手法の nRF、RF、実行時間が記録された生データ CSV。
2. **`results_large_scale/nrf_vs_taxa.png` (300 DPI)**:
   Divergence（0.5, 1.0, 2.0）ごとに、Taxa 数（32, 64, 128）に対する各手法のトポロジー誤差（nRF 平均 $\pm$ SEM）の推移をプロットした比較図。
3. **`results_large_scale/runtime_vs_taxa.png` (300 DPI)**:
   Taxa 数 $N$ に対する各手法の計算時間（秒）を両対数スケール（log-log）でプロットした計算効率スケーリング図。

---

## 7. 実装ファイル一覧

- [`benchmark_large_scale.py`](file:///Users/aibarakazuki/JoiningBasedOnAdditivity/benchmark_large_scale.py): 大規模ベンチマーク実行エントリポイント
- [`src/ml_wrapper.py`](file:///Users/aibarakazuki/JoiningBasedOnAdditivity/src/ml_wrapper.py): MAFFT + FastTree 連携パイプライン
- [`src/add_nj.py`](file:///Users/aibarakazuki/JoiningBasedOnAdditivity/src/add_nj.py): $N=128$ に耐えうる高速スクリーニング付き `add_nj2`
- [`src/fastme_wrapper.py`](file:///Users/aibarakazuki/JoiningBasedOnAdditivity/src/fastme_wrapper.py): FastME (BioNJ, FastME-SPR) ラッパー
- [`src/topology_search.py`](file:///Users/aibarakazuki/JoiningBasedOnAdditivity/src/topology_search.py): 4点条件ペナルティ評価・NNI局所探索
- [`environment.yml`](file:///Users/aibarakazuki/JoiningBasedOnAdditivity/environment.yml): 動作保証環境定義ファイル

