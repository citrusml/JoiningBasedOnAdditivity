# Method: 加法性に基づく系統樹推定アルゴリズムと実験手法

## 1. 理論的背景：4点条件（Bunemanの定理）と加法性欠損

距離行列 $D \in \mathbb{R}^{n \times n}$ が系統樹上の枝長最短経路距離（パトリスティック距離）と一致する必要十分条件は、任意の4つのタクソン（カルテット） $(i, j, k, l)$ に対し、ペア間距離の3つの和
$$
S = \{ D_{ij} + D_{kl}, \; D_{ik} + D_{jl}, \; D_{il} + D_{jk} \}
$$
を降順ソートした $A \ge B \ge C$ において、上位2つの和が等しい（$A = B$）ことです（4点条件、Four-Point Condition）。

1つのカルテットに対する相対加法性欠陥（Relative Additivity Defect） $\delta_{\text{ADD}}$ は以下のように定義されます：
$$
\delta_{\text{ADD}}(i, j, k, l) = \frac{A - B}{A + \varepsilon} \quad (\varepsilon = 10^{-7})
$$
全カルテット集合 $\mathcal{Q}$ に対する平均相対加法性欠損は以下の通りです：
$$
\Delta_{\text{ADD}}(D) = \frac{1}{|\mathcal{Q}|} \sum_{(i,j,k,l) \in \mathcal{Q}} \delta_{\text{ADD}}(i, j, k, l)
$$

---

## 2. 手法1: 加法性最小化貪欲結合法 (AddNJ)

従来のNeighbor-Joining (NJ) 法は、各ステップで推論される総枝長を最小化するペア $(i, j)$ を結合しますが、**AddNJ** では結合後の縮約行列 $D'$ における加法性欠損 $\Delta_{\text{ADD}}(D')$ を最小化するペアを貪欲に選択します。

### アルゴリズムの流れ
1. 現在のノード数を $m$ とする（初期値 $m = n$）。
2. $m > 4$ の間、以下を繰り返す：
   - 全ての候補ペア $(i, j)$ ($1 \le i < j \le m$) に対し、ペアを結合した仮想の親ノード $u$ を作成し、縮約行列 $D^{(i, j)} \in \mathbb{R}^{(m-1) \times (m-1)}$ を計算する。
     $$
     D_{uk} = \frac{1}{2} (D_{ik} + D_{jk} - D_{ij}) \quad (\forall k \ne i, j)
     $$
   - 縮約後行列の全カルテット平均加法性欠損 $\Delta_{\text{ADD}}(D^{(i, j)})$ を計算する。
   - $\Delta_{\text{ADD}}(D^{(i, j)})$ を最小化するペア $(i^*, j^*)$ を正式に結合対象として選択する。
   - 親ノード $u$ への枝長 $v_{i^*}, v_{j^*}$ を標準NJ法に従って算出する：
     $$
     r_x = \sum_{k=1}^m D_{xk}
     $$
     $$
     v_{i^*} = \frac{1}{2} D_{i^* j^*} + \frac{1}{2(m-2)} (r_{i^*} - r_{j^*}), \quad v_{j^*} = D_{i^* j^*} - v_{i^*}
     $$
   - 距離行列を $m-1$ 次元に更新し、$m \leftarrow m - 1$ とする。
3. $m = 4$ に到達した時：
   - 唯一のカルテット $(1, 2, 3, 4)$ における3つの和 $S_1 = D_{12}+D_{34}, S_2 = D_{13}+D_{24}, S_3 = D_{14}+D_{23}$ を計算する。
   - 最小の和を与えるペア（最も支持されるトポロジー）を結合する。
4. $m = 3$ に到達した時：
   - 残る3つのノードを中央ノードに結合し、標準NJの3点枝長計算を行って系統樹（Newick形式）を完成させる。

---

## 3. 手法2: 加法性に基づくトポロジー探索 (Add-Search: NNI / SPR / Hybrid)

初期系統樹（NJまたはBioNJ）を出発点として、加法性違反度を目的関数とする局所/広域探索を行います。

### 目的関数（Tree-induced Quartet Additivity Penalty）
系統樹 $T$ が与えられたとき、各カルテット $q = (i, j, k, l)$ は木上で唯一のトポロジー $q_T \in \{ij|kl, ik|jl, il|jk\}$ を誘導します。
観測距離行列 $D$ における対応するペア和を $S(q_T)$、3つのペア和の最小値を $S_{\min} = \min(S)$、最大値を $S_{\max} = \max(S)$ としたとき、カルテット加法性ペナルティを以下のように定義します：
$$
P(q \mid T) = \frac{S(q_T) - S_{\min}}{S_{\max} + \varepsilon}
$$
系統樹 $T$ の総合加法性ペナルティ $\mathcal{L}(T)$ は、全カルテットでの平均値です：
$$
\mathcal{L}(T) = \frac{1}{|\mathcal{Q}|} \sum_{q \in \mathcal{Q}} P(q \mid T)
$$
木 $T$ のトポロジーが全てのカルテットで4点条件の最善解と一致するとき、$\mathcal{L}(T) = 0$ となります。

### 探索アルゴリズム
1. **NNI (Nearest Neighbor Interchange)**:
   - 各内部枝に対して2つの代替トポロジーを生成（近傍サイズ $2(n-3)$）。
   - ペナルティ $\mathcal{L}(T)$ が最小となる近傍を採用する最急降下法を行い、局所最適解に達するまで反復。
2. **SPR (Subtree Pruning and Regrafting)**:
   - 内部枝を切り離し、他の枝へ再結合させる（広域探索）。
   - 最急降下法によりスコア改善がなくなるまで反復。
3. **Hybrid (NNI $\to$ SPR)**:
   - まずNNIで素早く局所最適化を行い、その後SPRで局所解からの脱出とさらなるスコア改善を図る多段階探索。

---

## 4. プロジェクトファイル構成と役割

- `src/add_nj.py`: 手法1（AddNJ）のコアアルゴリズム実装。ベクトル化されたカルテット計算と木構築。
- `src/topology_search.py`: 手法2（加法性NNI/SPR/Hybrid探索）の実装。DendroPyを用いた木構造操作。
- `src/fastme_wrapper.py`: 外部ツール FastME（BioNJ, FastME-NNI, FastME-SPR）の実行ラッパー。
- `src/sequence_simulator.py`: WAG+Gammaモデルによる配列進化シミュレータおよびBLOSUM62ポアソン距離計算。
- `src/tree_simulator.py`: 真の系統樹生成（Backward Yule過程）とパトリスティック距離行列の計算。
- `src/evaluation.py`: 評価指標（正規化Robinson-Foulds距離: nRF, カルテット加法性欠損 $\Delta_{\text{ADD}}$, 実行時間）。
- `main.py`: ベンチマーク実行用エントリポイント（30反復の自動実行、集計、可視化出力）。
- `tests/`: 各モジュールの pytest 単体テストコード群。
- `RawResults.md`: 全実験レプリケートの生データ記録。
- `status.md`: 開発および実験の進捗記録。

