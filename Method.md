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

### 手法1の課題（消去バイアス）
縮約行列全体の平均欠損 $\Delta_{\text{ADD}}(D')$ を評価すると、全体の約80%を占める $(i, j)$ と無関係な葉同士のカルテットの影響により、「真の姉妹群を選ぶ」ことよりも「ノイズの大きい葉を結合（除外）する」ことが優先される消去バイアスが生じることが判明した。

---

## 2.2 手法1改善: 姉妹群カルテット整合度法 (AddNJ2)

上記の消去バイアスを克服するため、**AddNJ2** では縮約行列全体ではなく、**結合候補 $(i, j)$ を真の姉妹群と仮定した際に成り立つべき4点条件の整合度**のみを目的関数として直接評価する。

### 数理モデル
2つのノード $(i, j)$ が系統樹上の真の姉妹群であるならば、残る任意のペア $(k, l)$ ($k, l \ne i, j$) に対し、カルテット $(i, j, k, l)$ は必ず $ij \mid kl$ のトポロジーを形成しなければならない。
したがって、ペア和 $S_1 = D_{ij} + D_{kl}, S_2 = D_{ik} + D_{jl}, S_3 = D_{il} + D_{jk}$ に対し、常に $S_1 \le \min(S_2, S_3)$ が成り立つ必要がある。
候補ペア $(i, j)$ の不整合ペナルティを以下のように定義する：
$$
\text{SisterDefect}(i, j) = \frac{1}{\binom{m-2}{2}} \sum_{k < l, \; k,l \notin \{i, j\}} \frac{(D_{ij} + D_{kl}) - \min(S)}{\max(S) + \varepsilon}
$$
各ステップで $\text{SisterDefect}(i, j)$ を最小化するペア $(i^*, j^*)$ を選択して結合する。
無関係なカルテットが評価から完全に除外されるため消去バイアスが発生せず、真の姉妹群が高精度に抽出される。


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
- `benchmark_lba.py`: 長枝牽引（LBA）効果検証の並列ベンチマーク実行スクリプト。
- `analysis/plot_lba_results.py`: LBA解析（DS推移、偽姉妹群率、nRF、実行時間）の論文品質可視化スクリプト。
- `RawResults.md`: 全実験レプリケートの生データ記録。
- `status.md`: 開発および実験の進捗記録。

---

## 5. 長枝牽引（LBA: Long Branch Attraction）のモデル化と評価手法

### 5.1 LBAアーティファクトの発生メカニズム
長枝牽引（Felsenstein 1978）は、真の系統樹上で遠く離れた2つの系統（Taxa $T_1$ および $T_N$）の末端枝長が極端に伸長した際に生じます。
ポアソン距離補正 $d = -\ln(1 - p)$ において、$p$ が飽和領域に近づくと分散が急増し、非姉妹群である長枝同士の収束変異やノイズによって、距離法（特に総枝長最小化を狙うNJ系手法）がこれらを誤って「姉妹群」として結合してしまう現象です。

### 5.2 系統樹シミュレーションにおけるLBA付与
BioPerl（`Bio::Tree::RandomFactory` による Backward Yule 過程）にて生成された樹形に対し、端点系統 $T_1$ および $T_N$ の末端枝長に倍率 `lba_ratio` を乗算します：
$$
l_{\text{scaled}}(T) = 
\begin{cases}
l(T) \times \text{lba\_ratio} & (T \in \{T_1, T_N\}) \\
l(T) & (\text{otherwise})
\end{cases}
$$
本研究では、対照群として `lba_ratio = 1.0`（LBAなし）、LBA負荷条件として `lba_ratio = 3.0`（先行研究 Matsui & Iwasaki 2020 準拠）を設定します。

### 5.3 LBA特化の評価指標

1. **最短パス枝数 $DS(T_1, T_N)$ (Shortest Path Distance)**:
   無向系統樹グラフにおいて、幅優先探索（BFS）により $T_1$ と $T_N$ を結ぶパスに含まれる枝数（エッジ数）を算出します。
   - 真の木における距離: $DS_{\text{true}}$
   - 推定木における距離: $DS_{\text{est}}$
   - 枝数短縮度: $\Delta DS = DS_{\text{true}} - DS_{\text{est}}$
2. **偽姉妹群形成率 (False Sister Rate, %)**:
   推定木において $DS_{\text{est}} = 2$（共通の直近内部ノードに直接結合した姉妹群）となってしまった試行の割合：
   $$
   \text{FalseSisterRate} = \frac{1}{R} \sum_{r=1}^R \mathbb{I}(DS_{\text{est}}^{(r)} = 2) \times 100\%
   $$
3. **加法性基準の耐性仮説**:
   加法性に基づく手法（AddNJ2）は、4点条件の不等式整合度（$S_1 \le \min(S_2, S_3)$）のみで結合ペアを選択するため、長枝の絶対値増大に惑わされにくく、従来のNJやBioNJよりも高いLBA耐性を発揮します。



### 5.4 LBA拡張ベンチマーク v2: 距離行列・SSS・4点条件違反度
- **距離行列3種**: `blosum62`（BLOSUM62 Global PSA + Poisson補正 $d=-\ln(1-p)$）、`wag_gamma`（PSAの20×20置換カウントから WAG+Γ(α=1.0) 下で最尤距離をグリッド探索＋放物線補間で推定）、`true`（真の木のパトリスティック距離）。
- **SSS**: MMseqs2 (`-s 7.5 -e 10`) の bit score $S$ から $w_{ij} = \frac{S(i,j)+S(j,i)}{S(i,i)+S(j,j)}$。MMseqs2 不在時は Smith-Waterman 正規化スコアで代替（`psa_sss_*` 列は常に併記）。`sss_lba_pair` は長枝ペア $w(T_1,T_N)$。
- **4点条件違反度**: 行列レベル $\Delta_{\text{ADD}}(D)$（`delta_add_*`）および木レベル $\mathcal{L}(T\mid D)$（`*_penalty_obs`, `*_penalty_true`）。
- **系統樹保存**: `results_lba/trees/N{N}_D{D}_LBA{r}/rep_XXX/{true_tree.nwk, sequences.fasta, MSA_ML.nwk, <matrix>/<method>.nwk}`。
- 実行: `python benchmark_lba.py --threads 6` → `python analysis/plot_lba_results.py`
