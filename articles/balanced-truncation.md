---
title: "Balanced truncation法による複素指数関数和のモデル縮約"
emoji: "🍡"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: ["アルゴリズム","julia"]
published: false
---

## 概要

モデル縮約（Model reduction）はシステム制御理論における。


## アルゴリズムの導入

すでに複素指数関数和

$$
f(t) = \sum_{k=1}^M c_k \mathrm{e}^{-a_k t}, \quad a_k,c_k \in \mathbb{C}
$$

が与えられているとする。Balanced truncation法の目的は、与えられた精度 $\epsilon$ に対して、より少ない項数を持つ別の最適な指数和を用いて $f(t)$ を近似することであり、すなわち

$$
C_{\mathrm{BTM}}(t)=\sum_{k=1}^{M'} c_k' \mathrm{e}^{-a_k' t}, 
\quad
\left\|\tilde{C}(t)-C_{\mathrm{BTM}}(t)\right\| < \epsilon 
$$

が成り立つような $C_{\mathrm{BTM}}(t)$ を考える。

$C(t)$ をラプラス変換すると、

$$
\tilde{C}(s) 
= \mathcal{L}\bigl[\tilde{C}(t)\bigr] 
= \int_{0}^\infty \mathrm{d}t\, \tilde{C}(t) \, \mathrm{e}^{-st} 
= \sum_{k=1}^{M} \frac{c_k}{s+a_k}
$$

となる。すると問題は、高次の有理関数を低次の有理関数で近似する最適化問題に変形され、その形は

$$
C_\mathrm{BTM}(s)
= \sum_{k=1}^{M} \frac{c_k'}{s+a_k'}, 
\quad
\left\|\tilde{C}(s)-C_{\mathrm{BTM}}(s)\right\| < \epsilon'
$$

となる。これはすべての $s > 0$ について成立する。  
この問題は線形時不変 (LTI) システムのモデル縮約として理解でき、$C(s)$ は LTI システムの伝達関数に対応する（らしい）。BTM の利点は、近似誤差に対して $L^\infty$ の上界を与える点にある。

次に$c_k,a_k$ を用いて、可制御性グラミアン $\mathbf{W}_\mathrm{c} \in \mathbb{C}^{M\times M}$ を構成する。これは

$$
[\mathbf{W}_\mathrm{c}]_{ij} = \frac{\sqrt{c_i c_j^*}}{a_i + a_j^*}
$$

で与えられる。このようにして得られる $\mathbf{W}_\mathrm{c}$ は自己随伴な準コーシー行列であり、次のように分解できる。

$$
\mathbf{W}_\mathrm{c} 
= \bar{\mathbf{U}} \,\boldsymbol{\Sigma}\, \mathbf{U}^{T}
$$

ここで、$\mathbf{U} \in \mathbb{C}^{M\times M}$ はユニタリ行列、$\boldsymbol{\Sigma} = \mathrm{diag}(\sigma_1,\dots,\sigma_M)$ は $\sigma_1 \ge \sigma_2 \ge \cdots \ge \sigma_M > 0$ をみたす対角行列となる。 これはcon-eigenvalue decompositionと呼ばれる（）。

整数 $M'$ が 

$$
2\sum_{i=M'+1}^{M} \sigma_i \leq \epsilon
$$

をみたすとき、$\mathbf{A}' \in \mathbb{C}^{M'\times M'}$ を

$$
\mathbf{A}' = \mathbf{U}_{M'}^{*}\,\mathbf{A}\,\bar{\mathbf{U}}_{M'}
$$

および、$\mathbf{B}' \in \mathbb{C}^{M'}$ を

$$
\mathbf{B}' = \mathbf{U}_{M'}\,\mathbf{B}
$$

によって構成する。ここで、$\mathbf{A} = \operatorname{diag}(a_1, a_2, \ldots, a_M)$、$\mathbf{B} = (\sqrt{c_1}, \sqrt{c_2}, \ldots, \sqrt{c_M})^T$、 $\mathbf{U}_{M'} = \mathbf{U}(1\!:\!M,\,1\!:\!M')$ である。

最後に、$\mathbf{A}'$ の固有値分解

$$
\mathbf{A}' = \mathbf{X} \,\boldsymbol{\Lambda}\,\mathbf{X}^{-1}
$$

を計算すると、

$$
a_{i}'=\boldsymbol{\Lambda}_{ii}, 
\quad
c_i'=(\mathbf{B}_i^{''})^{2} 
\quad
(i=1, \ldots, M')
$$

が得られる。ここで、$\mathbf{B}^{''} = \mathbf{X}^{-1}\mathbf{B}'$ である。


## ソースコード

C++版は
https://github.com/hydeik/mxpfit
で実装されている。

---

[1^] 



