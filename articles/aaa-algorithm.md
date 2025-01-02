---
title: "AAAアルゴリズムによる"
emoji: 🍙
type: "tech" # tech: 技術記事 / idea: アイデア
topics: ["アルゴリズム"]
published: false
---

## 概要

AAAアルゴリズムは、与えられた精度あるいは与えられた近似次数に対して有理近似を求めるために用いられるアルゴリズムである。ここではアルゴリズムの主な特徴について簡単に概説し、数値的な詳細は原論文[^1][^2]を参照していただきたい。
とある用途を想定しているため、本記事では関数$F(\omega)$を有理関数

で近似することを考える。このような近似を行えると、例えばフーリエ変換等の積分が以下のように解析的に行えるため嬉しい。


## アルゴリズムの導入

区間 $[-\omega_\mathrm{min}, \omega_\mathrm{max}]$ を $N$ 個の点に離散化する。こうして

$$
\boldsymbol{\Omega} = (\omega_1, \dots, \omega_N)^\mathrm{T}, 
\quad
\mathbf{F} = (F(\omega_1), \dots, F(\omega_N))^\mathrm{T}
$$

を得る。

$$
\| F(\omega) - \tilde{F}(\omega)\| < \epsilon', \quad \omega \in \boldsymbol{\Omega}
$$

[^1]。ここで

$$
\tilde{F}(\omega) = \sum_{k=1}^{M'} \frac{\eta_k}{\omega - z_k},
\quad
z_k, \eta_k \in \mathbb{C}
$$

である。

## AAAアルゴリズムの手順

AAAアルゴリズムは反復型アルゴリズムの一種である[^1]。ステップ $m$ において、入力データをサポート点

$$
\boldsymbol{\Omega}^{(m)} = (\omega_1, \dots, \omega_m)^\mathrm{T}
$$

と、その補集合

$$
\overline{\boldsymbol{\Omega}^{(m)}} = \boldsymbol{\Omega} \backslash \boldsymbol{\Omega}^{(m)} 
\equiv (\omega^{(m)}_1, \dots, \omega^{(m)}_{N-m})^\mathrm{T}
$$

に分割する。そして、バリセントリック形式で定義される有理関数

$$
F^{(m-1)}(\omega)
=
\frac{\displaystyle \sum_{k=1}^m \frac{w_k F(\omega_k)}{\omega - \omega_k}}
{\displaystyle \sum_{k=1}^m \frac{w_k}{\omega - \omega_k}}
$$
（式4）

を考える。次に、重み $\mathbf{w}=(w_1,\dots,w_m)^\mathrm{T}$ に関して、

$$
\text{minimize}
\quad
\| F(\omega) - F^{(m-1)}(\omega)\|_{\boldsymbol{\Omega}}
\quad
\mathrm{s.t.}
\quad
\|\mathbf{w}\|_{m}=1
$$

という線形最小二乗問題を解く。ここで $\|\cdot\|_{\Omega}$ は $\boldsymbol{\Omega}$ 上のL2-ノルム、$\|\cdot\|_{m}$
は $m$-次元ベクトル上の離散2-ノルムを表す。簡単な操作により、この問題は行列形式

$$
\text{minimize}
\quad
\|\mathbf{A}^{(m)}\mathbf{w}\|_{M-m}
\quad
\mathrm{s.t.}
\quad
\|\mathbf{w}\|_m=1
$$

へと書き換えられる。ここで、ロイナー行列 \(\mathbf{A}^{(m)}\) は

$$
\mathbf{A}^{(m)} =
\begin{pmatrix}
\dfrac{F(\omega^{(m)}_1)-F(\omega_1)}{\omega^{(m)}_1-\omega_1} & \cdots & \dfrac{F(\omega^{(m)}_1)-F(\omega_m)}{\omega^{(m)}_1-\omega_m}\\[6pt]
\vdots & \ddots & \vdots \\[6pt]
\dfrac{F(\omega^{(m)}_{M-m})-F(\omega_1)}{\omega^{(m)}_{M-m}-\omega_1} & \cdots & \dfrac{F(\omega^{(m)}_{M-m})-F(\omega_m)}{\omega^{(m)}_{M-m}-\omega_m}
\end{pmatrix}
$$
（式5）

と定義される。次に、この行列 \(\mathbf{A}^{(m)}\) に対して特異値分解 (SVD)

$$
\mathbf{A}^{(m)} = \mathbf{U} \boldsymbol{\Sigma} \mathbf{V}^*
$$

を行い、最後の右特異ベクトル \(\mathbf{V}\) を用いて

$$
\mathbf{w} = \mathbf{V}(1:m, m)
$$

と重み \(\mathbf{w}\) を求める。

ステップ \(m\) において、もし

$$
\| F(\omega) - F^{(m-1)}(\omega)\|_{\boldsymbol{\Omega}} < \epsilon'
$$

が満たされたならば反復を終了し、そのときの次数は \(M' = m-1\) となる。満たされない場合は、貪欲法[^1]によって補集合 \(\overline{\boldsymbol{\Omega}^{(m)}}\) から次のサポート点 \(\omega_{m+1}\) を選ぶ。

さらに、極 \(z_k\;(k=1,\dots,m-1)\) は、アローヘッド型の \((m+1)\times(m+1)\) 一般化固有値問題

$$
\begin{pmatrix}
0 & w_1 & w_2 & \cdots & w_m \\
1 & \omega_1 & & & \\
1 & & \omega_2 & & \\
\vdots & & & \ddots & \\
1 & & & & \omega_m
\end{pmatrix}
=
z_k
\begin{pmatrix}
0 & & & & \\
& 1 & & & \\
& & 1 & & \\
& & & \ddots & \\
& & & & 1
\end{pmatrix}
$$

の固有値として評価できる[^1]。また、留数 \(\eta_k \;(k=1,\dots,m-1)\) は、解析関数の商の単純極に対する公式を用いて

$$
\eta_k 
= \frac{\displaystyle \sum_{k'=1}^m \dfrac{w_{k'}\, F(\omega_{k'})}{z_k - \omega_{k'}}}
{\displaystyle \sum_{k'=1}^m \Bigl[-\dfrac{w_{k'}}{(z_k - \omega_{k'})^2}\Bigr]}
$$

によって計算される。

## 最終的な近似式

最後に基本定理を用いることで、\(M \sim \frac{M'}{2}\) 項をもつ最終的な式

$$
\tilde{C}(t)
= \frac{1}{2\pi} \int_{-\infty}^{\infty} \mathrm{d}\omega\, F(\omega)\,\mathrm{e}^{-i \omega t}
= -i \sum_{k=1}^{M} \eta_k\, \mathrm{e}^{-i t z_k^{\prime}}
= \sum_{k=1}^{M} c_k\, \mathrm{e}^{-a_k t}
$$
（式6）

を得る。ここで、\(z_k^{\prime}\) は複素平面の下半分における \(F(\omega)\) の極、\(c_k \equiv -i \eta_k\)、\(a_k \equiv i z_k^{\prime}\) である。

## コード

MATLAB版は

Python版は

で実装されている。Python版をJuliaに移植したので、そのコードを以下に示す。

---

[^1]: Nakatsukasa, Yuji, et al. *SIAM Journal on Scientific Computing* (2018).  
[^2]: Nakatsukasa, Yuji, et al. *SIAM Review* (2020).


