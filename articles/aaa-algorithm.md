---
title: "AAAアルゴリズムによる有理関数近似"
emoji: 🍙
type: "tech" # tech: 技術記事 / idea: アイデア
topics: ["アルゴリズム","julia"]
published: true
---

## 概要

AAA (adaptive Antoulas–Anderson) アルゴリズムは、与えられた精度あるいは与えられた近似次数に対して有理関数近似を求めるアルゴリズムである。[^1][^2]このアルゴリズムは、Antoulas–Anderson[^3]に倣って有理関数をbarycentric形式で表し、サポート点を段階的に追加して近似次数を高める適応的な戦略を採用することで、Padé近似等の有理近似で現れる偽の極ー零点対をほとんど回避することができるのが特徴である。
本稿では、AAAアルゴリズムを概観する。詳細は原論文[^1][^2]を参照のこと。


## 問題設定

実軸上のある区間における $N$ 個の点の集合

$$
\boldsymbol{\Omega} = (\omega_1, \dots, \omega_N)^\mathrm{T}
$$

および近似したい関数 $F(\omega)$ のその各点の値からなる集合

$$
\mathbf{F} = (F(\omega_1), \dots, F(\omega_N))^\mathrm{T}
$$

を入力データとする。

ここで考える問題は、$\|\cdot\|_{\boldsymbol{\Omega}}$ を $\boldsymbol{\Omega}$ 上のL2ノルムとしたとき、ある $\epsilon$ に対して、

$$
\| F(\omega) - \tilde{F}(\omega)\|_{\boldsymbol{\Omega}} < \epsilon, \quad \omega \in \boldsymbol{\Omega}
$$

を満たす

$$
\tilde{F}(\omega) =
\frac{\displaystyle \sum_{k=1}^M \frac{w_k F(\omega_k)}{\omega - \omega_k}} 
{\displaystyle \sum_{k=1}^M \frac{w_k}{\omega - \omega_k}}  
= \sum_{k=1}^{M} \frac{\eta_k}{\omega - z_k},
\quad
w_k \in \mathbb{R},\; z_k, \eta_k \in \mathbb{C}
$$

を見つけることである。上式右辺をbarycentric形式と呼ぶ。また本稿では、最右辺？は極ー留数形式と呼ぶことにする。


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

に分割する。そして、barycentric形式で定義される有理関数

$$
F^{(m-1)}(\omega)
=
\frac{\displaystyle \sum_{k=1}^m \frac{w_k F(\omega_k)}{\omega - \omega_k}}
{\displaystyle \sum_{k=1}^m \frac{w_k}{\omega - \omega_k}}
$$

を考え、重み $\mathbf{w}=(w_1,\dots,w_m)^\mathrm{T}$ に関して、

$$
\text{minimize}
\quad
\| F(\omega) - F^{(m-1)}(\omega)\|_{\boldsymbol{\Omega}}
\quad
\mathrm{s.t.}
\quad
\|\mathbf{w}\|_{m}=1
$$

という線形最小二乗問題を解く。$\|\cdot\|_{m}$ は $m$ 次元ベクトル上のL2ノルムを表す。簡単な操作により、この問題は行列形式

$$
\text{minimize}
\quad
\|\mathbf{A}^{(m)}\mathbf{w}\|_{M-m}
\quad
\mathrm{s.t.}
\quad
\|\mathbf{w}\|_m=1
$$

へと書き換えられる。ここで、Loewner行列 $\mathbf{A}^{(m)}$ は

$$
\mathbf{A}^{(m)} =
\begin{pmatrix}
\dfrac{F(\omega^{(m)}_1)-F(\omega_1)}{\omega^{(m)}_1-\omega_1} & \cdots & \dfrac{F(\omega^{(m)}_1)-F(\omega_m)}{\omega^{(m)}_1-\omega_m}\\[6pt]
\vdots & \ddots & \vdots \\[6pt]
\dfrac{F(\omega^{(m)}_{M-m})-F(\omega_1)}{\omega^{(m)}_{M-m}-\omega_1} & \cdots & \dfrac{F(\omega^{(m)}_{M-m})-F(\omega_m)}{\omega^{(m)}_{M-m}-\omega_m}
\end{pmatrix}
$$

と定義される。次に、$\mathbf{A}^{(m)}$ に対して特異値分解

$$
\mathbf{A}^{(m)} = \mathbf{U} \boldsymbol{\Sigma} \mathbf{V}^*
$$

を行い、右特異ベクトル $\mathbf{V}$ から重み $\mathbf{w}$ は

$$
\mathbf{w} = \mathbf{V}(1:m, m)
$$

のように得られる。

ステップ $m$ において、もし

$$
\| F(\omega) - F^{(m-1)}(\omega)\|_{\boldsymbol{\Omega}} < \epsilon
$$

が満たされたならば反復を終了し、そのときの次数は $M = m-1$ となる。満たされない場合は、貪欲法[^1]によって補集合 $\overline{\boldsymbol{\Omega}^{(m)}}$ から次のサポート点 $\omega_{m+1}$ を選ぶ。

さらに、極 $z_k\;(k=1,\dots,m-1)$ は、 $(m+1)\times(m+1)$ 一般化固有値問題

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

の固有値として評価できる[^1]。留数 $\eta_k \;(k=1,\dots,m-1)$ は、解析関数の商の単純極に対する公式を用いて

$$
\eta_k 
= \frac{\displaystyle \sum_{k'=1}^m \dfrac{w_{k'}\, F(\omega_{k'})}{z_k - \omega_{k'}}}
{\displaystyle \sum_{k'=1}^m \Bigl[-\dfrac{w_{k'}}{(z_k - \omega_{k'})^2}\Bigr]}
$$

によって計算される。


## ソースコード

MATLAB版は
https://github.com/chebfun/chebfun
Python版は
https://github.com/c-f-h/baryrat
Julia版は
https://github.com/complexvariables/RationalFunctionApproximation.jl
でそれぞれ実装されている。

ここでは、RationalFunctionApproximation.jlを使ってみる。例として、関数

$$
f(x) = \exp(-5(x+1)^2) + 2\exp(-10(x-2)^2)
$$

を $[-4,4]$ の範囲で近似することを考える。精度は $\epsilon = 1.0\times 10^{-8}$ とした。以下コードを示す。

```julia
using RationalFunctionApproximation
using CairoMakie

# 近似したい関数
f = x -> exp(-5(x+1)^2) + 2exp(-10(x-2)^2)

# 近似に用いるサンプル点を作成
x = LinRange(-4, 4, 200)

# AAAアルゴリズム
r = aaa(x, f.(x), tol=1e-8)

# 有理関数の次数、極、留数
deg = degree(r)
println("Degree: ", deg)
pol = poles(r)
res = residues(r)

# barycentric型有理関数の関数オブジェクト
f_approx1 = x -> r(x)

# 有理関数の極とその留数から有理関数を作成
f_approx2 = x -> sum(res[i] / (x - pol[i]) for i in 1:deg)

# CairoMakie を用いた可視化
xplot = LinRange(-4, 4, 500)
y_original = f.(xplot)
y_approx1 = f_approx1.(xplot)
y_approx2 = real.(f_approx2.(xplot))

fig = Figure()
ax = Axis(fig[1, 1],
    xlabel = "x",
    ylabel = "f(x)",
    title  = "Rational Approximation"
)

lines!(ax, xplot, y_original, label = "Original function", color = :blue)
lines!(ax, xplot, y_approx1, label = "AAA approximation (barycentric)", color = :orange, linestyle = :dash)
lines!(ax, xplot, y_approx2, label = "AAA approximation", color = :green, linestyle = :dot)
axislegend(ax, position = :lt)
save("aaa-figure.png", fig)

# 誤差を計算、図示
err1 = abs.(f.(xplot) - f_approx1.(xplot))
err2 = abs.(f.(xplot) - f_approx2.(xplot))
println("Max error (barycentric): ", maximum(err1))
println("Max error: ", maximum(err2))

fig = Figure()
ax = Axis(fig[1, 1],
    xlabel = "x",
    ylabel = "Error",
    title  = "Error of Rational Approximation"
)

lines!(ax, xplot, err1, label = "Error (barycentric)", color = :orange)
lines!(ax, xplot, err2, label = "Error", color = :green)
ax2 = Axis(
    fig;
    bbox = BBox(350, 550, 150, 350),
    title = "Enlarged view",
    xgridvisible = false,
    ygridvisible = false,
)
lines!(ax2, xplot, err1, label = "Error (barycentric)", color = :orange)
axislegend(ax, position = :lc)
save("aaa-error.png", fig)
```

結果は、barycentric形式（橙）と 極ー留数形式（緑）の両方の結果を図示した。与えた精度に対して、項数は $M=28$ となった。
![alt text](/images/aaa-figure.png)
どちらもよく近似できていることが分かる。しかし、誤差の絶対値を見てみると
![alt text](/images/aaa-error.png)
barycentric形式ではほとんど指定した精度に収まっているが、極と留数を計算する段階で誤差が大きくなるようなので注意が必要である。

---

[^1]: Y. Nakatsukasa, O. Sete, and L. N. Trefethen, SIAM J. Sci. Comp. 40, A1494-A1522 (2018).
[^2]: Y. Nakatsukasa and L. N. Trefethen, SIAM J. Sci. Comp. 42, A3157-A3179 (2020).
[^3]: A. C. Antoulas and B. D. Q. Anderson, IMA J. Math. Control Inform. 3, 61–88 (1986).


