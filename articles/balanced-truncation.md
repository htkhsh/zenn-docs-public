---
title: "Balanced Truncation法による複素指数関数和のモデル縮約"
emoji: "🍡"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: ["アルゴリズム","julia"]
published: false
---

## 概要

本稿では、システム制御理論におけるモデル縮約（Model Reduction）を複素指数関数和に適用し、Balanced Truncation法によってより少ない次数の複素指数関数和で元の関数を近似する手法を紹介する。まずアルゴリズムを概観し、実装例を示す。


## アルゴリズムの導入

以下のアルゴリズムは、論文[1^]を参照した。

すでに複素指数関数和

$$
f(t) = \sum_{k=1}^M c_k \mathrm{e}^{-a_k t}, \quad c_k,a_k\in\mathbb{C},\; \operatorname{Re}(a_k)>0.
$$

が与えられているとする。Balanced truncation法の目的は、与えられた精度 $\epsilon$ に対して、より少ない項数を持つ別の最適な指数和を用いて $f(t)$ を近似することであり、すなわち

$$
g(t)=\sum_{k=1}^{M'} c_k' \mathrm{e}^{-a_k' t}, 
\quad
\left\|f(t)-g(t)\right\| < \epsilon 
$$

が成り立つような $g(t)$ を考える。

$f(t)$ をラプラス変換すると、

$$
f(s) 
= \mathcal{L}\bigl[f(t)\bigr] 
= \int_{0}^\infty \mathrm{d}t\, f(t) \, \mathrm{e}^{-st} 
= \sum_{k=1}^{M} \frac{c_k}{s+a_k}
$$

となる。すると問題は、高次の有理関数を低次の有理関数で近似する最適化問題に変形され、その形は

$$
g(s)
= \mathcal{L}\bigl[g(t)\bigr] = \sum_{k=1}^{M'} \frac{c_k'}{s+a_k'}, 
\quad
\left\|f(s)-g(s)\right\| < \epsilon'
$$

となる。これはすべての $s > 0$ について成立する。BTM の利点は、近似誤差に対して $L^\infty$ の上界を与える点にある。

次に$c_k,a_k$ を用いて、可制御性グラミアン $\mathbf{W}_\mathrm{c} \in \mathbb{C}^{M\times M}$ を構成する。これは

$$
[\mathbf{W}_\mathrm{c}]_{ij} = \frac{\sqrt{c_i c_j^*}}{a_i + a_j^*}
$$

で与えられる。このようにして得られる $\mathbf{W}_\mathrm{c}$ は自己随伴準コーシー行列であり、con-eigenvalue分解

$$
\mathbf{W}_\mathrm{c} 
= \bar{\mathbf{U}} \,\boldsymbol{\Sigma}\, \mathbf{U}^{T}
$$

をもつ。ここで、$\mathbf{U} \in \mathbb{C}^{M\times M}$ はユニタリ行列、$\bar{\mathbf{U}}$はその複素共役、$\boldsymbol{\Sigma} = \mathrm{diag}(\sigma_1,\dots,\sigma_M)$ は $\sigma_1 \ge \sigma_2 \ge \cdots \ge \sigma_M > 0$ をみたす対角行列となる。 これは decompositionと呼ばれる。アルゴリズムは論文[3^]で詳しく紹介されている。

整数 $M'$ が 

$$
2\sum_{i=M'+1}^{M} \sigma_i \leq \epsilon
$$

をみたすとき、$\mathbf{A}' \in \mathbb{C}^{M'\times M'}$ を

$$
\mathbf{A}' = \mathbf{U}_{M'}^{*}\,\mathbf{A}\,\bar{\mathbf{U}}_{M'}
$$

および、$\mathbf{b}' \in \mathbb{C}^{M'}$ を

$$
\mathbf{b}' = \mathbf{U}_{M'}\,\mathbf{b}
$$

によって構成する。ここで、$\mathbf{A} = \operatorname{diag}(a_1, a_2, \ldots, a_M)$、$\mathbf{b} = (\sqrt{c_1}, \sqrt{c_2}, \ldots, \sqrt{c_M})^T$、 $\mathbf{U}_{M'} = \mathbf{U}(1\!:\!M,\,1\!:\!M')$ である。

最後に、$\mathbf{A}'$ の固有値分解

$$
\mathbf{A}' = \mathbf{X} \,\boldsymbol{\Lambda}\,\mathbf{X}^{-1}
$$

を計算すると、

$$
a_{i}'=\boldsymbol{\Lambda}_{ii}, 
\quad
c_i'=(\mathbf{b}_i^{''})^{2} 
\quad
(i=1, \ldots, M')
$$

が得られる。ここで、$\mathbf{b}^{''} = \mathbf{X}^{-1}\mathbf{b}'$ である。


## ソースコード

C++実装は
https://github.com/hydeik/mxpfit
で公開されている。

また筆者によって
https://github.com/DOC-Package/ExpFit.jl
Julia版も実装されている。

ここではExpFit.jlを使ってみる。例として、乱数で生成した複素パラメータ$\{(c_k,a_k)\}$をもつ関数を近似することを考える。精度は $\epsilon = 1.0\times 10^{-3}$ とした。以下コードを示す。

```julia
using LinearAlgebra
using ExpFit
using Random
using CairoMakie
using LaTeXStrings 

function generate_exponent_coefficient_pairs(n::Int)
    Random.seed!(123)
    exponents = Vector{ComplexF64}(undef, n)
    coefficients = Vector{ComplexF64}(undef, n)
    for i in 1:n
        re_exp = rand() * 9.9 + 0.1
        im_exp = randn()  
        exponents[i] = re_exp + im_exp*im
        re_coeff = randn()
        im_coeff = randn()
        coefficients[i] = re_coeff + im_coeff*im
    end
    return exponents, coefficients
end

tmin = 0.0
tmax  = 20.0      
eps = 1e-3
N = 200
t = range(tmin, tmax, length=N)

# 近似する指数関数の生成
a, c = generate_exponent_coefficient_pairs(200)
f = Exponentials(a,c)

# Balanced Truncation法による近似
er = expred(a, c, eps)
print("Approximation order = ", length(er.coeff), "\n")

# 近似値と誤差の計算
fv = f.(t)
erv = er.(t)
err = abs.(erv .- fv)

# プロット
fig = Figure(size = (640, 560))
ax1 = Axis(fig[1, 1], ylabel = L"f(t)", ylabelsize = 25)
ax2 = Axis(fig[2, 1], xlabel = L"t", ylabel = L"\delta f(t)", xlabelsize = 25, ylabelsize = 25)
for (F, col, lbl) in zip((real, imag), (:orangered, :royalblue), (L"\mathrm{Re} f(t)", L"\mathrm{Im} f(t)"))
    lines!(ax1, t, F.(erv), label = lbl, color = col, linewidth = 2.5)
    lines!(ax2, t, F.(erv .- fv),      color = col, linewidth = 2.5)
end
lines!(ax1, t, real.(fv), label = L"Reference", color = :black, linestyle = :dash, linewidth = 2.5)
lines!(ax1, t, imag.(fv), color = :black, linestyle = :dash, linewidth = 2.5)
axislegend(ax1, position = :rb, labelsize = 25)
save("result.png", fig)
```
与えた精度に対して、項数は $M'=14$ となった。以下に結果を図示した。
![alt text](/images/btm.png)
誤差が許容範囲内に収まっており、よく近似できていることが確認できる。

---

[1^] H. Ikeno, Comput. Phys. Commun. 230, 135–144 (2018).
[2^] B. Moore, IEEE Trans. Autom. Control 26, 17–32 (1981).
[3^] T. Haut and G. Beylkin, SIAM J. Matrix Anal. Appl. 33, 1101–1125 (2012). 



