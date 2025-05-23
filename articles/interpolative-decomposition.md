---
title: "補間分解（Interpolative Decomposition）と積分への応用"
emoji: "🍡"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: ["アルゴリズム","julia"]
published: false
---

## 概要

本稿では、補間分解（Interpolative Decomposition, ID）について解説を行う。


## 補間分解

ネット上の日本語の記事だとこちらで詳しく解説されている。
https://qiita.com/shachah-svaahaa/items/6a6530b273a9d5846d69



### ソースコード

Fortranで

実装されている。
SciPy版
https://github.com/hydeik/mxpfit
もあるが、内部では上記のFortranコードを呼んでいるようである。

## ガウス求積法

例えば、できるだけ少ない次数の離散和で分解したいとする

実関数 $f(t,\omega): \mathbb{R}\times\mathbb{R}\rightarrow \mathbb{R}$ について $\omega$ に関する定積分、より詳細には

$$
f(t) = \int_a^b \mathrm{d}\omega\; G(t,\omega),\quad t \in [c,d]
$$

を考える。

$$
G(t,\omega) = \sum_{i,j} g(t_i,\omega_j) h_i(t) \ell_j(\omega)
$$

これを元の式に代入すると

$f(t) = $

$$
\mathbf{f} = \mathbf{G}\mathbf{w}
$$

---

[^1] H. Cheng, Z. Gimbutas, P. G. Martinsson, and V. Rokhlin, “On the compression of low rank matrices,” SIAM J. Sci. Comput. 26, 1389–1404 (2005).
[2^] P.-G. Martinsson, V. Rokhlin, and M. Tygert, “On interpolation and integration in finite-dimensional spaces of bounded functions,” Commun. Appl. Math. Comput. Sci. 1, 133–142 (2006).
[3^] E. Liberty, F. Woolfe, P.-G. Martinsson, V. Rokhlin, and M. Tygert, “Randomized algorithms for the low-rank approximation of matrices,” Proc. Natl. Acad. Sci. U. S. A. 104, 20167–20172 (2007).
[4^] F. Woolfe, E. Liberty, V. Rokhlin, and M. Tygert, “A fast randomized algorithm for the approximation of matrices,” Applied and Computational Harmonic Analysis 25, 335–366 (2008).
[5^] N. Halko, P. G. Martinsson, and J. A. Tropp, “Finding Structure with Randomness: Probabilistic Algorithms for Constructing Approximate Matrix Decompositions,” SIAM Review 53, 217–288 (2011).
[6^] J. Kaye, K. Chen, and O. Parcollet, “Discrete Lehmann representation of imaginary time Green’s functions,” Phys. Rev. B 105, 235115 (2022).


