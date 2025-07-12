---
title: "扁長回転楕円体波動関数（prolate spheroidal wave function）の理論と応用"
emoji: "🍡"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: ["アルゴリズム","julia"]
published: false
---

## 概要

本稿では、


## PSWFの特徴

PSWFは、扁長回転楕円体座標におけるヘルムホルツ方程式を変数分離した微分方程式
$$
L_c[\varphi](x)=-\frac{d}{d x}\left(\left(1-x^2\right) \cdot \frac{d \varphi}{d x}(x)\right)+c^2 x^2 \varphi(x) = 0,\quad -1<x<1
$$
の解として得られる。ここで示した方程式は、回転楕円体座標における偏角の一般的な方程式の特殊な場合であり、解は0次 (order zero) のPSWFとなる。演算子$L_c$は、Sturm-Liouville演算子であり、PSWFはベッセル関数など他の多くのクラスの特殊関数と同様に、古典Sturn-Liouville理論で取り扱われる。

一方で、PSWFの特筆すべき点は演算子$L_c$だけでなく、積分演算子
$F_c: L^2[-1,1] \rightarrow L^2[-1,1]$, 
$$
F_c[\varphi](x)=\int_{-1}^1 \mathrm{d} t\; \varphi(t) \mathrm{e}^{\mathrm{i} c x t} 
$$
の固有関数にもなっている点にある。

演算子の固有値を$\lambda_j,\;(j\in\mathbb{N})$と記すことにすると、各固有値に対して対応する固有関数は
$$
\lambda_j\psi_j(x)=\int_{-1}^1 \mathrm{d} t\;\psi_j(t) \mathrm{e}^{\mathrm{i} c x t} 
$$
を満たす。

## ルジャンドル関数に基づくPSWFの評価

# 固有値の評価

## 求積法


## ソースコード



---

[^1]: B. Moore, IEEE Trans. Autom. Control 26, 17–32 (1981).
[^2]: H. Ikeno, Comput. Phys. Commun. 230, 135–144 (2018).
[^3]: T. Haut and G. Beylkin, SIAM J. Matrix Anal. Appl. 33, 1101–1125 (2012). 



