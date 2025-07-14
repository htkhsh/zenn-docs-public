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
L_c[\varphi](x)=-\frac{d}{d x}\left(\left(1-x^2\right) \cdot \frac{d \varphi}{d x}(x)\right)+c^2 x^2 \varphi(x)=0,\quad -1<x<1
$$

の解として得られる。ここで示した方程式は、回転楕円体座標における偏角の一般的な方程式の特殊な場合であり、解は0次 (order zero) のPSWFとなる。演算子 $L_c$ は、Sturm-Liouville演算子であり、PSWFはベッセル関数など他の多くのクラスの特殊関数と同様に、古典Sturn-Liouville理論で取り扱われる。したがって、演算子 $L_c$ の固有値を $\chi_n\;(n\in\mathbb{N})$ とおくと、PSWFは

$$
\left(1-x^2\right) \psi_n''(x)-2 x \cdot \psi_n'(x)+\left(\chi_n-c^2 x^2\right) \psi_n(x)=0
$$

を満たす。

一方で、PSWFの特筆すべき点は演算子 $L_c$ だけでなく、積分演算子
$F_c: L^2[-1,1] \rightarrow L^2[-1,1]$, 
$$
F_c[\varphi](x)=\int_{-1}^1 \mathrm{d} t\; \varphi(t) \mathrm{e}^{\mathrm{i} c x t} 
$$
の固有関数にもなっている点にある。演算子 $F_c$ の固有値を $\lambda_j\;(j\in\mathbb{N})$ と記すことにすると、各固有値に対して対応するPSWFは
$$
\lambda_j\psi_j(x)=\int_{-1}^1 \mathrm{d} t\;\psi_j(t) \mathrm{e}^{\mathrm{i} c x t} 
$$
を満たす。ここで、$c$ は帯域幅（band limit）と呼ばれるパラメータである。

:::message
**定理1** $c>0$ を実数とし、演算子 $F_c$ が上記の式で定義されるとする。このとき、$F_c$ の固有関数 $\psi_0, \psi_1, \ldots$ は純実数であり、$L^2[-1,1]$ において正規直交かつ完全である。偶数番号の関数は偶関数であり、奇数番号の関数は奇関数である。各関数 $\psi_n$ は区間 $(-1,1)$ 内にちょうど $n$ 個の単純零点を持つ。$F_c$ のすべての固有値 $\lambda_n$ は非零かつ単純である。偶数番号の固有値は純実数であり、奇数番号の固有値は純虚数である。特に、すべての非負整数 $n$ に対して $\lambda_n=i^n\left|\lambda_n\right|$ が成り立つ。さらに、すべての正整数 $n>0$ に対して、関数 $\psi_0, \ldots, \psi_{n-1}$ は区間 $[-1,1]$ 上でチェビシェフ系を構成する。
:::

## ルジャンドル関数とPSWF

### （正規化された）ルジャンドル多項式

古典的なルジャンドル多項式を $P_n$ で表すことにし、3項漸化式

$$
P_{n+1}(x)=\frac{2 n+1}{n+1} x P_n(x)-\frac{n}{n+1} P_{n-1}(x)
$$

と初期条件

$$
P_0(x)=1, \quad P_1(x)=x
$$

によって定義される。またすべての $k=0,1,2, \ldots$ について

$$
P_k(1)=1
$$

が成り立ち、各多項式 $P_k$ は微分方程式

$$
\left(1-x^2\right) P_k''(x)-2 x P_k'(x)+k(k+1) P_k(x)=0
$$

を満たす。

式(8)と(9)で定義される多項式は区間 $[-1,1]$ で直交するが、正規直交ではない。ルジャンドル多項式の正規化版は $\overline{P_n}$ で表され、

$$
\overline{P_n}(x)=P_n(x) \cdot \sqrt{n+1 / 2}
$$

と定義される。すると、各 $n \geqslant 0$ に対して、

$$
\int_{-1}^1\left(\overline P_n(x)\right)^2 \mathrm{~d} x=1
$$

を満たす。

### ルジャンドル多項式によるPSWFの評価

正規化されたルジャンドル多項式は $L^2[-1,1]$ の正規直交基底を構成する。すべての実数 $c>0$ とすべての非負整数 $n \geq 0$ に対して、帯域幅 $c$ に対応する扁長回転楕円体波動関数 $\psi_n$ は以下の級数に展開できる：

$$
\psi_n(x)=\sum_{k=0}^{\infty} \beta_k^{(n)} \cdot \overline{P_k}(x)=\sum_{k=0}^{\infty} \alpha_k^{(n)} \cdot P_k(x),
$$

ここで $-1 \leq x \leq 1$ である。$\beta_0^{(n)}, \beta_1^{(n)}, \ldots$ は次式で定義される：

$$
\beta_k^{(n)}=\int_{-1}^1 \psi_n(x) \cdot \overline{P_k}(x) d x
$$

また、$\alpha_0^{(n)}, \alpha_1^{(n)}, \ldots$ は次式で定義される：

$$
\alpha_k^{(n)}=\beta_k^{(n)} \cdot \sqrt{k+1 / 2}=(k+1 / 2) \cdot \int_{-1}^1 \psi_n(x) \cdot P_k(x) d x
$$

ここで $k=0,1,2, \ldots$ である。

$$
A_{k, k-2} \cdot \beta_{k-2}^{(n)}+A_{k, k} \cdot \beta_k^{(n)}+A_{k, k+2} \cdot \beta_{k+2}^{(n)}=\chi_n \cdot \beta_k^{(n)}
$$

がすべての $k=2,3, \ldots$ について成り立つ。ここで、$A_{k, k}, A_{k+2, k}, A_{k, k+2}$ は以下の式で定義される：

$$
\begin{aligned}
& A_{k, k}=k(k+1)+\frac{2 k(k+1)-1}{(2 k+3)(2 k-1)} \cdot c^2 \\
& A_{k, k+2}=A_{k+2, k}=\frac{(k+2)(k+1)}{(2 k+3) \sqrt{(2 k+1)(2 k+5)}} \cdot c^2
\end{aligned}
$$

ここで $k=0,1,2, \ldots$ である。言い換えると、無限ベクトル $\left(\beta_0^{(n)}, \beta_1^{(n)}, \ldots\right)$ は恒等式

$$
\left(A-\chi_n I\right) \cdot\left(\beta_0^{(n)}, \beta_1^{(n)}, \ldots\right)^T=0
$$

を満たす。ここで、$I$ は無限単位行列であり、無限対称行列 $A$ の非零成分は式(2.53)で与えられる。

行列 $A$ は自然に2つの無限対称三重対角行列 $A^{\text {even }}$ と $A^{\text {odd }}$ に分割される。前者は $A$ の偶数インデックスの行と列からなる要素で構成され、後者は $A$ の奇数インデックスの行と列からなる要素で構成される。さらに、すべての非負整数の組 $n, k \geq 0$ に対して、

$$
\beta_k^{(n)}=0, \quad \text { if } k+n \text { is odd }
$$

が第2.4節の定理2.3と式(2.48)の組み合わせにより成り立つ。以下の定理（文献[73]に少し異なる形で現れる）では、PSWFの評価のための数値アルゴリズムにつながるこれらの観察のいくつかの含意をまとめる。

**定理2.12** $c>0$ を実数とし、無限三重対角対称行列 $A^{\text {even }}$ と $A^{\text {odd }}$ がそれぞれ以下で定義されるとする：

$$
A^{\text {even }}=\left(\begin{array}{ccccc}
A_{0,0} & A_{0,2} & & & \\
A_{2,0} & A_{2,2} & A_{2,4} & & \\
& A_{4,2} & A_{4,4} & A_{4,6} & \\
& & \ddots & \ddots & \ddots
\end{array}\right)
$$

$$
A^{\text {odd }}=\left(\begin{array}{ccccc}
A_{1,1} & A_{1,3} & & & \\
A_{3,1} & A_{3,3} & A_{3,5} & & \\
& A_{5,3} & A_{5,5} & A_{5,7} & \\
& & \ddots & \ddots & \ddots
\end{array}\right),
$$

ここで、成分 $A_{k, j}$ は式(2.53)で定義される。また、無限ベクトル $\beta_{\text {even }}^{(n)} \in l^2$ と $\beta_{\text {odd }}^{(n)} \in l^2$ がそれぞれ以下の式で定義されるとする：

$$
\beta_{\text {even }}^{(n)}=\left(\beta_0^{(n)}, \beta_2^{(n)}, \ldots\right)^T, \quad \beta_{\text {odd }}^{(n)}=\left(\beta_1^{(n)}, \beta_3^{(n)}, \ldots\right)^T
$$

ここで、$\beta_0^{(n)}, \beta_1^{(n)}, \ldots$ は式(2.48)で定義される。$n$ が偶数の場合、

$$
A^{\text {even }} \beta_{\text {even }}^{(n)}=\chi_n \cdot \beta_{\text {even }}^{(n)}
$$

$n$ が奇数の場合、

$$
A^{\text {odd }} \beta_{\text {odd }}^{(n)}=\chi_n \cdot \beta_{\text {odd }}^{(n)}
$$


**備考3.** 無限ベクトル $\beta^{(n)} \in l^2$ を、$n$ が偶数の場合は $\beta_{\text {even }}^{(n)}$ に等しく、$n$ が奇数の場合は $\beta_{\text {odd }}^{(n)}$ に等しいと定義する。この記法では、$\beta^{(0)}, \beta^{(2)}, \ldots$ は $A^{\text {even }}$ の固有ベクトルであり、$\beta^{(1)}, \beta^{(3)}, \ldots$ は $A^{\text {odd }}$ の固有ベクトルである。

**備考4.** 行列 $A^{\text {even }}$ と $A^{\text {odd }}$ は無限行列であり、その成分は行番号や列番号の増加とともに減衰しないが、各固有ベクトル $\beta^{(n)}$ の各成分は超指数的に高速に減衰することが知られている。

## 固有値の評価

上記のPSWF評価アルゴリズムは微分演算子(24)の固有値 $\left\{\chi_j\right\}$ も生成するが、積分演算子 $F_c$（式(17)で定義）の固有値 $\left\{\lambda_j\right\}$ は生成しない。これらの固有値の一部は以下の式を用いて計算できる：

$$
\lambda_j \psi_j(x)=\int_{-1}^1 \mathrm{e}^{\mathrm{i} c x t} \psi_j(t) \mathrm{d} t
$$

右辺の積分を数値的に評価するが、この評価の条件数は明らかに $1 / \lambda_j$ 程度であり、小さな $\lambda_j$ の計算には不適切である。良条件な手順は以下の通りである：

- 式(69)を用いて $\lambda_0$ を計算する。右辺を数値的に評価し、$x=0$ とする（これにより $\psi_0(x)$ が小さくならない）。
- 計算された $\lambda_0$ と系3.2を用いて、$j=1,2, \ldots, m$ に対する絶対値 $\left|\lambda_j\right|$ を計算する。各 $\left|\lambda_j\right|$ を $\left|\lambda_{j-1}\right|$ から計算する（そして再び、必要な積分を数値的に評価する）。
- $\lambda_j=\mathrm{i}^j\left|\lambda_j\right|$ という事実（定理2.4を参照）を用いて計算を完成させる。

## 求積法

## フーリエ変換への応用

フーリエ変換に応用する

## ソースコード



---

[^1]: B. Moore, IEEE Trans. Autom. Control 26, 17–32 (1981).
[^2]: H. Ikeno, Comput. Phys. Commun. 230, 135–144 (2018).
[^3]: T. Haut and G. Beylkin, SIAM J. Matrix Anal. Appl. 33, 1101–1125 (2012). 



