# Mersenne Twister

## 2.1 MTの説明
この論文では, 太字の文字, 例えば $\mathbf{\mathbf{\mathbf{x}}}$ や $\mathbf{a}$ などは, $2$ 要素体 $\mathbb{F_2} = {\{0,1}\}$ 上の $w$ 次元の行ベクトルを示す. これらは, $w$ サイズのマシンワードと同一視される（最下位ビットが右にある）. 
MT（メルセンヌ・ツイスタ）アルゴリズムは, $0$ から $2^{w} - 1$ までの一様擬似乱数整数の系列を生成する. $2^{w} - 1$ で割ることで, 各単語ベクトルを $[0,1]$ の実数と見なす. 
このアルゴリズムは, 以下の線形再帰に基づいている. 

$$
\mathbf{\mathbf{\mathbf{x}}}_{k+n} := \mathbf{\mathbf{\mathbf{x}}}_{k+m} \oplus (\mathbf{\mathbf{\mathbf{x}}}^u_{k} | \mathbf{\mathbf{\mathbf{x}}}^l_{k+1})A\quad(k = 0, 1, \ldots)　\qquad (2.1)
$$

記号の説明：いくつかの定数がある. 再帰の次数を示す整数 $n$, $r$ ($\mathbf{\mathbf{\mathbf{x}}}^{u}_{k}$ の定義に隠されている) という整数, $0 \le r \le w - 1$, 整数 $m$, $1 \le m \le n$, および $\mathbb{F}_2$ 内のエントリを持つ定数 $w \times w$ 行列 $A$ がある. 初期シードとして $\mathbf{\mathbf{x}}_0, \mathbf{\mathbf{x}}_1, \ldots , \mathbf{\mathbf{x}}_{n-1}$ を与える. その後, ジェネレータは $k = 0$ の上記の再帰によって $\mathbf{\mathbf{x}}_n$ を生成する. $k = 1, 2, \ldots$ と設定することで, ジェネレータは $\mathbf{\mathbf{x}}_{n+1}, \mathbf{\mathbf{x}}_{n+2}, \ldots$ を決定する. 再帰の右辺において, $\mathbf{\mathbf{x}}^{u}_k$ は $\mathbf{\mathbf{x}}_k$ の上位 $w - r$ ビットを意味し, $\mathbf{\mathbf{x}}^l_{k+1}$ は $\mathbf{\mathbf{x}}_{k+1}$ の下位 $r$ ビットを意味する. したがって, もし $\mathbf{\mathbf{x}}$ が $(x_{w-1}, x_{w-2}, \ldots , x_0)$ の場合, 定義により, $\mathbf{\mathbf{x}}^u$ は $w-r$ ビットのベクトル $(x_{w-1}, \ldots , x_r)$ であり, $\mathbf{\mathbf{x}}^l$ は $r$ ビットのベクトル $(x_{r-1} , \ldots , x_0)$ である. したがって, ($\mathbf{\mathbf{x}}^u | \mathbf{\mathbf{x}}^{l}_{k+1}$) は単に連結です. すなわち, $\mathbf{\mathbf{x}}_k$ の上位 $w - r$ ビットと $\mathbf{\mathbf{x}}_{k+1}$ の下位 $r$ ビットをその順序で連結して得られる 単語ベクトルである. その後, このベクトルに行列 $A$ を右から乗算する. 最後に, このベクトルに $\mathbf{\mathbf{x}}_{k+m}$ を加算する（$\oplus$ は二進数の加算を示し, $2$ で割った余りが結果になる）. そして次のベクトル $\mathbf{\mathbf{x}}_{k+n}$ を生成する. 
複雑な再帰式 (1) を選択した理由は, 第3.1節で明らかになる. ここでは, もし $r = 0$ ならば, この再帰は, 松本と栗田 [1992; 1994] で提案された以前の TGFSR に簡約され, $r = 0$ かつ $A = I$ ならば, GFSR [Lewis and Payne 1973] に簡約される. 
行列 $A$ の形式を選択する理由は, $A$ による乗算が非常に速いためである. 以下は候補である：

$$
A = \begin{pmatrix}
& 1 &  &  &  & \\
&  & 1 &  &  & \\
&  &  & \ddots &  & \\
&  &  &  & 1 \\
 a_{w-1} & a_{w-2} & \cdots & \cdots & a_0
\end{pmatrix}
$$

そのため, $\mathbf{\mathbf{x}}A$ の計算はビット演算だけで行うことができる. 

$$
\begin{equation}
\mathbf{x}A =
    \begin{cases}
      \textrm{shiftright}(\mathbf{x}) & \text{if $x_{0} = 0$} \\
      \textrm{shiftright}(\mathbf{x}) \oplus \mathbf{a} & \text{if $x_{0} = 1$}
    \end{cases}
\end{equation}
$$

ここで, $\mathbf{a} = (a_{w-1},a_{w-2},\ldots,a_0)$ , $\mathbf{\mathbf{x}} = (x_{w-1},x_{w-2}, \ldots,x_0)$ である. また, 再帰式 (2.1) の $\mathbf{\mathbf{x}}^{u}_k$ と $\mathbf{\mathbf{x}}^{l}_{k+1}$ はビットごとの AND 演算で計算できる. したがって, 再帰式 (2.1) の計算は, ビットシフト, ビットごとの排他的論理和, ビットごとの論理和, およびビットごとの AND 演算で実現される. 
$k -$ 分布を $v$ ビットの精度で改善するために, 各生成された単語に適切な $w \times w$ の可逆行列 $T$ を右から掛ける（松本と栗田 [1994] でテンパリングと呼ばれている）. テンパリング行列 $\mathbf{\mathbf{x}} \mapsto \mathbf{z} = \mathbf{\mathbf{x}}T$ について, 次の連続的な変換を選択した：

$$
\begin{align}
\mathbf{y} &:= \mathbf{x} \oplus (\mathbf{x} \gg u) \\
\mathbf{y} &:= \mathbf{y} \oplus ((\mathbf{y} \ll s) ~ \textrm{AND} ~ \mathbf{b}) \\
\mathbf{y} &:= \mathbf{y} \oplus ((\mathbf{y} \ll t) ~ \textrm{AND} ~ \mathbf{c}) \\
\mathbf{y} &:= \mathbf{y} \oplus (\mathbf{y} \gg l)
\end{align}
$$

ここで, $l$ , $s$ , $t$ , および $u$  は整数であり, $\mathbf{b}$ と $\mathbf{c}$ は適切なワードサイズのビットマスクである. また, $(\mathbf{\mathbf{x}} \gg u)$ は $u$ ビット右シフトを表し, $(\mathbf{\mathbf{x}} \ll u)$ は $u$ ビット左シフトを表す. 変換 $(2.3)$ と $(2.4)$ は松本と栗田 [1994] で使用されたものと同じである. 変換 $(2.2)$ と $(2.5)$ は, MT が最下位ビットを改善できるように追加されている. 
再帰式 (2.1) を実行するために, 以下のように, $n$ ワードの配列を作業領域として取ることが十分である. $\mathbf{\mathbf{x}}[0 : n - 1]$ をワードサイズの $n$ 個の符号なし整数の配列とし, $i$ を整数変数, $\mathbf{u}$ , $\mathbf{ll}$ , $\mathbf{a}$ をワードサイズの符号なし定数整数とします. 

$\textbf{Step 0.}$

$$
\mathbf{u} \gets \underbrace{1 \cdots 1}_{w-r} \> \underbrace{0 \cdots 0}_{r} \quad ;(上位の w-r ビット用のビットマスク) \\
\qquad \mathbf{ll} \gets \underbrace{0 \cdots 0}_{w-r} \> \underbrace{1 \cdots 1}_{r} \quad ;(下位の r ビット用のビットマスク) \\
\qquad \mathbf{a} \gets a_{w-1} a_{w-2} \cdots a_{1} a_{0} \quad ;(行列 A の最終行) \\
$$

$\textbf{Step 1.}$

$$
i \gets 0 \qquad \mathbf{x}[0], \mathbf{x}[1], \ldots , \mathbf{x}[n-1] \gets "任意の非ゼロの初期値"
$$

$\textbf{Step 2.}$

$$
\mathbf{y} \gets (\mathbf{x}[i] ~ \textrm{AND} ~ \mathbf{u}) ~ \textrm{OR} ~ (\mathbf{x}[(i + 1) ~ \textrm{mod}] ~ \textrm{AND} ~ \mathbf{ll});((\mathbf{x}^u_{i} | \mathbf{x}^l_{i + 1})を計算する)
$$

$\textbf{Step 3.}$

$$
 \mathbf{x}[i] \gets \mathbf{x}[(i + m) ~ \textrm{mod} ~ n] ~ \textrm{XOR} ~ (y \gg 1) \\
\begin{equation}
\textrm{XOR} =
    \begin{cases}
      0 \qquad \mathbf{y} の最下位ビットが ~0~ の場合 \\
      \mathbf{a} \qquad \mathbf{y} の最下位ビットが ~1~ の場合
    \end{cases}
\end{equation}
$$

$\textbf{Step 4.} \quad (\mathbf{x}[i]T)を計算する$

$$
\mathbf{y} \gets \mathbf{x}[i] \\
\mathbf{y} \gets \mathbf{y} ~ \textrm{XOR} ~ (\mathbf{y} \gg u) ~ ;(\mathbf{y} を u ビット右にシフトし \mathbf{y} に加算) \\
\mathbf{y} \gets \mathbf{y} ~ \textrm{XOR} ~ (\mathbf{y} \ll s) ~ \textrm{AND} ~ \mathbf{b} \\
\mathbf{y} \gets \mathbf{y} ~ \textrm{XOR} ~ (\mathbf{y} \ll t) ~ \textrm{AND} ~ \mathbf{c} \\
\mathbf{y} \gets \mathbf{y} ~ \textrm{XOR} ~ (\mathbf{y} \gg l) \\
\mathbf{y} を出力
$$

$\textbf{Step 5.} \quad i \gets (i + 1) ~ \textrm{mod} ~ n$

$\textbf{Step 6.} \quad \textrm{Step.2} に進む $

一度に配列全体を書き換えることで, モジュロ $n$ 演算を省略することができる. したがって, 非常に高速な操作のみが必要（付録 $C$ のコードを参照してください）. 
以下の2つのパラメータクラスがあります： $(1)$ 周期を決定する周期パラメータ：整数パラメータ $w$ (ワードサイズ) , $n$ (再帰の次数), $m$ (中間項), $r$ (1つのワードの分離点), およびベクトルパラメータ $\mathbf{a}$ (行列 $A$ )； $(2)$ $v$ ビット精度の $k$-分布のテンパリングパラメータ：整数パラメータ $l$ , $u$ , $s$ , $t$ とベクトルパラメータ $\mathbf{b}$ , $\mathbf{c}$. 

[参考ページ](https://dl.acm.org/doi/pdf/10.1145/272991.272995)

## memo

メルセンヌ・ツイスタ法の乱数の周期は, $2^{nw - r} - 1$ の最大周期を得ることができる.
ここで $2^{nw -r} - 1 = 2^{624\times 32 - 31} - 1 = 2^{19937} - 1$ のようにメルセンヌ素数を用いることで,超長周期を実現できる.

### メルセンヌ数・メルセンヌ素数
#### 定義
$$
M_n = 2^n - 1, \quad n \in \mathbb{N}
$$

上記のように自然数 $n$ を用いて $M_n = 2^n - 1$ で表される数をメルセンヌ数という.具体的には, 

$$
2^1 - 1 = 1 \\
2^2 - 1 = 3 \\
2^3 - 1 = 7 \\
2^4 - 1 = 15 \\
2^5 - 1 = 31 \\
2^6 - 1 = 63 \\
2^7 - 1 = 127 \\
2^8 - 1 = 255 \\
2^9 - 1 = 511 \\
2^{10} - 1 = 1023 \\
2^{11} - 1 = 2047 \\
\vdots
$$

上記に基づいて $M_n$ に素数が出てくる場合をメルセンヌ素数という. $M_n$ が素数である場合, $n$ は素数であるが, 逆は成り立たない.

### 完全数

メルセンヌ数 $M_n = 2^n - 1$ が素数の時, $2^{n-1}(2^n - 1)$ は完全数である. 完全数は, 正の約数の和に等しい自然数を定めた数である. 

$n = 2, M_n = 3$
$$
\begin{align}
2^{n-1}(2^n - 1)
&= 2^{2-1}(2^2 - 1) \\
&= 2 \times 3 \\
&= 6 \\
1 + 2 + 3 &= 6
\end{align}
$$

$n = 3, M_n = 7$
$$
\begin{align}
2^{n-1}(2^n - 1)
&= 2^{3-1}(2^3 - 1) \\
&= 4 \times 7 \\
&= 28 \\
1 + 2 + 4 + 7 + 14 &= 6
\end{align}
$$

[参考サイト](https://www.hello-statisticians.com/explain-terms-cat/mersenne_twister1.html)

### メルセンヌツイスタの構造

$\mathbf{x}A$ の計算は, 有限体 $\mathbb{F}_2$ で行う. すなわち, 足し算は $\textrm{XOR}$ 演算であり, 掛け算は $\textrm{AND}$ 演算である. $\mathbf{x}A$ の計算は, $w^2$ 回ほどの演算が必要になるが, $A$ の構造を利用して次のように高速に計算することが可能.

$$
\begin{equation}
\mathbf{x}A =
    \begin{cases}
      (\mathbf{x} \gg 1) & \text{$x_{0} = 0$} \\
      (\mathbf{x} \gg 1) \oplus \mathbf{a} & \text{$x_{0} = 1$}
    \end{cases}
\end{equation}
$$

ただし, 

$$
\begin{align}
\mathbf{x} &= (x_{w-1}, x_{w-2}, \ldots ,x_0) \\
\mathbf{a} &= (a_{w-1}, a_{w-2}, \ldots ,a_0)
\end{align}
$$

である.
ここで, $\mathbf{x}_i$ の分布を改善するために,  の正則行列 $T$ を右から掛ける. この工程はTemperingと呼ばれる. ここでも, $T$ を愚直に掛けることはせず, $T$ を掛けることに相当する演算として, $u, s, t, l, \mathbf{b}, \mathbf{c}$ パラメータ  に対して以下を行う.

$$
\begin{align}
\mathbf{y} &:= \mathbf{x} \oplus (\mathbf{x} \gg u) \\
\mathbf{y} &:= \mathbf{y} \oplus ((\mathbf{y} \ll s) ~ \textrm{AND} ~ \mathbf{b}) \\
\mathbf{y} &:= \mathbf{y} \oplus ((\mathbf{y} \ll s) ~ \textrm{AND} ~ \mathbf{c}) \\
\mathbf{y} &:= \mathbf{y} \oplus (\mathbf{y} \gg l)
\end{align}
$$

例えば, $\mathbf{y} := \mathbf{x} \oplus (\mathbf{x} \gg u)$ は $\mathbf{x}$ に

$$
\begin{pmatrix}
1 & 0 & \cdots & 0 &  & \cdots & 0 \\
0 & 1 & \cdots & 0 &  & \cdots & 0 \\
  & \vdots & \ddots & 0 &  & \cdots & 0 \\
1 & 0 &  & 1 & 0 & \cdots & 0 \\
0 & 1 & \cdots & 0 & \ddots &  & 0 \\
  & \vdots & \ddots & 0 &  & \ddots &  0 \\
0 & 0  & \cdots & 1 & 0 & \cdots & 1 \\
\end{pmatrix}
$$

という形式の行列を右からかけることに対応している. 他の指揮に関しても, 左下の $1$ の列が右上に移動したり, この $1$ の列に $0$ が混じるような列となる. 各行列は明らかに $\textrm{rank}$ が $w$ で正則なので, $T$ は正則である.

[参考サイト](https://acompany.tech/privacytechlab/pseudorandom-number_1)

----

### GFSR
$\text{LFSR}$ の数列をベクトル列に変えた手法.
$\text{GFSR}$ は $\text{General Feedbacked Shift Register}$ の略.
整定値 $p,q$ をうまく選ぶと最大周期 $2^p - 1$　となる.
乱数性に問題あり(特にランダムウォークで)
$$
\mathbf{x}_{n+p} \equiv \mathbf{x}_{n+q} + \mathbf{x}_n \qquad (\text{mod 2}, ~ q<p)
$$

$p = 5, q = 2$ の時, 以下のよに定める.

$$
\begin{align}
\mathbf{x}_{0} &=
\begin{pmatrix}
1 & 1 & 0 & 1 & 0 \\
\end{pmatrix} \\
\mathbf{x}_{1} &=
\begin{pmatrix}
1 & 1 & 0 & 1 & 0 \\
\end{pmatrix} \\
\mathbf{x}_{2} &=
\begin{pmatrix}
1 & 1 & 0 & 1 & 1 \\
\end{pmatrix}\\
\mathbf{x}_{3} &=
\begin{pmatrix}
1 & 1 & 1 & 0 & 0 \\
\end{pmatrix}\\
\mathbf{x}_{4} &=
\begin{pmatrix}
1 & 0 & 0 & 1 & 1 \\
\end{pmatrix}\\
\end{align}
$$

$$
\begin{align}
\mathbf{x}_{0+5}
&= \mathbf{x}_{0+2} + \mathbf{x}_{0} \\
&= \begin{pmatrix} 1 & 1 & 0 & 1 & 1 \end{pmatrix} + \begin{pmatrix} 1 & 1 & 0 & 1 & 0 \\ \end{pmatrix} \\
&= \begin{pmatrix} 0 & 0 & 0 & 0 & 1 \end{pmatrix}
\end{align}
$$

$$
\begin{align}
\mathbf{x}_{1+5}
&= \mathbf{x}_{1+2} + \mathbf{x}_{1} \\
&= \begin{pmatrix} 1 & 0 & 0 & 0 & 1 \end{pmatrix} + \begin{pmatrix} 1 & 1 & 1 & 0 & 0 \\ \end{pmatrix} \\
&= \begin{pmatrix} 0 & 1 & 1 & 0 & 1 \end{pmatrix}
\end{align}
$$


### TGFSR
$\text{Twisted GFSR}$ 法では, $\text{GFSR}$ 法の漸化式の一部を $\text{Twister}$ と呼ぶ $\mathbb{F}_2$ 係数正方行列 $A$ をかけることで最大周期が $2^p - 1$ である $\text{GFSR}$ 法に比べて, 最大周期は $2^{wp} - 1$ にすることができる.

$$
\mathbf{x}_{n+p} \equiv \mathbf{x}_{n+q} + \mathbf{x}_{n}A \qquad (\text{mod 2}, ~ q<p)
$$

$$
A = \begin{pmatrix}
& 1 &  &  &  & \\
&  & 1 &  &  & \\
&  &  & \ddots &  & \\
&  &  &  & 1 \\
 a_{w-1} & a_{w-2} & \cdots & \cdots & a_0
\end{pmatrix}
$$

### Mersenne Twister
$$
\mathbf{x}_{n+p} = \mathbf{x}_{n+q} + \mathbf{x}_{n+1}B + \mathbf{x}_{n}C
$$
$B, C$ は, $\mathbb{F}_2$ 係数の $(32\times 32)$ 行列
$B = O, C = I$ の場合, $\text{GFSR}$ となる.
$B = O, C = A$ の場合, $\text{TGFSR}$ となる.

- $\text{MT}$ の方が帰規則における状態ベクトル間の複雑な相互作用により周期が長い.
- $\text{MT}$ の方がランダムな数列を生成する. 再帰規則が状態ベクトルをより複雑な方法で混合するため
- $\text{MT}$ の方が計算効率がいい.

再帰規則は, ある要素を, それよりも前の要素を用いて定義する数学的な規則
例:フィボナッチ数列

$F_{0} = 0 \\
F_{1} = 1 \\
F_{n+2} = F_{n} + F_{n+1} \quad(n \geq 0)$