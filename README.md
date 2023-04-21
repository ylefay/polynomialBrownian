# A polynomial decomposition of Brownian motion (Corollary 4.1.4, Theorem 4.1.6, Theorem 4.1.9, Theorem 4.1.10 from Foster. 2020)
Let $W$ be a standard Brownian motion on $[0, 1]$. We have
$$W^n = W_1 e_0 + \sum_{k=1}^{n-1}I_k e_k$$

where $e_0(t) := t$ and the random variables $\{I_k\}$ are independant of $W_1$.

For $1\leq k\leq n-1$, 
$$I_k = \int_0^1 (W_t-W_1e_0(t))\frac{e_k(t)}{t(1-t)}\mathrm{d}t$$

where $e_k(t) = 1/k \sqrt{k(k+1)(2k+1)}P_{k+1}(2t-1)$.

For each $k\geq 2$, we have
$$k(k+2)P_{k+2}(X) = (k+1)(2k+1)XP_{k+1}(X)-k(k+1)P_k(X)$$
and $P_2(X) = 1/4(X^2-1)$, $P_3(X)=1/2X(X^2-1)$.

# Parabola-ODE method
Consider the following Stratonovich SDE:
$$\mathrm{d}y_t = f_0(y_t)\mathrm{d}t + f_1(y_t)\circ \mathrm{d}W_t.$$

The Parabola-ODE method (see theorem 4.3.11) gives the following numerical scheme:
Let $h>0$, $n$ the number of points we consider inside every interval $[kh,(k+1)h]$.
We define $(Y_k)$ the approximation of the solution $(y_0, y_h, \ldots)$ by

$$Y_0 := y_0.
$$

and for each $k\geq 0$,

$$Y_{k+1} := z_n.
$$

For each $k$, we define $z$ (we omit the dependency on $k$) by
$$z_0 = Y_k.$$
For each $0\leq i \leq n-1$,

$$z_{i+1} = z_i + \frac{h}{n}f_0(z_{i})+\frac{1}{n}f_1(z_{i})(W_{kh,(k+1)h}+(6-12\frac{i}{n})H_{kh,(k+1)h}).$$

In practice, we split our brownian motion and normalize each part, that is we consider the mapping

$$W \mapsto ((h^{-1/2}(W_{kh+s}-W_{kh}))_{0\leq s\leq h})_{k}.$$

Then, we use the scaling property of brownian motion. We also use the following equality $H_{t_k,t_{k+1}} = \sqrt{\frac{h}{6}}I_1$.

# Polynomial-ODE method
Same numerical scheme as the previous one but we replace 
$$(W_{kh,(k+1)h}+(6-12\frac{i}{n})H_{kh,(k+1)h})$$ 

by

$$\sqrt{h}(W^n_k((i+1)/n)-W^n_k(i/n)),$$
where $W^n_k$ is the $n$-th degree polynomial corresponding to the $(k+1)-th$ standardized part of our brownian motion.

Since $W^n$ is a polynomial, it is easy to derivate it. One may want to implement composition and derivation of polynomials to compute $\frac{\mathrm{d}W^n_k}{\mathrm{d}u}(i/n)$. 

We have $|W-W^n|_{L^2(\mathbb{P})} = O(n^{-1/2})$.

# Inhomogeneous Geometric Brownian Motion
Consider the following Stratonovich SDE:

$$\mathrm{d}y_t = a(b-y_t)\mathrm{d}t + \sigma y_t\circ\mathrm{d}W_t$$

It is known that it admits a unique strong solution satisfying:

$$y_t = e^{-ah+\sigma W_{s,t}}\bigg(y_s+ab\int_s^t e^{a(u-s)-\sigma W_{s,u}}\mathrm{d}u\bigg)$$

The Parabola-ODE method gives

$$
Y_0 := y_0$$

$$
Y_{k+1} = e^{-ah+\sigma W_{s,t}}(Y_k+ab\int_{t_k}^{t_{k+1}} e^{a(u-t_k)-\sigma \tilde{W}_{t_k,u}}\mathrm{d}u)
$$

The Log-ODE method gives

$$
Y_0 := y_0$$

$$
Y_{k+1} = Y_k e^{-\tilde{a}h + \sigma W_{t_k,t_{k+1}}} + abh\bigg(1-\sigma H_{t_k,t_{k+1}}+\sigma^2\bigg(3/5h H_{t_k,t_{k+1}}^2+1/30h\bigg)\bigg)\frac{e^{-\tilde{a}h+\sigma W_{t_k,t_{k+1}}}-1}{-\tilde{a}h+\sigma W_{t_k,t_{k+1}}}
$$

# Compilation and performance profiling
Please see `http://ocamlverse.net/content/optimizing_performance.html`
```
ocamlopt -g -O3 utils.ml normal.ml brownian.ml polynomial.ml polynomialKarhunenLoeveBrownian.ml parabola_method.ml igbm.ml main.ml -o main
perf record --call-graph=dwarf -- ./main
perf report

hyperfine './main 1000' './main 10000' './main 100000' './main 1000000'
```

1000 points (Y_1, ..., Y_1000), with a Brownian motion sampled {1000, 10000, 100000, 1000000} times.

```
Benchmark 1: ./main 1000
  Time (mean ± σ):       7.5 ms ±   0.6 ms    [User: 5.2 ms, System: 2.5 ms]
  Range (min … max):     6.7 ms …  13.6 ms    270 runs
 
  Warning: Statistical outliers were detected. Consider re-running this benchmark on a quiet PC without any interferences from other programs. It might help to use the '--warmup' or '--prepare' options.
 
Benchmark 2: ./main 10000
  Time (mean ± σ):      26.0 ms ±   0.7 ms    [User: 22.3 ms, System: 3.5 ms]
  Range (min … max):    24.9 ms …  28.7 ms    101 runs
 
Benchmark 3: ./main 100000
  Time (mean ± σ):     335.7 ms ±   3.5 ms    [User: 322.8 ms, System: 10.7 ms]
  Range (min … max):   329.2 ms … 342.1 ms    10 runs
 
Benchmark 4: ./main 1000000
  Time (mean ± σ):     15.826 s ±  0.232 s    [User: 15.631 s, System: 0.085 s]
  Range (min … max):   15.487 s … 16.248 s    10 runs
 
Summary
  './main 1000' ran
    3.49 ± 0.32 times faster than './main 10000'
   44.97 ± 3.92 times faster than './main 100000'
 2120.08 ± 186.04 times faster than './main 1000000'
```