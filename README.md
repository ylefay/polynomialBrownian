# A polynomial decomposition of Brownian motion (Corollary 4.1.4, Theorem 4.1.6, Theorem 4.1.9, Theorem 4.1.10, Foster. 2020)
Let $W$ be a standard Brownian motion on $[0, 1]$. Using the Karhunen-Loève theorem for some well-chosen $L^2$ space, we have
$$W^n = W_1 e_0 + \sum_{k=1}^{n-1}I_k e_k$$

where $e_0(t) := t$ and the random variables $\{I_k\}$ are independent centred normal variables with variance $1/(k+1)$ satisfying,

$$I_k = \int_0^1 (W_t-W_1e_0(t))\frac{e_k(t)}{t(1-t)}\mathrm{d}t$$

where $e_k(t) = 1/k \sqrt{k(k+1)(2k+1)}P_{k+1}(2t-1)$,

and where for any $k\geq 2$,
$$k(k+2)P_{k+2}(X) = (k+1)(2k+1)XP_{k+1}(X)-k(k+1)P_k(X)$$
and $P_2(X) = 1/4(X^2-1)$, $P_3(X)=1/2X(X^2-1)$.

# From Stratonovich SDE to Ito SDE
Consider the following Ito SDE:
$$\mathrm{d}y_t = f_0(y_t)\mathrm{d}t + f_1(y_t)\mathrm{d}W_t.$$

Its Stratonovich counterpart is
$$\mathrm{d}y_t = \bar{f_0}(y_t)\mathrm{d}t + f_1(y_t)\circ\mathrm{d}W_t.$$

where 
$$\bar{f_0} = f_0 - \frac{1}{2}f_1 \partial_x f_1.$$

# Parabola-ODE method
Consider the following Stratonovich SDE:
$$\mathrm{d}y_t = f_0(y_t)\mathrm{d}t + f_1(y_t)\circ \mathrm{d}W_t.$$

Let $t_k = kh$ for some $h>0$, the corresponding sequence of parabola ODEs (Theorem 4.3.11) is

$$
\frac{\mathrm{d}z_k}{\mathrm{d}u} = f_0(z_k)h + f_1(z)(W_{t_k,t_{k+1}} + (12u-6)H_{t_k,t_{k+1}}),
$$

where $W$ and $H$ are obtained either by simulating random normal variables or by computing pathwise integrals.
This gives the following numerical scheme, let $n>0$,
We obtain $(Y_k)$ an approximation of $y$ at $0, h, \ldots$ by applying the Euler scheme to the previous sequence of ODEs, i.e

$$Y_0 := y_0.
$$

and for each $k\geq 0$,

$$Y_{k+1} := z_{k,n}
$$

where 
$$z_{k,0} = Y_k.$$
and each $0\leq i \leq n-1$,

$$z_{k,i+1} = z_{k,i} + \frac{h}{n}f_0(z_{k,i})+\frac{1}{n}f_1(z_{k,i})(W_{kh,(k+1)h}+(12\frac{i}{n}-6)H_{kh,(k+1)h}).$$

# Polynomial-ODE method
Same numerical scheme as the previous one but we replace 
$$(W_{kh,(k+1)h}+(12\frac{i}{n}-6)H_{kh,(k+1)h})$$ 

by

$$\sqrt{h}(\frac{\mathrm{d}W^n_k}{\mathrm{d}u}(i/n))$$
where $W^n_k$ is the $n$-th degree polynomial corresponding to the $(k+1)-th$ standardized part of our brownian motion.


# Log-ODE method
The Log-ODE numerical scheme stems from applying Euler scheme to the following sequence of ODEs, for $k\geq 0$,

$$\frac{\mathrm{d}z}{\mathrm{d}u} = f_0(z)h + f_1(z)W_{t_{k},t_{k+1}} + \[f_1,f_0\](z) h H_{t_k,t_{k+1}} + \[f_1, \[f_1,f_0\]\](z)(0.6hH^2_{t_k,t_{k+1}}+1/30h^2).$$

# Inhomogeneous Geometric Brownian Motion
Consider the following Stratonovich SDE:

$$\mathrm{d}y_t = a(b-y_t)\mathrm{d}t + \sigma y_t\circ\mathrm{d}W_t.$$

It is known that it admits a unique strong solution satisfying:

$$y_t = e^{-ah+\sigma W_{s,t}}\bigg(y_s+ab\int_s^t e^{a(u-s)-\sigma W_{s,u}}\mathrm{d}u\bigg).$$

The parabola ODE method gives

$$
Y_0 := y_0,$$

$$
Y_{k+1} = e^{-ah+\sigma W_{s,t}}(Y_k+ab\int_{t_k}^{t_{k+1}} e^{a(u-t_k)-\sigma \tilde{W}_{t_k,u}}\mathrm{d}u).
$$

The Log-ODE method gives

$$
Y_0 := y_0,$$

$$
Y_{k+1} = Y_k e^{-\tilde{a}h + \sigma W_{t_k,t_{k+1}}} + abh\bigg(1-\sigma H_{t_k,t_{k+1}}+\sigma^2\bigg(3/5h H_{t_k,t_{k+1}}^2+1/30h\bigg)\bigg)\frac{e^{-\tilde{a}h+\sigma W_{t_k,t_{k+1}}}-1}{-\tilde{a}h+\sigma W_{t_k,t_{k+1}}}.
$$

# Compilation and performance profiling
Please see `http://ocamlverse.net/content/optimizing_performance.html`
```
ocamlopt -g -O3 utils.ml rand.ml polynomial.ml polynomialKarhunenLoeveBrownian.ml parabola_method.ml log_method.ml igbm.ml gbm.ml main.ml -o main
perf record --call-graph=dwarf -- ./main
perf report
```


```
hyperfine './main 1000 10' './main 10000 100' './main 100000 100' './main 1000000 1000' './main 10000000 1000' './main 10000000 100' 
Benchmark 1: ./main 1000 10
  Time (mean ± σ):       6.3 ms ±   0.4 ms    [User: 4.6 ms, System: 2.0 ms]
  Range (min … max):     5.7 ms …   8.2 ms    286 runs
 
Benchmark 2: ./main 10000 100
  Time (mean ± σ):      43.4 ms ±   2.2 ms    [User: 40.5 ms, System: 2.7 ms]
  Range (min … max):    40.6 ms …  51.8 ms    66 runs
 
Benchmark 3: ./main 100000 100
  Time (mean ± σ):     399.5 ms ±   9.1 ms    [User: 392.5 ms, System: 4.4 ms]
  Range (min … max):   386.4 ms … 416.2 ms    10 runs
 
Benchmark 4: ./main 1000000 1000
  Time (mean ± σ):      3.954 s ±  0.062 s    [User: 3.909 s, System: 0.020 s]
  Range (min … max):    3.843 s …  4.036 s    10 runs
 
Benchmark 5: ./main 10000000 1000
  Time (mean ± σ):     41.077 s ±  2.743 s    [User: 40.583 s, System: 0.155 s]
  Range (min … max):   38.234 s … 47.733 s    10 runs
 
Benchmark 6: ./main 10000000 100
  Time (mean ± σ):     48.495 s ±  5.526 s    [User: 47.326 s, System: 0.381 s]
  Range (min … max):   42.192 s … 57.308 s    10 runs
 
Summary
  './main 1000 10' ran
    6.85 ± 0.56 times faster than './main 10000 100'
   62.98 ± 4.26 times faster than './main 100000 100'
  623.31 ± 40.89 times faster than './main 1000000 1000'
 6475.67 ± 597.51 times faster than './main 10000000 1000'
 7645.08 ± 997.98 times faster than './main 10000000 100'
``` 

# Reference
Foster, JM. 2020. “Numerical Approximations for Stochastic Differential Equations.” PhD thesis, University of Oxford.
Foster, JM. Lyons, T. H. Oberhauser. 2019. "An optimal polynomial approximation of Brownian motion."