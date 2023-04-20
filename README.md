# A polynomial decomposition of Brownian motion (Corollary 4.1.4, Theorem 4.1.6, Theorem 4.1.9, Theorem 4.1.10 from Foster. 2020)
Let $W$ be a standard Brownian motion on $[0, 1]$. We have
$$W^n = W_1 e_0 + \sum_{k=1}^{n-1}I_k e_k$$

where $e_0(t) := t$ and the random variables $\{I_k\}$ are independant of $W_1$.

For $1\leq k\leq n-1$, 
$$I_k = \int_0^1 (W_t-W_1e_0(t))\frac{e_k(t)}{t(1-t)}\mathrm{d}t$$

where $e_k(t) = 1/k \sqrt{k(k+1)(2k+1)}P_{k+1}(2t-1)$.

For each $k\geq 2$, we have
$$k(k+2)P_{k+2}(X) = (k+1)(2k+1)XP_{k+1}(X)-k(k+1)P_k(X)$$

# Parabola-ODE method
Consider the following Stratonovich SDE:
$$\mathrm{d}y_t = f_0(y_t)\mathrm{d}t + f_1(y_t)\circ \mathrm{d}W_t$$.

The Parabola-ODE method (see theorem 4.3.11) gives the following numerical scheme:
Let $h>0$, $n$ the number of points we consider inside every interval $[kh,(k+1)h]$.
We define $(Y_k)$ the approximation of the solution $(y_0, y_h, \ldots)$ by

$$Y_0 := y_0$$
and for each $k\geq 0$,
$$Y_{k+1} = z_n$$.

For each $k$, we define $z$ by
$$z_0 = Y_k$$
For each $0\leq i \leq n-1$,

$$z_{i+1} = z_i + f_1(z_{i})(W_{kh,(k+1)h}+(6-12\frac{i}{n})H_{kh,(k+1)h}) + \frac{h}{n}f_0(z_{i})$$

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

In practice, we normalize the different part of our brownian motion and use the scaling property of brownian motion as well as
the following equality $H_{t_k,t_{k+1}} = \sqrt{\frac{h}{6}}I_1$.
# Compilation and performance profiling
Please see `http://ocamlverse.net/content/optimizing_performance.html`
```
ocamlopt -g ./polynomial.ml ./brownian.ml ./polynomialKarhunenLoeveBrownian.ml ./igbm.ml 
perf record --call-graph=dwarf -- ./a.out
perf report
```