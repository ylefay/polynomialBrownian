(*
For benchmarking purpose
Geometric Brownian Motion:
\mathrm{d}S_t = S_t \mut + S_t \sigma\mathrm{d}W_t
or equivalently
\mathrm{d}S_t =  S_t(\mu-\sigma^2/2)\mathrm{d}t + S_t \sigma \circ \mathrm{d}W_t
*)
open Parabola_method;;

(*
Parabola-ODE method for GBM.
*)
let parabola_gbm_given_path mu sigma s0 path n_t t_max =
    let f0 st = st*.(mu-.sigma*.sigma*.0.5) and f1 st = sigma *. st in
    parabola_given_path path f0 f1 s0 n_t t_max;;
é

