(*
For benchmarking purpose
Geometric Brownian Motion:
\mathrm{d}S_t = S_t \mut + S_t \sigma\mathrm{d}W_t
or equivalently
\mathrm{d}S_t =  S_t(\mu-\sigma^2/2)\mathrm{d}t + S_t \sigma \circ \mathrm{d}W_t
*)
open List;;
open Parabola_method;;

(*
Parabola-ODE method for GBM.
*)
let parabola_gbm_given_path mu sigma s0 path n_t t_max =
    let f0 st = st*.(mu-.sigma*.sigma*.0.5) and f1 st = sigma *. st in
    parabola_given_path path f0 f1 s0 n_t t_max;;

let gbm_given_path mu sigma s0 path t_max =
    (*
    S_{t+\mathrm{d}t} = S_t e^{(\mu-\sigma^2/2)\mathrm{d}t + \sigma \mathrm{d}W_t}
    *)
    let dt = t_max/.(float_of_int (length path)) in
    let rec aux brownian_path accu prevwt = match brownian_path with
        | wt::rest_brownian_path -> let st = match accu with
            | st::q -> st
            | [] -> s0
            in
            aux rest_brownian_path ((st*. exp ((mu-.sigma*.sigma*.0.5)*.dt)+.sigma*.(wt-.prevwt))::accu) wt
        | [] -> rev accu
    in aux path [s0] 0.0;;


