open Brownian;;
open PolynomialKarhunenLoeveBrownian;;
open List;;
(*
Parabola-ODE method for Inhomogeneous Brownian Motion
\mathrm{d}y_{t} = a(b-y_t)\mathrm{d}t+\sigma y_t\mathrm{d}W_t

Y_0 := y_0
Y_{k+1} := e^{-\tilde{a}h+\sigma W_{t_k,t_{k+1}}} (Y_k + ab\int_{t_k}^{t_{k+1}} e^{\tilde{a}(s-t_k)-\sigma \tilde{W}_{t_k,s}}\mathrm{d}s

Return [Y_0, Y_h, Y_2h,...] where h = t_max / n_t
*)

let parabola_igbm a b sigma y0 t_max n_t n_int =
    let tildea = a+.0.5*.sigma*.sigma in
    let h = t_max /. float_of_int n_t in
    let ds = h /. float_of_int n_int in
    let grid = range 0. ds h n_int in
    let paths = bm_paths 0. 1. 1. n_int n_t (*((W_s)_{s\in[0,1]})_{k<=n_int}, will be rescaled*) in
    let parabolas = map (parabola_brownian n_int) paths in
    let rec aux k parabolas paths accu =
        if k <= n_t then
            match accu,paths,parabolas with
                | yk::_,path::q_paths,parabola::q_parabolas ->
                    let w1 = path |> rev |> hd in
                    let integrands = map (fun s -> exp (tildea*.s-.sigma*.(sqrt h)*.(parabola s)) *. ds) grid in
                    let ykp1 = exp (-1.*.tildea*.h+.(sqrt h)*.sigma*.w1)*.(yk+.a*.b*.(fold_left (+.) 0.0 integrands)) in
                aux (k+1) q_parabolas q_paths ([ykp1]@accu)
        else
            rev accu
    in aux 1 parabolas paths [y0]
    ;;

parabola_igbm 0.1 0.04 0.6 0.06 5. 200 10000 |> iter pprint;;