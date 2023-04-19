open Brownian;;
open PolynomialKarhunenLoeveBrownian;;
open List;;

(*
Parabola-ODE method for Inhomogeneous Brownian Motion
\mathrm{d}y_{t} = a(b-y_t)\mathrm{d}t+\sigma y_t\mathrm{d}W_t

Y_0 := y_0
Y_{k+1} := e^{-\tilde{a}h+\sigma W_{t_k,t_{k+1}}} (Y_k + ab\int_{t_k}^{t_{k+1}} e^{\tilde{a}(s-t_k)-\sigma \tilde{W}_{t_k,s}}\mathrm{d}s

Return [Y_0, Y_h, Y_2h,...] where h = t_max / n_t

Generate on the fly the brownian motion, no memory overflow but not reproducible.
*)

let parabola_igbm a b sigma y0 t_max n_t n_int =
    let tildea = a+.0.5*.sigma*.sigma and h = t_max /. float_of_int n_t in
    let sqrth = sqrt h and ds = h /. float_of_int n_int in
    let grid = range 0. ds 1. n_int in
    let rec aux k accu =
        if k <= n_t then
            let path = bm_paths 0. 1. 1. n_int 1 |> hd (*(W_s)_{s\in[0,1]} will be rescaled*) in
            let w1 = path |> rev |> hd in
            let parabola = parabola_brownian n_int path ~w1bool:true ~w1:w1 in
            match accu with
                | yk::_->
                    let integrands = map (fun s -> exp (tildea*.s*.h-.sigma*.sqrth*.(parabola s)) *. ds) grid in
                    let ykp1 = exp (-1.*.tildea*.h+.sqrth *.sigma*.w1)*.(yk+.a*.b*.(fold_left (+.) 0.0 integrands)) in
                aux (k+1) (ykp1::accu)
        else
            rev accu
    in aux 1 [y0]
    ;;



(*
Log-ODE method
Y_{k+1} = Y_k e^{-\tilde{a}h + \sigma W_{t_k,t_{k+1}}} + abh\bigg(1-\sigma H_{t_k,t_{k+1}}+\sigma^2\bigg(3/5h H_{t_k,t_{k+1}}^2+1/30h\bigg)\bigg)\frac{e^{-\tilde{a}h+\sigma W_{t_k,t_{k+1}}}-1}{-\tilde{a}h+\sigma W_{t_k,t_{k+1}}}
*)
let oversqrt6 = 0.40824829046386301637;;
let log_ode_igbm a b sigma y0 t_max n_t n_int =
    let tildea = a+.0.5*.sigma*.sigma and h = t_max /. float_of_int n_t in
    let sqrth = sqrt h and eigen_func = jacobi 2. |> eigen 1. in
    let rec aux k accu =
        if k <= n_t then
            let path = bm_paths 0. 1. 1. n_int 1 |> hd (*(W_s)_{s\in[0,1]} will be rescaled*) in
            let w1 = path |> rev |> hd in
            let i1 = eigen_func |> basis 1. n_int path |> hd in
            let space_time_levy_area = oversqrt6 *. i1 in
            match accu with
                | yk::_->
                    let ykp1 = yk*.(exp (-1.*.tildea*.h+.sigma*.sqrth*.w1)) +. a*.b*.h*.(
                        1. -. sigma *. space_time_levy_area *. sqrth +. sigma*.sigma*.(0.6*.h*.space_time_levy_area*.space_time_levy_area+.1./.30.*.h))*.
                        (exp (-1.*.tildea*.h+.sigma*.sqrth*.w1) -. 1.) /. (-1. *.tildea*.h+.sigma*.sqrth*.w1) in
                aux (k+1) (ykp1::accu)
        else
            rev accu
    in aux 1 [y0]
    ;;

(*
reproducible (given path), can cause stack overflow
slower than log_ode_igbm given equals n_int & n_t
*)
let log_ode_igbm_given_path a b sigma y0 path n_t t_max =
    let tildea = a+.0.5*.sigma*.sigma and n_int = length path / n_t and h = t_max /. float_of_int n_t in
    let sqrth = sqrt h in
    let oversqrth = 1. /. sqrth in
    let splitted_brownian = split_f path n_t in
    let standardized_brownians =
        map (fun w -> let w0 = hd w in map (fun wt -> oversqrth *. (wt -. w0)) w) splitted_brownian in
    let eigen_func = jacobi 2. |> eigen 1. in
    let rec aux accu ongoing_standardized_brownians k =
        if k <= n_t then
            match ongoing_standardized_brownians with
                | current_path::other_paths ->
                    let w1 = current_path |> rev |> hd in
                    let i1 = eigen_func |> basis 1. n_int current_path |> hd in
                    let space_time_levy_area = oversqrt6 *. i1 in
                    match accu with
                        | yk::_->
                            let ykp1 = yk*.(exp (-1.*.tildea*.h+.sigma*.sqrth*.w1)) +. a*.b*.h*.(
                                1. -. sigma *. space_time_levy_area *. sqrth +. sigma*.sigma*.(0.6*.h*.space_time_levy_area*.space_time_levy_area+.1./.30.*.h))*.
                                (exp (-1.*.tildea*.h+.sigma*.sqrth*.w1) -. 1.) /. (-1. *.tildea*.h+.sigma*.sqrth*.w1) in
                            aux (ykp1::accu) other_paths (k+1)
        else
            rev accu
    in aux [y0] standardized_brownians 1
    ;;

(*
reproducible (given path), can cause memory overflow
slower than parabola_igbm given equals n_int & n_t
*)
let parabola_igbm_given_path a b sigma y0 path n_t t_max =
    let tildea = a+.0.5*.sigma*.sigma in
    let n_int = length path / n_t and h = t_max /. float_of_int n_t in
    let ds = h /. float_of_int n_int and sqrth = sqrt h in
    let oversqrth = 1. /. sqrth in
    let splitted_brownian = split_f path n_t and grid = range 0. ds 1. n_int in
    let standardized_brownians =
        map (fun w -> let w0 = hd w in map (fun wt -> oversqrth *. (wt -. w0)) w) splitted_brownian in
    let rec aux accu ongoing_standardized_brownians k =
        if k <= n_t then
            match ongoing_standardized_brownians with
                | current_path::other_paths ->
                let w1 = current_path |> rev |> hd in
                                       let parabola = parabola_brownian n_int current_path ~w1bool:true ~w1:w1 in
                                       match accu with
                                        | yk::_->
                                            let integrands = map (fun s -> exp (tildea*.s*.h-.sigma*.sqrth*.(parabola s)) *. ds) grid in
                                            let ykp1 = exp (-1.*.tildea*.h+.sqrth *.sigma*.w1)*.(yk+.a*.b*.(fold_left (+.) 0.0 integrands)) in
                                        aux (ykp1::accu) other_paths (k+1)
        else
            rev accu
    in aux [y0] standardized_brownians 1
    ;;

(*on the fly*)
Printf.printf "\n parabola_ode \n";;
parabola_igbm 0.1 0.04 0.6 0.06 1. 10000 1000 |> iter pprint;;
Printf.printf "\n log_ode \n";;
log_ode_igbm 0.1 0.04 0.6 0.06 1. 10000 1000 |> iter pprint;;
Printf.printf "\n parabola_ode_given_path \n";;

(*given path*)
let mypath = bm_paths 0.0 1.0 1. 10000000 1 |> hd in (*n_int*n_t*)
parabola_igbm_given_path 0.1 0.04 0.6 0.06 mypath 10000 1. |> iter pprint;
Printf.printf "\n log_ode_given_path \n";
log_ode_igbm_given_path 0.1 0.04 0.6 0.06 mypath 10000 1. |> iter pprint;;