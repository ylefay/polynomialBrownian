open Rand;;
open Utils;;
open PolynomialKarhunenLoeveBrownian;;
open Parabola_method;;

(*
Parabola-ODE method for Inhomogeneous Brownian Motion
\mathrm{d}y_{t} = a(b-y_t)\mathrm{d}t+\sigma y_t\mathrm{d}W_t

Y_0 := y_0,
Y_{k+1} := e^{-\tilde{a}h+\sigma W_{t_k,t_{k+1}}} (Y_k + ab\int_{t_k}^{t_{k+1}} e^{\tilde{a}(s-t_k)-\sigma \tilde{W}_{t_k,s}}\mathrm{d}s.

Return [Y_0, Y_h, Y_{2h}, \ldots] where h = t_max / n_t.

Generate on the fly the brownian motion, no memory overflow but not reproducible.
Ues standardized brownian motion then scale it.
*)
let parabola_igbm a b sigma y0 t_max n_t n_int =
    let tildea = a+.0.5*.sigma*.sigma and h = t_max /. float_of_int n_t in
    let sqrth = sqrt h and ds = h /. float_of_int n_int in
    let grid = range 0. (1./.float_of_int n_int) n_int in
    let res = Array.make (n_t + 1) 0.0 in
    let parabola_fun = parabola_brownian n_int ~grid:grid in
    let aux res =
        res.(0) <- y0;
        for i = 1 to n_t do
            let path = bm_paths 0. 1. 1. n_int 1 |> first (*(W_s)_{s\in[0,1]} will be rescaled*) in
            let w1 = path |> last in
            let parabola = parabola_fun path ~w1:w1 in
            let sum_integrand = fun pre s -> pre +. exp (tildea*.s*.h-.sigma*.sqrth*.(parabola s)) *. ds in
            res.(i) <- exp (-1.*.tildea*.h+.sqrth *.sigma*.w1)*.(res.(i-1)+.a*.b*.(Array.fold_left sum_integrand 0.0 grid));
        done;
    in let _  = aux res in
    res
    ;;



(*
Log-ODE method
Y_{k+1} = Y_k e^{-\tilde{a}h + \sigma W_{t_k,t_{k+1}}} + abh\bigg(1-\sigma H_{t_k,t_{k+1}}+\sigma^2\bigg(3/5h H_{t_k,t_{k+1}}^2+1/30h^2\bigg)\bigg)\frac{e^{-\tilde{a}h+\sigma W_{t_k,t_{k+1}}}-1}{-\tilde{a}h+\sigma W_{t_k,t_{k+1}}}
*)
let log_ode_igbm a b sigma y0 t_max n_t n_int =
    let tildea = a+.0.5*.sigma*.sigma and h = t_max /. float_of_int n_t in
    let sqrth = sqrt h and h2 = h*.h and space_time_levy_area_fun = space_time_levy_area_fun n_int ?grid:None in
    let res = Array.make (n_t + 1) 0.0 in
    let aux res =
        res.(0) <- y0;
        for i = 1 to n_t do
            let path = bm_paths 0. 1. 1. n_int 1 |> first (*(W_s)_{s\in[0,1]} will be rescaled*) in
            let w1 = path |> last in
            let space_time_levy_area = space_time_levy_area_fun path (Some w1) in
            res.(i) <- res.(i-1)*.(exp (-1.*.tildea*.h+.sigma*.sqrth*.w1)) +. a*.b*.h*.(
                        1. -. sigma *. space_time_levy_area *. sqrth +. sigma*.sigma*.(0.6*.h2*.space_time_levy_area*.space_time_levy_area+.1./.30.*.h2))*.
                        (exp (-1.*.tildea*.h+.sigma*.sqrth*.w1) -. 1.) /. (-1. *.tildea*.h+.sigma*.sqrth*.w1);
        done;
        in let _ = aux res in
        res
    ;;

(*
Construct the log-ODE solution given given a path, is reproducible, can cause stack overflow.
Slower than log_ode_igbm given n_int == n_t.
*)
let log_ode_igbm_given_path a b sigma y0 path n_t t_max =
    let tildea = a+.0.5*.sigma*.sigma and n_int = Array.length path / n_t and h = t_max /. float_of_int n_t in
    let sqrth = sqrt h and h2 = h*.h in
    let standardized_brownians = split_and_normalize_brownian path n_t h in
    let space_time_levy_area_fun = space_time_levy_area_fun n_int ?grid:None in
    let res = Array.make (n_t + 1) 0. in
    let rec aux res =
        res.(0) <- y0;
        for i = 1 to n_t do
            let current_path = standardized_brownians.(i-1) in
            let w1 = current_path |> last in
            let space_time_levy_area = space_time_levy_area_fun current_path (Some w1) in
            res.(i) <- res.(i-1)*.(exp (-1.*.tildea*.h+.sigma*.sqrth*.w1)) +. a*.b*.h*.(
                                1. -. sigma *. space_time_levy_area *. sqrth +. sigma*.sigma*.(0.6*.h2*.space_time_levy_area*.space_time_levy_area+.1./.30.*.h2))*.
                                (exp (-1.*.tildea*.h+.sigma*.sqrth*.w1) -. 1.) /. (-1. *.tildea*.h+.sigma*.sqrth*.w1)
        done;
    in let _ = aux res in
    res
    ;;

(*
Construct the parabola-ODE solution given given a path, is reproducible, can cause stack overflow.
Slower than log_ode_igbm given n_int == n_t.
*)
let parabola_igbm_given_path a b sigma y0 path n_t t_max =
    let tildea = a+.0.5*.sigma*.sigma in
    let n_int = Array.length path / n_t and h = t_max /. float_of_int n_t in
    let ds = h /. float_of_int n_int and sqrth = sqrt h in
    let grid = range 0. (1./. float_of_int n_int) n_int in
    let standardized_brownians = split_and_normalize_brownian path n_t h in
    let res = Array.make (n_t + 1) 0. in
    let parabola_fun = parabola_brownian n_int ~grid:grid in
    let rec aux res =
        res.(0) <- y0;
        for i = 1 to n_t do
            let current_path = standardized_brownians.(i-1) in
            let w1 = current_path |> last in
            let parabola = parabola_fun current_path ~w1:w1 in
            let sum_integrand = fun pre s -> pre +. exp (tildea*.s*.h-.sigma*.sqrth*.(parabola s)) *. ds in
            res.(i) <- exp (-1.*.tildea*.h+.sqrth *.sigma*.w1)*.(res.(i-1)+.a*.b*.(Array.fold_left sum_integrand 0.0 grid))
        done;
    in let _ = aux res in
    res
    ;;
