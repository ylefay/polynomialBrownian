open PolynomialKarhunenLoeveBrownian;;
open Utils;;

(*
Lie bracket
Let U, V : R -> R
[U, V](f) := (UV - VU)(f) = (U(x)∂_x(V)(x) - V(x)∂_x(U)(x))∂_x(f)
*)
let lie dx u v =
  let overdx = 1. /. dx in
  let du = fun x -> (u (x +. dx) -. u x) *. overdx in
  let dv = fun x -> (v (x +. dx) -. v x) *. overdx in
  let duv = fun x -> (u x *. dv x) -. (v x *. du x) in
  fun f x -> duv x *. ((f (x +. dx) -. f x) *. overdx);;

(*
Log-ODE method for Stratovitch SDE:
\mathrm{d}y_{t} = f_0(y_t)\mathrm{d}t+f_1(y_t)\circ\mathrm{d}W_t

See theorem 4.3.8, the log-ODE method gives
Y_{k+1} to be the solution at u=1 of the following ODE:
\mathrm{d}z/\mathrm{d}u = f_0(z)h+f_1(z)W_{t_k,t_{k+1}} + [f_1,f_0](z)hH_{t_k,t_{k+1}} + [f_1,[f_1,f_0]](z)(3/5hH_{t_k,t_{k+1}}^2+1/30h^2)
z_0 = Y_k

Let W be a standardized brownian motion, let h > 0, and n_int be the number of points we consider inside every intervals [kh,(k+1)h]
The corresponding numerical scheme to obtain Y_{k+1} ~= y_{(k+1}h} from Y_k is

z_0 = Y_k
z_{i+1} = z_i + f_0(z_{i})ds + \sqrt{h}(f_1(z_{i}) W1+ + [f_1,f_0](z)h I_1/sqrt(6)
    + [f_1,[f_1,f_0]](z)(1/10hI1^2+1/30 h^1.5))
z_{n_int} = Y_{k+1}
where ds = h/n_int.
*)
let log_given_path path f0 f1 y0 n_t t_max =
    let n_int = Array.length path / n_t and h = t_max /. float_of_int n_t in
    let ds = h /. float_of_int n_int and sqrth = sqrt h in
    let h2 = h*.h in
    let du = 1./.(float_of_int n_int) and space_time_levy_area_fun = space_time_levy_area_fun n_int in
    let grid = range 0. du n_int in
    let standardized_brownians = split_and_normalize_brownian path n_t h in
    let lie_bracket_10 = lie du f1 f0 (fun z -> z) in let lie_bracket_110 = lie du f1 lie_bracket_10 (fun z -> z) in
    let res = Array.make (n_t + 1) 0.0 in
    let rec aux res =
        res.(0) <- y0;
        for i = 1 to n_t do
            let current_path = standardized_brownians.(i-1) in
            let w1 = current_path |> last in
            let space_time_levy_area = space_time_levy_area_fun current_path (Some w1) in
            (* numerical scheme *)
            let sum_integrand = fun pre u ->
                pre +. (f0 pre)*.ds +.
                (sqrth) *.
                (
                (f1 pre) *. w1 +.
                (lie_bracket_10 pre) *. h *. space_time_levy_area +.
                (lie_bracket_110 pre) *. (0.6*.h*.space_time_levy_area*.space_time_levy_area +. 1./.30. *. h2)
                )*.du
            in
            res.(i) <- Array.fold_left sum_integrand res.(i-1) grid
        done;
        in let _ = aux res in
        res
    ;;