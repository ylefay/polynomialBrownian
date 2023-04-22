open Brownian;;
open Utils;;
open PolynomialKarhunenLoeveBrownian;;
open List;;

(*
Parabola-ODE method for Stratovitch SDE:
\mathrm{d}y_{t} = f_0(y_t)\mathrm{d}t+f_1(y_t)\circ\mathrm{d}W_t

See theorem 4.3.11, the parabola-ODE method gives
dX_u = f_0(X_u)du + f_1(X_u)d\tilde{W}_u
X_s = y_s

Let W be a standardized brownian motion, let h > 0, and n_int be the number of points we consider inside every intervals [kh,(k+1)h]
The corresponding numerical scheme to obtain Y_{k+1} ~= y_{(k+1}h} from Y_k is

z_0 = Y_k
z_{i+1} = z_i + \sqrt{h}f_1(z_{i}})(\tilde{W}_{i/n_int}-\tilde{W}_{(i-1)/n_int}) + f_0(z_{i})ds
z_{n_int} = Y_{k+1}
where ds = h/n_int.

See the diffusion term is in L^2(P) ~ sqrt(h/n_int) = sqrt(ds), over n_int step, sqrt(h).
Equivalently, we can replace
    (\tilde{W}_{i/n_int}-\tilde{W}_{(i-1)/n_int})~dW=((6-12u)I1/sqrt6 + W1)du
     by setting du = 1/n_int, u = i/n_int
    1/n_int*((6-12i/n_int)I1/sqrt6 + W1)
We still obtain a diffusion term prop. to sqrt(h) after n_int steps.
*)
let parabola_given_path path f0 f1 y0 n_t t_max =
    let n_int = length path / n_t and h = t_max /. float_of_int n_t in
    let ds = h /. float_of_int n_int and sqrth = sqrt h in
    let du = 1./.(float_of_int n_int) and space_time_levy_area_fun = space_time_levy_area_fun n_int in
    let grid = range 0. du 1. n_int in
    let standardized_brownians = split_and_normalize_brownian path n_t h in
    let rec aux accu ongoing_standardized_brownians k =
        if k <= n_t then
            match ongoing_standardized_brownians with
                | current_path::other_paths ->
                let w1 = current_path |> rev |> hd in
                let space_time_levy_area = space_time_levy_area_fun current_path ?w1:(Some w1) in
                                       match accu with
                                        | yk::_-> (*numerical scheme*)
                                            let sum_integrand = fun pre u -> pre +. (sqrth) *. (f1 pre) *. (w1 +. (6. -. 12. *. u)*.space_time_levy_area)*.du +. (f0 pre)*.ds in
                                            let ykp1 = (fold_left sum_integrand yk grid) in
                                        aux (ykp1::accu) other_paths (k+1)
        else
            rev accu
    in aux [y0] standardized_brownians 1
    ;;


(*
Polynomial-ODE method for Stratovitch SDE:
\mathrm{d}y_{t} = f_0(y_t)\mathrm{d}t+f_1(y_t)\circ\mathrm{d}W_t

The Polynomial-ODE method is similar to the Parabola-ODE method where
we replace \tilde{W} by W^n, the n-th degree polynomial approximation of W.
Or equivalently, we replace \circ\mathrm{d}W_t by \mathrm{d}W^n_t.

We have ||W-W^n||_{L_2(P)} = O(n^{-1/2}).

Here we use the derivation of polynomials to compute dW^n/du then multiply it by du.
*)
let polynomial_given_path deg path f0 f1 y0 n_t t_max =
    let n_int = length path / n_t and h = t_max /. float_of_int n_t in
    let ds = h /. float_of_int n_int and sqrth = sqrt h in
    let du = 1./.(float_of_int n_int) in
    let jac = jacobi (float_of_int deg) in
    let deigen_list =
        jac
        |> d_eigen (float_of_int deg -. 1.)
        and
        eigen_list =
        jac
        |> eigen (float_of_int deg -. 1.) in (*de_k/dX, e_k*)
    let grid = range 0. du 1. n_int in
    let standardized_brownians = split_and_normalize_brownian path n_t h in
    let rec aux accu ongoing_standardized_brownians k =
        if k <= n_t then
            match ongoing_standardized_brownians with
                | current_path::other_paths ->
                let w1 = current_path |> rev |> hd in
                let basis_coefficients = basis (float_of_int deg -. 1.) n_int current_path ~w1:(w1) eigen_list in
                let dbrownian_pol_fun = fun t ->
                    map2 (fun coeff deigen_fun -> coeff *. (deigen_fun t)) basis_coefficients deigen_list
                    |> fold_left (+.) (w1) in
                                       match accu with
                                        | yk::_-> (*numerical scheme*)
                                            let sum_integrand = fun pre u -> pre +. (sqrth) *. (f1 pre) *. (dbrownian_pol_fun u) *. du+. (f0 pre)*.ds in
                                            let ykp1 = (fold_left sum_integrand yk grid) in
                                        aux (ykp1::accu) other_paths (k+1)
        else
            rev accu
    in aux [y0] standardized_brownians 1
    ;;

