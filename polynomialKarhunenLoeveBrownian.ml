open Polynomial;;
open Brownian;;
open List;;
open Utils;;


(*
Numerical approximations for stochastic differential equations, Foster
*)

(* Theorem 4.1.10 (A recurrence relation for constructing Jacobi-like polynomials)
Assuming the product and the sum of two polynomials are done in O(1), we obtain
the sequence of Jacobi polynomials of degree between 2 and n in O(n)

*)
let jacobi deg =
    let rec aux n accu =
        if n >= deg -.1. then
            accu
        else
            match accu with
                | pm1::pm2::_ -> aux (n+.1.) ([produit (somme (produit [{coeff=(n+.1.)*.(2.*.n+.1.); deg=1}] pm1) (produit [{coeff=(-.1.)*.n*.(n+.1.); deg=0}] pm2)) ([{coeff=1./.(n*.(n+.2.)); deg=0}])] @ accu)
                | _ -> accu
    in match deg with
    | 2. -> [[{coeff=0.25; deg=2}; {coeff=(-0.25); deg=0}]]
    | 3. -> [[{coeff=0.5; deg=3}; {coeff=(-0.5); deg=1}]; [{coeff=0.25; deg=2}; {coeff=(-0.25); deg=0}]]
    | _ -> aux 2. [[{coeff=0.5; deg=3}; {coeff=(-0.5); deg=1}]; [{coeff=0.25; deg=2}; {coeff=(-0.25); deg=0}]]
;;

(* Theorem 4.1.0, expression for e_k
O(n)
*)

let eigen deg jacobi_list =
    let aux k pol_jac_kp1 =
        let jac_before_composition = (produit [{coeff=1./.k *. (sqrt (k*.(k+.1.)*.(2.*.k+.1.)) ); deg=0}] pol_jac_kp1) in
        fun t -> evaluer jac_before_composition (2.*.t-.1.)
    in map2 aux (init (length jacobi_list) (fun x -> float_of_int ((length jacobi_list) - x))) jacobi_list
;;

let d_eigen deg jacobi_list =
    let aux k pol_jac_kp1 =
        let d_jac_before_composition =
            (produit [{coeff=1./.k *. (sqrt (k*.(k+.1.)*.(2.*.k+.1.)) ); deg=0}] pol_jac_kp1)
            |> derive in
        fun t -> 2.*. (evaluer d_jac_before_composition (2.*.t-.1.)) (*F = P°(2X-1), F' = 2P'°(2X-1)*)
    in map2 aux (init (length jacobi_list) (fun x -> float_of_int ((length jacobi_list) - x))) jacobi_list



let basis deg n path ?w1 eigen_list =
    let w1 = match w1 with
        | None -> path |> rev |> hd
        | Some w1 -> w1 in
    let dt = 1. /. float_of_int n in
    let grid = range 0. dt 1. n in
    let first_term = map2 (fun wt t ->
    match t with
        | 0. | 1. -> 0. (*div by zero*)
        | _ -> (wt -. w1 *. t) *. dt *. 1. /. (t *. (t -. 1.)))
    path grid in
    List.map (fun eigen_fun -> map2 (fun x y -> x *. y) (map eigen_fun grid) first_term |> fold_left (+.) 0.0) eigen_list
;;

(*
A Polynomial based Karhunen-Loève theorem
See proof of theorem 4.1.6

I_1, ..., I_deg s.t
W = W_1*t+I_1 e_1(t)+...
*)
let compute_basis deg n standard_path =
    let fdeg = float_of_int deg in
    jacobi (fdeg+.1.)
    |> eigen fdeg
    |> basis fdeg n standard_path ?w1:None
    |> iter pprint
;;

let sqrt6 = 2.4494897427831780982;;
let parabola_brownian n path ?(w1) =
    let w1 = match w1 with
        | None -> path |> rev |> hd
        | Some w1 -> w1 in
    let eigen_func = jacobi 2. |> eigen 1. in
    let i1 = eigen_func |> basis 1. n path ~w1:w1 |> hd in
    (* let e1 = eigen_func |> hd in
    fun t -> (w1*.t +. i1*.(e1 t));;*)  (*Wpara(t) = W_1t + I_2sqrt(6)t(t-1) *)
    fun t -> w1*.t +. i1*.sqrt6*.t*.(t-.1.)

(*
Space time levy area, see Definition 4.2.1
*)
let oversqrt6 = 0.40824829046386301637;;
let space_time_levy_area_fun n_int =
    let eigen_func = jacobi 2. |> eigen 1. in
    fun path ?w1 -> begin
        let w1 = match w1 with
                    | None -> path |> rev |> hd
                    | Some w1 -> w1 in
        let i1 = eigen_func |> basis 1. n_int path ~w1:w1 |> hd in
        let space_time_levy_area_value = oversqrt6 *. i1 in
        space_time_levy_area_value
        end;;
(*let path = bm_paths 0. 1. 1. 50 1 |> hd |> iter pprint
    ;;*)
(*map (fun f -> f 0.75) (eigen 1. (jacobi 2.)) |> iter pprint;;*)
(*
compute_basis 10 500;;*)

(*parabola_brownian 500 (bm_paths 0. 1. 1. 500 1 |> hd) 5.;;*)
