open Polynomial;;
open Rand;;
open Utils;;
open Array;;


(*
Numerical approximations for stochastic differential equations, Foster 2020.
*)

(* Theorem 4.1.10 (A recurrence relation for constructing Jacobi-like polynomials)
Assuming the product and the sum of two polynomials are done in O(1), we obtain
the sequence of Jacobi polynomials of degree between 2 and n in O(n)
*)
let jacobi deg =
      let res = Array.make (deg - 1) [] in
      let aux res =
          match deg with
          | 1 -> ();
          | 2 -> res.(0) <- [{coeff=0.25; deg=2}; {coeff=(-0.25); deg=0}];
          | 3 -> res.(0) <- [{coeff=0.25; deg=2}; {coeff=(-0.25); deg=0}]; res.(1) <- [{coeff=0.5; deg=3}; {coeff=(-0.5); deg=1}];
          | _ -> begin
                    res.(0) <- [{coeff=0.25; deg=2}; {coeff=(-0.25); deg=0}];
                    res.(1) <- [{coeff=0.5; deg=3}; {coeff=(-0.5); deg=1}];
                    for i = 2 to (deg - 2) do
                        let n = float_of_int i in
                        let pm1 = res.(i-1) in
                        let pm2 = res.(i-2) in
                        let p1 = produit [{coeff=(n+.1.)*.(2.*.n+.1.); deg=1}] pm1 in
                        let p2 = produit [{coeff=(-.1.)*.n*.(n+.1.); deg=0}] pm2 in
                        let p3 = somme p1 p2 in
                        let p4 = produit [{coeff=1./.(n*.(n+.2.)); deg=0}] p3 in
                        res.(i) <- p4;
                    done
                end
      in
      let _ = aux res in res;;

(* Theorem 4.1.0, expression for e_k
O(n)
*)
let eigen deg jacobi_list =
    let aux k pol_jac_kp1 =
        let jac_before_composition = (produit [{coeff=1./.k *. (sqrt (k*.(k+.1.)*.(2.*.k+.1.)) ); deg=0}] pol_jac_kp1) in
        fun t -> evaluer_horner_jacobi jac_before_composition (2.*.t-.1.)
    in map2 aux (init (length jacobi_list) (fun x -> float_of_int ((length jacobi_list) - x))) jacobi_list
;;

let d_eigen deg jacobi_list =
    let aux k pol_jac_kp1 =
        let d_jac_before_composition =
            (produit [{coeff=1./.k *. (sqrt (k*.(k+.1.)*.(2.*.k+.1.)) ); deg=0}] pol_jac_kp1)
            |> derive in
        fun t -> 2.*. (evaluer_horner_jacobi d_jac_before_composition (2.*.t-.1.)) (*F = P°(2X-1), F' = 2P'°(2X-1)*)
    in map2 aux (init (length jacobi_list) (fun x -> float_of_int ((length jacobi_list) - x))) jacobi_list



let basis deg n ?grid eigen_list path w1 =
    let w1 = match w1 with
        | None -> path |> last
        | Some w1 -> w1 in
    let dt = 1. /. float_of_int n in
    let grid = match grid with
        | None ->  range 0. dt n
        | Some grid -> grid in
    let first_term = map2 (fun wt t ->
    match t with
        | 0. | 1. -> 0. (*div by zero*)
        | _ -> (wt -. w1 *. t) *. dt *. 1. /. (t *. (t -. 1.)))
    path grid in
    Array.map (fun eigen_fun -> map2 (fun x y -> x *. y) (map eigen_fun grid) first_term |> fold_left (+.) 0.0) eigen_list
;;

(*
A Polynomial based Karhunen-Loève theorem
See proof of theorem 4.1.6

I_1, ..., I_deg s.t
W = W_1*t+I_1 e_1(t)+...
*)
let compute_basis deg n =
    let fdeg = float_of_int deg in
    let fun_compute_basis =
    jacobi (deg+1)
    |> eigen fdeg
    |> basis fdeg n ?grid:None
    in
    fun standard_path -> fun_compute_basis standard_path None;;


let eigen_func = jacobi 2 |> eigen 1.
let sqrt6 = 2.4494897427831780982;;
let parabola_brownian n ?grid path ?w1 =
    let w1 = match w1 with
        | None -> path |> last
        | Some w1 -> w1 in
    let grid = match grid with
        | None ->  range 0. (1./.float_of_int n) n
        | Some grid -> grid in
    let i1 = basis 1. n ~grid:grid eigen_func path (Some w1) |> first in
    (* let e1 = eigen_func |> hd in
    fun t -> (w1*.t +. i1*.(e1 t));;*)  (*Wpara(t) = W_1t + I_2sqrt(6)t(t-1) *)
    fun t -> w1*.t +. i1*.sqrt6*.t*.(t-.1.)
;;

(*
Space time levy area, see Definition 4.2.1
*)
let oversqrt6 = 0.40824829046386301637;;
let space_time_levy_area_fun n_int ?grid =
    let grid = match grid with
        | None ->  range 0. (1./.float_of_int n_int) n_int
        | Some grid -> grid in
    let basis_coeff_fun = eigen_func |> basis 1. n_int ~grid:grid in
    fun path w1 -> begin
        let w1 = match w1 with
                    | None -> path |> last
                    | Some w1 -> w1 in
        let i1 = basis_coeff_fun path (Some w1) |> first in
        let space_time_levy_area_value = oversqrt6 *. i1 in
        space_time_levy_area_value
        end;;
