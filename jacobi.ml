open Polynomial
open Brownian
open Float
open List
open Printf

(*
Jacobi-like polynomials & eigen functions
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
    | 2. -> [[{coeff=0.5; deg=3}; {coeff=(-0.5); deg=1}]]
    | 3. -> [[{coeff=0.5; deg=3}; {coeff=(-0.5); deg=1}]; [{coeff=0.25; deg=2}; {coeff=(-0.25); deg=0}]]
    | _ -> aux 2. [[{coeff=0.5; deg=3}; {coeff=(-0.5); deg=1}]; [{coeff=0.25; deg=2}; {coeff=(-0.25); deg=0}]]
;;

(* Theorem 4.1.0, expression for e_k
O(n)
*)

let eigen deg jacobi_list =
    let aux k pol_jac_kp1 =
        let p = (produit [{coeff=1./.k *. (sqrt (k*.(k+.1.)*.(2.*.k+.1.)) ); deg=0}] pol_jac_kp1) in
        fun t -> evaluer p (2.*.t-.1.)
    in map2 aux (init (length jacobi_list) (fun x -> float_of_int ((length jacobi_list) - x))) jacobi_list
;;

let pprint x = print_float x; print_string ",";;

(*print_float (evaluer (List.hd (jacobi 3.)) 4.);;*)
(*print_float (aux 2. (List.hd (jacobi 3.)) 4.);;*)
(*print_float (List.hd (e 2. (jacobi 3.)) 4.);;*)
(*let a = map (List.hd (e 2. (jacobi 3.))) (init 50000 (function x -> float_of_int (x+1)))
let () = List.iter print_float a;;*)

(*let paths = bm_paths 0. 1. 1. 1000 1;;
map (fun path -> List.iter print_float path) paths;;*)

let basis deg n path eigen_list =
    let dt = 1. /. float_of_int n in
    let grid = range 0. dt 1. n in
    let w1 = path |> rev |> hd in
    let first_term = map2 (fun x y -> x -. w1 *. y) path grid in
    let aux x y =
        match x with
        | 0. | 1. -> 0.
        | _ -> y *. dt *. 1. /. (x *. (x -. 1.))
    in
    let second_term = map2 aux grid first_term in
    List.map (fun eigen_fun -> List.map2 (fun x y -> x *. y) (List.map eigen_fun grid) second_term |> fold_left (+.) 0.0) eigen_list
;;


let compute_basis deg n =
    let fdeg = float_of_int deg in
    let path = bm_paths 0. 1. 1. n 1 |> hd in
    iter pprint path;
    jacobi (fdeg+.1.)
    |> eigen fdeg
    |> basis fdeg n path
    |> iter pprint
;;


(*let path = bm_paths 0. 1. 1. 50 1 |> hd |> iter pprint
    ;;*)
(*map (fun f -> f 0.75) (eigen 1. (jacobi 2.)) |> iter pprint;;*)
compute_basis 2 500;;
(*let evalp t f = evaluer f t |> pprint;;*)
(*
(jacobi 4.) |> iter (evalp 10.);;
*)