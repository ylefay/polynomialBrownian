open Polynomial;;
open Float;;
open List;;
open Printf;;


(* Theorem 4.1.10 (A recurrence relation for constructing Jacobi-like polynomials)
P_2(x) = 1/4(x-1)(x+1)
P_3(x) = 1/2x(x-1)(x+1)
n(n+2)P_(n+2)(x) = (n+1)(2n+1)xP_(n+1)(x) - n(n+1)P_n(x)
*)

(*let rec jacobi deg =
    match deg with
    | 2. -> [{coeff=0.25; deg=2}; {coeff=(-0.25); deg=0}]
    | 3. -> [{coeff=0.5; deg=3}; {coeff=(-0.5); deg=1}]
    | _ -> let n = deg -. 2. in  produit ([{coeff=1./.(n*.(n+.2.)); deg=0}]) (somme (produit [{coeff=(n+.1.)*.(2.*.n+.1.); deg=1}] (jacobi (n+.1.))) (produit [{coeff=n*.(n+.1.); deg=0}] (jacobi n)))
;;*)

let jacobi deg =
    let rec aux n accu =
        if n = deg -.2. then
            accu
        else
            match accu with
                | [pm1; pm2; _] -> aux (n+.1.) ([produit ([{coeff=1./.(n*.(n+.2.)); deg=0}]) (somme (produit [{coeff=(n+.1.)*.(2.*.n+.1.); deg=1}] pm1) (produit [{coeff=n*.(n+.1.); deg=0}] pm2))] @ accu)
                | _ -> accu
    in aux 1. [[{coeff=0.5; deg=3}; {coeff=(-0.5); deg=1}]; [{coeff=0.25; deg=2}; {coeff=(-0.25); deg=0}]];;

(* Theorem 4.1.0, expression for e_k
e_k(t) = 1/k sqrt(k(k+1)(2k+1)) P_(k+1)(2t-1)
*)

let e deg jacobi_list =
    let aux k pol_jac_kp1 =
        let p = (produit [{coeff=1./.k *. 1./.(sqrt (k*.(k+.1.)*.(2.*.k+.1.)) ); deg=0}] pol_jac_kp1) in
        function t -> evaluer p (2.*.t-.1.)
    in map2 aux (init (length jacobi_list) (function x -> float_of_int ((length jacobi_list) - x))) jacobi_list
;;

(*
print_float (evaluer (List.hd (jacobi 3.)) 4.)
;;
print_float (aux 2. (List.hd (jacobi 3.)) 4.);*)