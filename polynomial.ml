type monome = {coeff: float; deg: int};;
type polynome = monome list;;
open Printf

let rec somme p1 p2 =
  match p1, p2 with
  | [], _ -> p2
  | _, [] -> p1
  | ({coeff=c1; deg=n} as m1)::p1p, ({coeff=c2; deg=m} as m2)::p2p ->
      match compare n m with
      | 1 -> m1::somme p1p p2
      | -1 -> m2::somme p1 p2p
      | 0 -> if c1 +. c2 = 0.0 then somme p1p p2p else {coeff=c1 +. c2; deg=n}::somme p1p p2p

let rec produit p1 p2 =
    match p1, p2 with
    | [], _ | _, [] -> []
    | ({coeff=c1; deg=n})::p1p, p2 ->
        (List.map (fun {coeff=c2; deg=m} -> {coeff=c1*.c2; deg=m+n}) p2)
        |> somme (produit p1p p2)

(*
Faster than built-in exponentiation
*)
let rec puissance_opti_square x n = match n with
    | 0 -> 1.
    | 1 -> x
    | n -> let b = puissance_opti_square x (n / 2) in b*.b*. (if n mod 2 = 0 then 1. else x)


let rec evaluer p x =
    match p with
    | [] -> 0.
    | {deg=n; coeff=c}::pp -> (c *. puissance_opti_square x n ) +. evaluer pp x
;;


let afficher_monome {coeff=c ; deg=d} =
  if c > 0. then printf "+" ;
  printf "%fX^%d" c d
;;

let afficher p = match p with
| [] -> printf "0" (* un cas particulier, quand même *)
| _  -> List.iter afficher_monome p
;;

(* Afficher avec retour à la ligne *)
let afficher_bis p = afficher p ; print_newline ()
;;