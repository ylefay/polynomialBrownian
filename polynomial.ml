type monome = {coeff: float; deg: int};;
type polynome = monome list;;

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
    | ({coeff=c1; deg=n})::p1p, p2p ->
        (List.map (fun {coeff=c2; deg=m} -> {coeff=c1*.c2; deg=m+n}) p2)
        |> somme (produit p1p p2p)

let rec puissance x n =
    if n = 0 then 1. else x *. puissance x (n-1)

let rec evaluer p x =
    match p with
    | [] -> 0.
    | {deg=n; coeff=c}::pp -> (c *. puissance x n ) +. evaluer pp x
;;

