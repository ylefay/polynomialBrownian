type monome = {coeff: float; deg: int};;
type polynome = monome list;;

let rec somme p1 p2 =
    match p1, p2 with
    | [], _ -> p2
    | _, [] -> p1
    | ({deg=n} as m1)::p1p, ({deg=m} as m2)::p2p ->
        if n > m then
            m1::somme p1p p2
        else if n < m then
            m2::somme p1 p2p
        else
            let c = m1.coeff +. m2.coeff in
            if c == 0. then
                somme p1p p2p
            else
                {coeff=c; deg=n}::somme p1p p2p
;;
let rec produit p1 p2 =
        match p1, p2 with
        | [], _ -> []
        | _, [] -> []
        | ({deg=n} as m1)::[], ({deg=m} as m2)::p2p ->
            {coeff=m1.coeff*.m2.coeff; deg=m+n}::produit p1 p2p
        | ({deg=n} as m1)::p1p, p2p ->
            somme (produit [m1] p2p) (produit p1p p2p)
;;

let rec puissance x n =
    match n with
    | 0 -> 1.
    | _ -> x*. puissance x (n-1)

let rec evaluer p x =
    match p with
    | [] -> 0.
    | ({deg=n; coeff=c})::pp -> (c *. puissance x n ) +. evaluer pp x
;;

