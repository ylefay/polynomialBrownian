open Polynomial;;
open Float;;
(* Theorem 4.1.10 (A recurrence relation for constructing Jacobi-like polynomials
P_2(x) = 1/4(x-1)(x+1)
P_3(x) = 1/2x(x-1)(x+1)
n(n+2)P_(n+2)(x) = (n+1)(2n+1)xP_(n+1)(x) - n(n+1)P_n(x)
*)
let rec jacobi deg =
    match deg with
    | 2. -> [{coeff=0.25; deg=2}; {coeff=(-0.25); deg=0}]
    | 3. -> [{coeff=0.5; deg=3}; {coeff=(-0.5); deg=1}]
    | _ -> let n = deg -. 2. in  produit ([{coeff=1./.(n*.(n+.2.)); deg=0}]) (somme (produit [{coeff=(n+.1.)*.(2.*.n+.1.); deg=1}] (jacobi (n+.1.))) (produit [{coeff=n*.(n+.1.); deg=0}] (jacobi n)))
;;

(* Theorem 4.1.0, expression for e_k
e_k(t) = 1/k sqrt(k(k+1)(2k+1)) P_(k+1)(2t-1)
*)
let ek k = let p = (produit [{coeff=1./.k *. (sqrt (k*.(k+.1.)*.(2.*.k+.1.)) ); deg=0}] (jacobi (k+.1.))) in
    function t -> evaluer p (2.*.t-.1.);;