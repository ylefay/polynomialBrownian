let pi = Float.pi
open List

(* random seed *)
let _ = Random.self_init();;
(* fix seed *)
(*let _ = Random.init 1;;*)


(*
cumulative sum function
*)
let c_sum init l =
  let rec sum pre acc = function
    | [] -> rev acc
    | hd::tl -> let tmp_sum = pre+.hd in sum tmp_sum (tmp_sum::acc) tl
  in sum init [] l

(*
range function, maybe use an iterator, iter from c-cube repo ?
*)
let range start step stop nn =
    let rec aux j n =
            if j < stop && n < nn then j :: aux (j+.step) (n+1) else []
    in aux start 0;;

(*
Box-Muller method
*)
let normal_gen mu sigma =
    let u = sqrt (-2. *. log (Random.float 1.)) *. cos (2. *. pi *. Random.float 1.)
    in mu +. sigma *. u;;

(*
Faster but imprecise
Rao et al.
*)
let normal_gen_2 mu sigma =
    let u = Random.float 1. in
    let z = -. log (1. /. u -. 1.) /. 1.702 in
    mu +. sigma *. z


(*
Generate normal increments
*)
let dW_gen n dt =
  let dW () = normal_gen_2 0. (sqrt dt) in
  fun () -> init n (fun _ -> dW ())

(*
Construct brownian paths
*)
let bm_paths x0 sigma t n m =
    let dt = t /. (float_of_int n) in
    let generate = dW_gen n (sigma*.sigma*.dt) in
    init m (fun _ ->
        generate ()
        |> c_sum x0
    );;

(*
split list in multiple sublists
*)
let split_f list n =
    let l = length list / n in
    let rec aux accu accu2 j =
        if j < l then
            match accu2 with
                | x::q -> match accu with
                    | y::qp -> aux ([y@[x]]@qp) q (j+1)
                    | _ -> aux [[x]] q (j+1)
                | [] -> rev accu
        else
            match accu2 with
                | x::q -> aux ([[x]]@accu) q 1
                | _ -> rev accu
    in aux [[]] list 0;;

