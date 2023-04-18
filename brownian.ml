let pi = Float.pi
open List

let _ = Random.self_init()

(*
cumulative sum function
*)
let c_sum l =
  let rec sum pre acc = function
    | [] -> rev acc
    | hd::tl -> let tmp_sum = pre+.hd in sum tmp_sum (tmp_sum::acc) tl
  in sum 0. [] l

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
    let generate = dW_gen n dt in
    init m (fun _ ->
        generate ()
        |> c_sum
        |> map (fun bm -> x0 +. sigma *. bm)
    );;

(*slower
let bm_paths_bis x0 sigma t n m =
  let dx = sigma*.(t/.(float_of_int n) |> sqrt) in
  let rec loop i x =
    if i = n then []
    else
      let x' = x +. dx *. (normal_gen_2 0.0 1.0) in
      x' :: loop (i+1) x'
  in
  init m (fun _ -> loop 0 x0);;
*)


let bm_paths_bis = bm_paths;;