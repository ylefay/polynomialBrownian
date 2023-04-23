open Array;;
open Utils;;
let pi = Float.pi;;
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
    mu +. sigma *. z;;


(*
Generate normal increments
*)
let dW_gen n dt =
  let dW () = normal_gen_2 0. (sqrt dt) in
  fun () -> init n (fun _ -> dW ());;

(*
let num_domains = try int_of_string Sys.argv.(1) with _ -> 1
*)
let num_domains = 4;;
(*
hyperfine './brownian_parallel.out' './brownian_parallel.out 4' :
...
Summary
    './brownian_parallel.out 4' ran
    1.13 pm 0.04 times faster than './brownian_parallel.out'
*)
(*
let dW_gen_par n dt =
    let dW () = normal_gen_2 0. (sqrt dt) in
        let list_dW size = init size (fun _ -> dW ()) in
        let list_increments = init num_domains (fun _ -> Domain.spawn (fun _ -> list_dW (n/num_domains)))
        and rest_increments = n mod num_domains |> list_dW in
        fun () -> rest_increments@(map Domain.join list_increments |> concat)
        ;;
*)



(* random seed *)
let _ = Random.self_init();;
(* fix seed *)
(*let _ = Random.init 1;;*)

(*
Construct brownian paths
dW_gen_par is not faster than dW_gen?
*)
let bm_paths x0 sigma t n m =
    let dt = t /. (float_of_int n) in
    let generate = dW_gen n (sigma*.sigma*.dt) in
    init m (fun _ ->
        generate ()
        |> c_sum x0
    );;
