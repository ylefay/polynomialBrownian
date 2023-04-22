open List;;
open Utils;;
open Normal;;


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

