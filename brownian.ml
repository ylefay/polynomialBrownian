open List;;
open Utils;;
open Normal;;


(* random seed *)
let _ = Random.self_init();;
(* fix seed *)
(*let _ = Random.init 1;;*)

(*
range function, maybe use an iterator, iter from c-cube repo ?
*)
let range start step stop nn =
    let rec aux j n =
            if j < stop && n < nn then j :: aux (j+.step) (n+1) else []
    in aux start 0;;

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

