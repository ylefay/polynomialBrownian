let pi = 3.14159265358979312
open Random
open List
open Seq
(*
Box-Muller method
*)
let normal_gen mu sigma = let u = sqrt (-2. *. log (Random.float 1.)) *. cos (2. *. pi *. Random.float 1.)
    in mu +. sigma *. u;;

let dW_gen ~n ~dt =
  let dW () = normal_gen 0. (sqrt dt) in
  fun () -> init n (fun _ -> dW ())

let bm_paths x0 sigma t n m =
    let dt = t /. (float_of_int n) in
    let generate = dW_gen ~n ~dt in
    init m (fun _ ->
        generate ()
        |> map (fun dW -> sigma *. dW)
        |> scan (+.) 0.
        |> mapi (fun i x -> float_of_int i *. dt, x0 +. x)
    )

