let pi = Float.pi
open Random
open List

let c_sum l =
  let rec sum pre acc = function
    | [] -> rev acc
    | hd::tl -> let tmp_sum = pre+.hd in sum tmp_sum (tmp_sum::acc) tl
  in sum 0. [] l

let sum l = fold_left (fun x y -> x +. y) 0. l

let range start step stop nn =
    let rec aux j n =
            if j < stop && n < nn then j :: aux (j+.step) (n+1) else []
    in aux start 0;;

(*
Box-Muller method
*)
let normal_gen mu sigma = let u = sqrt (-2. *. log (Random.float 1.)) *. cos (2. *. pi *. Random.float 1.)
    in mu +. sigma *. u;;

let dW_gen n dt =
  let dW () = normal_gen 0. (sqrt dt) in
  fun () -> init n (fun _ -> dW ())

let bm_paths x0 sigma t n m =
    let dt = t /. (float_of_int n) in
    let generate = dW_gen n dt in
    init m (fun _ ->
        generate ()
        |> c_sum
        |> map (fun bm -> x0 +. sigma *. bm)
    )
