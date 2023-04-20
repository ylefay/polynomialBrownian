open List;;
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
        else
            match accu2 with
                | x::q -> aux ([[x]]@accu) q 1
                | _ -> rev accu
    in aux [[]] list 0;;

(*
cumulative sum function
*)
let c_sum init l =
  let rec sum pre acc = function
    | [] -> rev acc
    | hd::tl -> let tmp_sum = pre+.hd in sum tmp_sum (tmp_sum::acc) tl
  in sum init [] l

