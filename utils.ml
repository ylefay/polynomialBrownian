open List;;
(*
split list in multiple sublists
*)
let split_f list n =
    let l = length list / n in
    let rec aux accu c_list j =
        if j < l then
            begin
            match c_list with
                | x::q ->
                    begin
                    match accu with
                    | y::qp -> aux ([y@[x]]@qp) q (j+1)
                    | [] -> aux [[x]] q (j+1)
                    end
                | [] -> rev accu
            end
        else
            match c_list with
                | x::q -> aux ([[x]]@accu) q 1
                | [] -> rev accu
    in aux [[]] list 0;;

(*
cumulative sum function
*)
let c_sum init l =
  let rec sum pre acc = function
    | [] -> rev acc
    | hd::tl -> let tmp_sum = pre+.hd in sum tmp_sum (tmp_sum::acc) tl
  in sum init [] l

(*
split a brownian motion into n_t parts
and normalize each part:
(W_t) -> ((h^{-1/2}(W_s-W_kh))_{kh<=s<=(k+1)h})
*)
let split_and_normalize_brownian path n_t h =
    let splitted_brownian = split_f path n_t in
    let standardized_brownians = map (fun w -> let w0 = hd w in map (fun wt -> 1./.(sqrt h) *. (wt -. w0)) w) splitted_brownian
    in standardized_brownians;;

let pprint x = print_float x; print_string ",";;
