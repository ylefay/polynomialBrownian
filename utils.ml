
let last arr =
    arr.(Array.length arr - 1);;
let first arr =
    arr.(0);;
(*
split arr in multiple subarr
*)
let split_f array n =
    let len = Array.length array in
    let l = len / n in
    let res = Array.make n [||] in
    let aux res =
        for i = 0 to n-2 do
            res.(i) <- Array.sub array (i*l) l;
        done;
        res.(n-1) <- Array.sub array ((n-1)*l) (len - (n-1)*l);
    in
    let _ = aux res in
    res;;


(*
cumulative sum function
*)
let c_sum init array =
    let n = Array.length array in
    let res = Array.make n 0.0 in
    let aux res =
        res.(0) <- init +. array.(1);
        for i = 1 to n-1 do
            res.(i) <- res.(i-1) +. array.(i);
        done;
    in
    let _ = aux res in
    res;;

(*
split a brownian motion into n_t parts
and normalize each part:
(W_t) -> ((h^{-1/2}(W_s-W_kh))_{kh<=s<=(k+1)h})
*)
let split_and_normalize_brownian path n_t h =
    let splitted_brownian = split_f path n_t in
    let standardized_brownians =
        splitted_brownian
        |> Array.map (fun w -> Array.map (fun wt -> 1./.(sqrt h) *. (wt -. w.(0))) w)
    in standardized_brownians;;

let pprint x = print_float x; print_string ",";;

(*
range function, maybe use an iterator, iter from c-cube repo ?
*)
let range start step nn =
  let arr = Array.make nn 0.0 in
  let aux arr =
      for i = 0 to nn-1 do
        arr.(i) <- start +. float_of_int i *. step;
      done;
  in
  let _ = aux arr in
  arr;;
