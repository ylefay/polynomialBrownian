open List;;
open Brownian;;
open Igbm;;
open Utils;;
open Parabola_method;;

let n = try int_of_string Sys.argv.(1) with _ -> 1000000;;

(*on the fly, not reproducible*)
(*
Printf.printf "\n parabola_ode \n";;
parabola_igbm 0.1 0.04 0.6 0.06 1. 1000 1000 |> iter pprint;;

Printf.printf "\n log_ode \n";;
log_ode_igbm 0.1 0.04 0.6 0.06 1. 1000 1000 |> iter pprint;;
*)

(*given path*)
let mypath = bm_paths 0.0 1.0 1. n 1 |> hd;; (*n_int*n_t*)
Printf.printf "\n parabola_ode_given_path \n";;
parabola_igbm_given_path 0.1 0.04 0.6 0.06 mypath 1000 1. |> iter pprint;;

Printf.printf "\n log_ode_given_path \n";;
log_ode_igbm_given_path 0.1 0.04 0.6 0.06 mypath 1000 1. |> iter pprint;;

(* defining IGBM stratovitch sde*)
let a = 0.1 and b = 0.04 and sigma=0.6 and y0 = 0.06;;
let atilde = a +. 0.5*.sigma*.sigma and btilde = 2.*.a*.b/.(2.*.a+.sigma*.sigma);;
let f0 yt = atilde*.(btilde-.yt) and f1 yt = sigma *. yt;;

Printf.printf "\n general_parabola_ode_given_path \n";;
parabola_given_path mypath f0 f1 y0 1000 1. |> iter pprint;;
Printf.printf "\n general_deg3_ode_given_path \n";;
polynomial_given_path 3 mypath f0 f1 y0 1000 1. |> iter pprint;;
Printf.printf "\n general_deg2_ode_given_path \n";;
polynomial_given_path 3 mypath f0 f1 y0 1000 1. |> iter pprint;;
Printf.printf "\n general_deg1_ode_given_path \n";;
polynomial_given_path 3 mypath f0 f1 y0 1000 1. |> iter pprint;;