open Array;;
open Rand;;
open Igbm;;
open Gbm;;
open Utils;;
open Parabola_method;;
open Log_method;;

let n = try int_of_string Sys.argv.(1) with _ -> 1000000;;
let n_int = try int_of_string Sys.argv.(2) with _ -> 1000;;
(*dt = t_max/n, h = n_int*dt*)

(*
IGBM simulations
*)

(*on the fly, not reproducible*)
(*
Printf.printf "\n parabola_ode \n";;
parabola_igbm 0.1 0.04 0.6 0.06 1. 1000 1000 |> iter pprint;;

Printf.printf "\n log_ode \n";;
log_ode_igbm 0.1 0.04 0.6 0.06 1. 1000 1000 |> iter pprint;;
*)

(*given path*)
let mypath = bm_paths 0.0 1.0 1. n 1 |> first;; (*n_int*n_t*)
Printf.printf "\n BROWNIAN\n";;
mypath |> iter pprint;;

Printf.printf "\n IGBM_Parabola\n";;
parabola_igbm_given_path 0.1 0.04 0.6 0.06 mypath n_int 1. |> iter pprint;;
Printf.printf "\n IGBM_Log_ODE\n";;
log_ode_igbm_given_path 0.1 0.04 0.6 0.06 mypath n_int 1. |> iter pprint;;

(* defining IGBM stratovitch sde*)
let a = 0.1 and b = 0.04 and sigma=0.6 and y0 = 0.06;;
let atilde = a +. 0.5*.sigma*.sigma and btilde = 2.*.a*.b/.(2.*.a+.sigma*.sigma);;
let f0 yt = atilde*.(btilde-.yt) and f1 yt = sigma *. yt;;

Printf.printf "\n polyparabola_integration_IGBM\n";;
parabola_given_path mypath f0 f1 y0 n_int 1. |> iter pprint;;
Printf.printf "\n poly1_integration_IGBM\n";;
polynomial_given_path 1 mypath f0 f1 y0 n_int 1. |> iter pprint;;
Printf.printf "\n poly2_integration_IGBM\n";;
polynomial_given_path 2 mypath f0 f1 y0 n_int 1. |> iter pprint;;
Printf.printf "\n poly3_integration_IGBM\n";;
polynomial_given_path 3 mypath f0 f1 y0 n_int 1. |> iter pprint;;

(*
GBM simulations

let mypath = bm_paths 0.0 1.0 252. n 1 |> hd;;
let mu = 0.04/.252. and sigma=0.22/.15.87 and s0 = 100.;;
Printf.printf "\n closed_form_GBM \n";;
gbm_given_path mu sigma s0 mypath 252. |> iter pprint;;
Printf.printf "\n general_deg2_ode_given_path \n";;
parabola_gbm_given_path mu sigma s0 mypath 252 252. |> iter pprint;;
*)

(*
Different SDEs
Formula to convert from Ito drift to Strat. drift:
a_strat = a_ito - 1/2b∂_x[b]
*)

(*
Case 1:
dy_t = y_t^2 dW_t, y_0 = 1
=> y_t = 1/(1-W_t)
*)
Printf.printf "\n poly2_integration_CASE1\n";;
polynomial_given_path 2 mypath (fun x -> 0.) (fun x -> x*.x) 1. n_int 1. |> iter pprint;;
Printf.printf "\n log_integration_CASE1\n";;
log_given_path mypath (fun x -> 0.) (fun x -> x*.x) 1. n_int 1. |> iter pprint;;
(*
Case 2:
dy_t = y_t dt + y_t dW_t, y_0 = 1
a_strat = 1/2y_t
=> y_t = exp(1/2t + W_t)
*)
Printf.printf "\n poly2_integration_CASE2\n";;
polynomial_given_path 2 mypath (fun x -> 0.5*.x) (fun x -> x) 2. n_int 1. |> iter pprint;;
Printf.printf "\n log_integration_CASE2\n";;
log_given_path mypath (fun x -> 0.5*.x) (fun x -> x) 2. n_int 1. |> iter pprint;;

Printf.printf "\n poly5_integration_CASE2\n";;
polynomial_given_path 5 mypath (fun x -> 0.5*.x) (fun x -> x) 2. n_int 1. |> iter pprint;;
(*
Case 3:
dy_t = sin(y_t)dW_t, y_0 = pi/2
*)
Printf.printf "\n poly2_integration_CASE3\n";;
polynomial_given_path 2 mypath (fun x -> 0.) (fun x -> sin x) (Float.pi /. 2.) n_int 1. |> iter pprint;;
Printf.printf "\n log_integration_CASE3\n";;
log_given_path mypath (fun x -> 0.) (fun x -> sin x) (Float.pi /. 2.) n_int 1. |> iter pprint;;

(*
Case 4:
Cox–Ingersoll–Ross model
dy_t = a(b-y_t)dt+sigma sqrt(y_t)dW_t
<=>
d_yt = (a(b-y_t)-sigma/4)dt + sigma sqrt(y_t)°dW_t
*)
let a = 1. and b = 2. and sigma = 0.2;;
Printf.printf "\n poly2_integration_CASE4\n";;
polynomial_given_path 2 mypath (fun x -> a*.(b-.x)-.sigma*.0.25) (fun x -> sigma*. (sqrt x)) 1. n_int 1. |> iter pprint;;
Printf.printf "\n log_integration_CASE4\n";;
log_given_path mypath (fun x -> a*.(b-.x)-.sigma*.0.25) (fun x -> sigma*. (sqrt x)) 1. n_int 1. |> iter pprint;;

(*
Case 5:
A modified Wright-Fishcer Model
dy_t = y_t(1-y_t)dt + sqrt(y_t(1-y_t))dW_t, y_0 € ]0,1[, y_0 = 0.3
<=>
dy_t = y_t(1-y_t) - (1-2y_t)/(4sqrt(y_t(1-y_t))) + sqrt(y_t(1-y_t))°dW_t
*)
Printf.printf "\n poly2_integration_CASE5\n";;
polynomial_given_path 2 mypath
    (fun x -> x*.(1.-.x)-.0.5*.(1.-.2.*.x)/.(4.*. (sqrt (x*.(1.-.x)))))
    (fun x -> sqrt (x*.(1.-.x)))
0.3 n_int 1. |> iter pprint;;
Printf.printf "\n log_integration_CASE5\n";;
log_given_path mypath
    (fun x -> x*.(1.-.x)-.0.5*.(1.-.2.*.x)/.(4.*. (sqrt (x*.(1.-.x)))))
    (fun x -> sqrt (x*.(1.-.x)))
0.3 n_int 1. |> iter pprint;;


