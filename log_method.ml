open Parabola_method;;
open PolynomialKarhunenLoeveBrownian;;

(*
Lie bracket
Let U, V : R -> R
[U, V](x) := (UV - VU)(x) = U(x)∂_x(V)(x) - V(x)∂_x(U)(x)
*)

(*
Log-ODE method for Stratovitch SDE:
\mathrm{d}y_{t} = f_0(y_t)\mathrm{d}t+f_1(y_t)\circ\mathrm{d}W_t

See theorem 4.3.8, the log-ODE method gives
Y_{k+1} to be the solution at u=1 of the following ODE:
\mathrm{d}z/\mathrm{d}u = f_0(z)h+f_1(z)W_{t_k,t_{k+1}} + [f_1,f_0](z)hH_{t_k,t_{k+1}} + [f_1,[f_1,f_0]](z)(3/5hH_{t_k,t_{k+1}}^2+1/30h^2)
z_0 = Y_k

Let W be a standardized brownian motion, let h > 0, and n_int be the number of points we consider inside every intervals [kh,(k+1)h]
The corresponding numerical scheme to obtain Y_{k+1} ~= y_{(k+1}h} from Y_k is

z_0 = Y_k
z_{i+1} = z_i + f_0(z_{i})ds + \sqrt{h}f_1(z_{i}})(\tilde{W}_{i/n_int}-\tilde{W}_{(i-1)/n_int}) + [f_1,f_0](z)h^{3/2}I_1/sqrt(6)
    + [f_1,[f_1,f_0]](z)(1/10h^{3/2}I1^2+1/30 h^2)
z_{n_int} = Y_{k+1}
where ds = h/n_int.
*)
