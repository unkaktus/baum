(* tensorrules.m *)
(* Bernd Bruegmann, 2/96, 10/02 *)

(* Mathematica script, to be read by TensorEquationsToC.m:

   Define rules for tensor formulas like  
   - covariant derivatives
   - determinants
   - factoring out common factors like the metric

   Several entries here are obsolete, make current.
*)  


(* keep track of time *)
timerevents = {}
timerold = 0
timer[a_] := Module[{},
  timernew = TimeUsed[];
  timerevents = Append[timerevents, {a, timernew, timernew-timerold}];
  timerold = timernew;
]
prtimer[] := 
  Print["Timer:  ", PaddedForm[TableForm[timerevents], {8, 2}]]

timer["starting tensorrules.m"]





(**************************************************************************)
(* special tensor knowledge *)

(* hidden symmetries for contraction! *)
gamma[c_,a_,b_] := gamma[c,b,a] /; !OrderedQ[{a,b}]
gammado[c_,a_,b_] := gammado[c,b,a] /; !OrderedQ[{a,b}]

(* compile-time substitutions *)
collectvar = {
    codel[a_, beta[b_]] -> codelbeta[a,b],
    codel[a_, K[b_,c_]] -> codelK[a,b,c],
    ginv[c_,d_] K[a_,c_] K[b_,d_] -> KK[a,b]
}

(* covariant derivatives *)
expcodel = {
  codel[a_, codel[b_,alpha]] -> 
    deldel[a, b, alpha] - gamma[c,a,b] del[c, alpha],

  codel[a_, beta[b_]] -> 
    del[a, beta[b]] + gamma[b,a,c] beta[c],

  codel[a_, K[b_,c_]] -> 
    del[a, K[b,c]] - gamma[d,a,b] K[d,c] - gamma[d,a,c] K[b,d]
}



(************************************************************************)
(* factorizations, mostly obsolete, make current *)
facterms = {
  KK[a_,b_] == x_ :> 
   KK[a,b] == Collect[x, {K[1,a], K[2,a], K[3,a]}],
  (* works: K[a_,b_] x_ + K[a_,b_] y_ -> K[a,b] (x + y), *)

  R[a_,b_] == x_ :> 
    R[a,b] == Collect[x, Union[Flatten[Array[ginv,{3,3}]]]], 

  delt[K[a_,b_]] == x_ :>
    delt[K[a,b]] == Collect[x, {alpha, beta[1], beta[2], beta[3]}]
}  
  


(* miscellaneous *)

(* this saves 11 multiplies or so in 3d but is cumbersome *)
(*
  gg = Array[g, {3,3}];
  gginvdet = Inverse[gg] Det[gg];
  ginvdet[a_Integer,b_Integer] := gginvdet[[a,b]]
  gginv[a,b] == ginvdet[a,b],
  detg == g[1,c] gginv[1,c],
  ginv[a,b] == gginv[a,b] / detg,
*)

(* compute determinant and inverse *)
(* the number of dimensions can be 3 or 4 *) 
If [TensorEquationsDim === 4, imin = 0, imin = 1];

matrixarray[g_] := Array[g, {4-imin,4-imin}, {imin,imin}]
matrixdet[g_] := Det[matrixarray[g]]
matrixinv[g_, a_Integer, b_Integer] := 
  Inverse[matrixarray[g]] [[a-imin+1, b-imin+1]]
matrixinvdet[g_, a_Integer, b_Integer] := 
  (Det[matrixarray[g]] Inverse[matrixarray[g]])[[a-imin+1, b-imin+1]]


(* Kronecker delta *)
delta[a_,b_] := 0 /; a != b
delta[a_,b_] := 1 /; a == b


timer["after tensorrules.m"]
