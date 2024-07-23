(* FiniteDifferences.m *)

(* Bernd Bruegmann 3/2006

   compute finite difference stencils using the unique interpolating polynomial
   script version of Mathematica notebook
*)

poly[x_, points_] := 
  Simplify[InterpolatingPolynomial[Map[{#, f[#]}&, points], x]]

dpoly[x_, points_] := 
  Simplify[D[poly[x, points /. x :> xdummy], x] /. xdummy :> x]

ddpoly[x_, points_] := 
  Simplify[D[poly[x, points /. x :> xdummy], x, x] /. xdummy :> x] 

pointlist[m_, n_] := Table[x + i*h, {i, m, n}]



(* first derivative *)

(* 3 point stencils, second order, O(h^2) *)
Print["\nfirst derivative, 3 point stencils\n"]
Print[MatrixForm[
  Table[dpoly[x, pointlist[center-2, center]], {center, 0, 2}]
]]


(* 4 point stencils, third order, O(h^3) *)
Print["\nfirst derivative, 4 point stencils\n"]
Print[MatrixForm[
  Table[dpoly[x, pointlist[center-3, center]], {center, 0, 3}]
]]


(* 5 point stencils, fourth order, O(h^4) *)
Print["\n\nfirst derivative, 5 point stencils\n"]
Print[MatrixForm[
  Table[dpoly[x, pointlist[center-4, center]], {center, 0, 4}]
]]



