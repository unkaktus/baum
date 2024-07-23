(* TensorEquationsToC.m *)
(* Bernd Bruegmann, 2/96, 10/02 *)

(* Mathematica script:
   Given equations in tensor notation, write corresponding C code.
*)  


(*************************)
(* apply tensor formulas *)

<< tensorrules.m



(************************************************)
(* replace abstract indices by explicit indices *)

<< expandindices.m

(*************************)
(* apply tensor formulas *)

<< jacobians.m

(***********************************************************************)
(* put variables on grid and replace derivatives by finite differences *)

<< discretize.m


(************************************************)
(* optimize for best floating point performance *)

<< optimize.m


(****************)
(* write C code *)

<< writeC.m


