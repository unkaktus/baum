(* TensorEquationsToC.m *)
(* Bernd Bruegmann, 2/96, 10/02 Nina Jansen 01.03*)

(* Mathematica script:
   Given equations in tensor notation, write corresponding C code. 
   Assumes that tonsors are 4 D
*)  


(*************************)
(* apply tensor formulas *)

<< tensorrules.m


(************************************************)
(* replace abstract indices by explicit indices *)

<< expandindices4.m


(***********************************************************************)
(* put variables on grid and replace derivatives by finite differences *)

<< discretize.m


(************************************************)
(* optimize for best floating point performance *)

<< optimize.m


(****************)
(* write C code *)

<< writeC.m


