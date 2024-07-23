(* VectorLaplaceFlat.m 
   Bernd Bruegmann 11/01, 3/03 *)

(* compute flat vector Laplace operator *)


(* variables *)
variables = {v[a], Lv[a]}


(* compute in this order *)
tocompute = {
  oot == 1/3,
  ddv[a,b,c] == deldel[a,b,v[c]],
  Lv[a] == ddv[b,b,a] + oot ddv[a,b,b]
}


(* symmetries *)
ddv[a_,b_,c_] := ddv[b,a,c] /; !OrderedQ[{a,b}]
deldel[a_,b_,c_] := deldel[b,a,c] /; !OrderedQ[{a,b}]




(************************************************************************)
(* information for C output *)

(* information about the function 
   the C file will be functionname.c 
*)
includeANDdefine[] := Module[{},

  pr["#include \"bam.h\"\n"];
  pr["\n"];
  pr["#define Power(x,y) (pow((double) (x), (double) (y)))\n"];
  pr["#define Log(x)     log((double) (x)))\n"];
  pr["#define pow2(x)    ((x)*(x))\n"];
  pr["#define pow2inv(x) (1.0/((x)*(x)))\n\n"];
  pr["#define Cal(x,y,z) ((x)?(y):(z))\n\n"];

  pr["\n\n\n"];
];

functionname = "VectorLaplaceFlat";

functionarguments = "tL *level, tVarList *vllu, tVarList *vlu";

functionloop    = "forinner19(level) {";
endfunctionloop = "} endfor;  /* loop i, j, k */";



(* custom variable declaration
   we have to translate between the tensor names and the C variables  
   we may need special variables, for example for parameters
*)
variabledeclarations[] := Module[{},

  prdecvl[{Lv[a]}, "vllu"];
  prdecvl[{v[a]}, "vlu"];
  pr["\n"];

];    




(************************************************************************)
(* now we are ready to go *)

<< "../../math/MathToC/TensorEquationsToC.m"
