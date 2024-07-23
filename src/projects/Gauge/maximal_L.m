(* maximal_L.m 
   Bernd Bruegmann 7/03 *)


(* variables *)
variables = {Lu, u, gi[a,b], G[a], KK}


(* compute in this order *)
tocompute = {

  du[a] == del[a,u],
  ddu[a,b] == deldel[a,b,u],

(* added 1 to u because someone set alpha -= 1. because of robin boundary condition *)
  Lu == gi[a,b] ddu[a,b] - G[a] du[a] - KK (1.+u)

}


(* symmetries *)
gi[a_,b_] := gi[b,a] /; !OrderedQ[{a,b}]
ddu[a_,b_] := ddu[b,a] /; !OrderedQ[{a,b}]
deldel[a_,b_,c_] := deldel[b,a,c] /; !OrderedQ[{a,b}]




(************************************************************************)
(* information for C output *)

(* information about the function 
   the C file will be functionname.c 
*)
includeANDdefine[] := Module[{},

  pr["#include \"bam.h\"\n"];
  pr["#include \"Gauge.h\"\n"];
  pr["\n"];
  pr["#define Power(x,y) (pow((double) (x), (double) (y)))\n"];
  pr["#define Log(x)     log((double) (x)))\n"];
  pr["#define pow2(x)    ((x)*(x))\n"];
  pr["#define pow2inv(x) (1.0/((x)*(x)))\n\n"];
  pr["#define Cal(x,y,z) ((x)?(y):(z))\n\n"];

  pr["\n\n\n"];
];

functionname = "maximal_L";

functionarguments = "tL *level, tVarList *vlLu, tVarList *vlu, tVarList *vlgi, tVarList *vlG, tVarList *vlKK";

functionloop    = "forinner19(level) {";
endfunctionloop = "} endfor;  /* loop i, j, k */";


(* custom variable declaration
   we have to translate between the tensor names and the C variables  
   we may need special variables, for example for parameters
*)
variabledeclarations[] := Module[{},

  prdecvl[{Lu}, "vlLu"];
  prdecvl[{u}, "vlu"];
  prdecvl[{gi[a,b]}, "vlgi"];
  prdecvl[{G[a]}, "vlG"];
  prdecvl[{KK}, "vlKK"];
  pr["\n"];
];    




(************************************************************************)
(* now we are ready to go *)

<< "../../math/MathToC/TensorEquationsToC.m"
