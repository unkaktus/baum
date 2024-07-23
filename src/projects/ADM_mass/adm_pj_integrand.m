(* adm_pj_integrand.m 
   Mark Hannam 02/06 
 Based on adm_mass_integrand.m *)

(* compute integrand of the ADM mass *)


(* variables *)
  variables = {x[a], g[a,b], K[a,b], TrK, IP[a], IJ[a]}



(* compute in this order *)
tocompute = {

  (* Determinant of 3-metric *)
  det == matrixdet[g],

  (* inverse physical metric *)
  detginvf == 1/(matrixdet[g]),
  ginv[a,b] == detginvf matrixinvdet[g,a,b],

  (* R *)
  R == sqrt [delta[a,b] x[a] x[b]],

  (* Normal *)
  NN[a] == x[a]/R,

  (* Integrand of the ADM P^i *)
  IP[a] ==  1.0/8.0/PI sqrt[det] (ginv[b,c] K[d,b] - delta[d,c] TrK ) delta[a,d] NN[c] R*R,

  (* Integrand of the ADM J^i *)
  IJ[a] == 1.0/8.0/PI sqrt[det] (ginv[e,d] K[c,e] - delta[c,d] TrK ) epsilon[a,b,c] x[b] NN[d] R*R 
}


(* symmetries *)

epsilon[1,2,3] := 1
epsilon[a_,b_,c_] := -epsilon[a,c,b] /; !OrderedQ[{b,c}]
epsilon[a_,b_,b_] := 0
epsilon[a_,b_,c_] := -epsilon[c,b,a] /; !OrderedQ[{a,c}]
epsilon[a_,b_,a_] := 0
epsilon[a_,b_,c_] := -epsilon[b,a,c] /; !OrderedQ[{a,b}]
epsilon[a_,a_,c_] := 0

g[a_,b_]       := g[b,a] /; !OrderedQ[{a,b}]
K[a_,b_]       := K[b,a] /; !OrderedQ[{a,b}]
ginv[a_,b_]    := ginv[b,a] /; !OrderedQ[{a,b}]
delg[c_,a_,b_] := delg[c,b,a] /; !OrderedQ[{a,b}]

(************************************************************************)
(* information for C output *)

(* information about the function 
   the C file will be functionname.c 
*)
includeANDdefine[] := Module[{},

  pr["#include \"bam.h\"\n"];
  pr["#include \"ADM_mass.h\"\n"];
  pr["\n"];
  pr["#define Power(x,y) (pow((double) (x), (double) (y)))\n"];
  pr["#define pow2(x)    ((x)*(x))\n"];
  pr["#define Sqrt(x)    sqrt(x)\n"];
  pr["#define PI  3.1415926535897932\n\n"];
  pr["\n\n\n"];

];

functionname = "adm_pj_integrand";

functionarguments = "tL *level, int i_x, int i_g, int i_K, int i_TrK";

(*functionloop = "forinner7(level) {";
endfunctionloop = "} endfor; /* loop i, j, k */";
 *)

functionloop = "forinner19(level) {";
endfunctionloop = "} endforinner; /* loop i, j, k */";


(* custom variable declaration
   we have to translate between the tensor names and the C variables  
   we may need special variables, for example for parameters
*)
variabledeclarations[] := Module[{},

  prdecvar[{x[a]},      "i_x"];
  prdecvar[{g[a,b]},    "i_g"];
  prdecvar[{K[a,n]},    "i_K"];
  prdecvar[{TrK},       "i_TrK"];

  pr["double *IP1 = level->v[Ind(\"ADM_mass_Pxint\")];\n"];
  pr["double *IP2 = level->v[Ind(\"ADM_mass_Pyint\")];\n"];
  pr["double *IP3 = level->v[Ind(\"ADM_mass_Pzint\")];\n"];
  pr["double *IJ1 = level->v[Ind(\"ADM_mass_Jxint\")];\n"];
  pr["double *IJ2 = level->v[Ind(\"ADM_mass_Jyint\")];\n"];
  pr["double *IJ3 = level->v[Ind(\"ADM_mass_Jzint\")];\n"];
];    


(************************************************************************)
(* now we are ready to go *)

<< "../../math/MathToC/TensorEquationsToC.m"
