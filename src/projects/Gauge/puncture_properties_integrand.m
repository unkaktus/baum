(* puncture_properties_Integrand.m 
   dtim 11/15 *)

(* variables *)
  variables = {x[a], g[a,b], K[a,b], TrK, IM, IP[a],IS[a]}

(* compute in this order *)
tocompute = {

  (* partial derivatives *)
  delg[c,a,b] == del[c,g[a,b]],

  (* Determinant of 3-metric *)
  det == matrixdet[g],

  (* inverse physical metric *)
  detginvf == 1/(matrixdet[g]),
  ginv[a,b] == detginvf matrixinvdet[g,a,b],

  R == 0,
  rho ==0, 
  Cinstruction ==  "R =sqrt((x1[ijk]-x0)*(x1[ijk]-x0)+(x2[ijk]-y0)*(x2[ijk]-y0)+(x3[ijk]-z0)*(x3[ijk]-z0));",
  Cinstruction ==  "rho =sqrt((x1[ijk]-x0)*(x1[ijk]-x0)+(x2[ijk]-y0)*(x2[ijk]-y0));",
  xp[a]    == 0,
  ones[a]  == 1,
  Cinstruction ==  "xp1 = (x1[ijk]-x0);",
  Cinstruction ==  "xp2 = (x2[ijk]-y0);",
  Cinstruction ==  "xp3 = (x3[ijk]-z0);",
  
  NN[a]       == xp[a]/R,
  phi[a,b]   == epsilon[a,c,e] delta[b,c] delta[d,e] xp[d],

(*
  dxdr[a,b] == 0,
  dxdr11 == xp1/R,
  dxdr12 == (xp1 xp3)/rho,
  dxdr13 == -xp2,
  dxdr21 == xp2/R,
  dxdr22 == (xp2 xp3)/rho,
  dxdr23 == xp1,
  dxdr31 == xp3/R,
  dxdr32 == -rho,
  dxdr33 == 0,

  gtt == dxdr[a,2] dxdr[b,2] g[a,b],
  gtp == dxdr[a,2] dxdr[b,3] g[a,b],
  gpp == dxdr[a,3] dxdr[b,3] g[a,b],
*)

  IM    ==  1.0/16.0/PI sqrt[det] ginv[a,b] (delg[b,a,c] - delg[c,a,b]) xp[c] R, 
  IP[a] ==  1.0/8.0/PI sqrt[det] (ginv[b,c] K[d,b] - delta[d,c] TrK ) delta[a,d] xp[c] R,
  IS[a] ==  1.0/8.0/PI sqrt[det] (ginv[e,d] K[c,e] - delta[c,d] TrK ) phi[c,a] xp[d] R
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
  pr["#include \"Gauge.h\"\n"];
  pr["\n"];
  pr["#define Power(x,y) (pow((double) (x), (double) (y)))\n"];
  pr["#define pow2(x)    ((x)*(x))\n"];
  pr["#define pow3(x)    ((x)*(x)*(x))\n"];
  pr["#define pow2inv(x) (1.0/((x)*(x)))\n\n"];
  pr["#define Sqrt(x)    sqrt(x)\n"];
  pr["#define PI  3.1415926535897932\n\n"];
  pr["\n\n\n"];

];

functionname = "puncture_properties_integrand";

functionarguments = "tL *level, int i_x, int i_g, int i_K, int i_TrK, double x0,double y0, double z0"

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

  pr["double *IM  = level->v[Ind(\"puncture_properties_Mint\")];\n"];

  pr["double *IP1 = level->v[Ind(\"puncture_properties_Pxint\")];\n"];
  pr["double *IP2 = level->v[Ind(\"puncture_properties_Pyint\")];\n"];
  pr["double *IP3 = level->v[Ind(\"puncture_properties_Pzint\")];\n"];

  pr["double *IS1 = level->v[Ind(\"puncture_properties_Sxint\")];\n"];
  pr["double *IS2 = level->v[Ind(\"puncture_properties_Syint\")];\n"];
  pr["double *IS3 = level->v[Ind(\"puncture_properties_Szint\")];\n"];

];    


(************************************************************************)
(* now we are ready to go *)

<< "../../math/MathToC/TensorEquationsToC.m"
