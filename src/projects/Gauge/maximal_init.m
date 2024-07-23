(* maximal_init.m 
   Bernd Bruegmann 7/03 *)


(* variables *)
variables = {g[a,b], K[a,b], rho, SS[a,b], ginv[a,b], G[a], KK, psi, dpop[a]}


(* compute in this order *)
tocompute = {

  (* partial derivatives *)
  delg[c,a,b] == del[c,g[a,b]],

  (* make transition to physical metric *)
  f == psi^4,
  delf[a] == 4 f dpop[a],
  delg[c,a,b] == f delg[c,a,b] + g[a,b] delf[c],

  (* inverse physical metric *)
  detginvf == 1/(f matrixdet[g]),
  ginv[a,b] == detginvf matrixinvdet[g,a,b],

  (* connection *)
  gammado[c,a,b] == 1/2 (delg[a,b,c] + delg[b,a,c] - delg[c,a,b]),
  gamma[c,a,b] == ginv[c,d] gammado[d,a,b], 
  G[a] == ginv[b,c] gamma[a,b,c],

  (* square of K_ab *)
  Cif == matter,
    KK == ginv[a,c] ginv[b,d] K[a,b] K[c,d] + 4*PI*(rho + ginv[b,a]*SS[a,b]),
  Cif == else,
    KK == ginv[a,c] ginv[b,d] K[a,b] K[c,d],
  Cif == end,

  (* normalize principal part *)
  ginv[a,b] == f ginv[a,b],
  G[a] == f G[a],
  KK == f KK
}




(* symmetries *)
g[a_,b_] := g[b,a] /; !OrderedQ[{a,b}]
K[a_,b_] := K[b,a] /; !OrderedQ[{a,b}]
SS[a_,b_] := SS[b,a] /; !OrderedQ[{a,b}]
ginv[a_,b_] := ginv[b,a] /; !OrderedQ[{a,b}]
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
  pr["#define Log(x)     log((double) (x)))\n"];
  pr["#define pow2(x)    ((x)*(x))\n"];
  pr["#define pow2inv(x) (1.0/((x)*(x)))\n\n"];
  pr["#define pow4(x)    ((x)*(x)*(x)*(x))\n"];
  pr["#define Cal(x,y,z) ((x)?(y):(z))\n\n"];

  pr["\n\n\n"];
];

functionname = "maximal_init";

functionarguments = "tL *level";

functionloop    = "forinner19(level) {";
endfunctionloop = "} endfor;  /* loop i, j, k */";


(* custom variable declaration
   we have to translate between the tensor names and the C variables  
   we may need special variables, for example for parameters
*)
variabledeclarations[] := Module[{},

  pr["int matter = Getv(\"physics\",\"matter\");\n\n"];

  prdecvarname[{g[a,b]},   "gxx"];
  prdecvarname[{K[a,b]},   "Kxx"];
  prdecvarname[{ginv[a,b]},"maximal_gixx"];
  prdecvarname[{G[a]},     "maximal_Gx"];
  prdecvarname[{KK},       "maximal_KK"];
  prdecvarname[{psi},      "psi"];
  prdecvarname[{dpop[a]},  "dpsiopsix"];

  pr["int index_rho,index_SSxx;\n"];
  pr["if (matter) {\n"];
    pr["index_rho = Ind(\"rho\");\n"];
    pr["index_SSxx = Ind(\"SSxx\");\n"];
  pr["} else {\n"];
    pr["printf(\" these are dummy variables, because c has some problems by initialicing variables inside an if statement\\n\");"];
    pr["index_rho = Ind(\"gxx\");\n"];
    pr["index_rho = Ind(\"gxx\");\n"];
  pr["}\n"];

  pr["double *rho = level->v[index_rho + 0];\n"];
  pr["double *SS11 = level->v[index_SSxx + 0];\n"];
  pr["double *SS12 = level->v[index_SSxx + 1];\n"];
  pr["double *SS13 = level->v[index_SSxx + 2];\n"];
  pr["double *SS22 = level->v[index_SSxx + 3];\n"];
  pr["double *SS23 = level->v[index_SSxx + 4];\n"];
  pr["double *SS33 = level->v[index_SSxx + 5];\n"];


  pr["\n"];
];    





(************************************************************************)
(* now we are ready to go *)

<< "../../math/MathToC/TensorEquationsToC.m"
