(* VectorLaplaceFlatGS.m 
   Bernd Bruegmann 11/01, 3/03 *)

(* compute Gauss-Seidel step with flat vector Laplace operator *)


(* variables *)
variables = {v[a], u[a], f[a]}


(* compute in this order *)
tocompute = {
  oot == 1/3,
  ddu[a,b,c] == deldel[a,b,u[c]],
  Lu[a] == ddu[b,b,a] + oot ddu[a,b,b],
  Lii == - 20 oot oo2dx,
  v[a] == u[a] + (f[a] - Lu[a])/Lii
}


(* symmetries *)
ddu[a_,b_,c_] := ddu[b,a,c] /; !OrderedQ[{a,b}]
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

functionname = "VectorLaplaceFlatGS";

functionarguments = "tL *level, tVarList *vlv, tVarList *vlu, tVarList *vlf";

functionloop    = "forinner19(level) {";
endfunctionloop = "} endfor;  /* loop i, j, k */";



(* custom variable declaration
   we have to translate between the tensor names and the C variables  
   we may need special variables, for example for parameters
*)
variabledeclarations[] := Module[{},

  prdecvl[{v[a]}, "vlv"];
  prdecvl[{u[a]}, "vlu"];
  prdecvl[{f[a]}, "vlf"];
  pr["\n"];

];    




(************************************************************************)
(* now we are ready to go *)

<< "../../math/MathToC/TensorEquationsToC.m"
