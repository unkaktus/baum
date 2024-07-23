(* VectorLieFlat.m 
   Bernd Bruegmann 3/03 *)

(* compute flat Lie operator "l" *)


(* variables *)
variables = {v[a], lv[a,b]}


(* compute in this order *)
tocompute = {
  dv[a,b] == del[a,v[b]],
  trdv    == 2/3 dv[a,a],
  lv[a,b] == dv[a,b] + dv[b,a] - delta[a,b] trdv
}


(* symmetries *)
ddv[a_,b_,c_] := ddv[b,a,c] /; !OrderedQ[{a,b}]
deldel[a_,b_,c_] := deldel[b,a,c] /; !OrderedQ[{a,b}]
lv[a_,b_] := lv[b,a] /; !OrderedQ[{a,b}]



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

functionname = "VectorLieFlat";

functionarguments = "tL *level, tVarList *vllu, tVarList *vlu";

functionloop    = "forinner19(level) {";
endfunctionloop = "} endfor;  /* loop i, j, k */";



(* custom variable declaration
   we have to translate between the tensor names and the C variables  
   we may need special variables, for example for parameters
*)
variabledeclarations[] := Module[{},

  prdecvl[{lv[a,b]}, "vllu"];
  prdecvl[{v[a]}, "vlu"];
  pr["\n"];

];    




(************************************************************************)
(* now we are ready to go *)

<< "../../math/MathToC/TensorEquationsToC.m"
