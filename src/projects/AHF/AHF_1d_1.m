(* AHF_1d_1.m 
   mth 06/09 *)

(* compute right hand side of BSSN equations *)


(* variables *)
variables = { g[a,b],A[a,b],K,phi, s[a],H,
             xp[a]
            }
              


(* compute in this order *)
tocompute = {
    
    S[a] == xp[a],
    
    norm == Exp[4*phi]*g[a,b]*S[a]*S[b],
    
    s[a] == xp[a]/sqrt[norm]
    
}



(* symmetries *)

g[a_,b_] :=  g[b,a] /; !OrderedQ[{a,b}]
A[a_,b_] :=  A[b,a] /; !OrderedQ[{a,b}]





(************************************************************************)
(* information for C output *)

(* information about the function 
   the C file will be functionname.c 
*)
includeANDdefine[] := Module[{},

  pr["#include \"bam.h\"\n"];
  pr["#include \"AHF.h\"\n"];
  pr["\n"];
  pr["#define Power(x,y) (pow((double) (x), (double) (y)))\n"];
  pr["#define Log(x)     log((double) (x)))\n"];
  pr["#define pow2(x)    ((x)*(x))\n"];
  pr["#define pow2inv(x) (1.0/((x)*(x)))\n\n"];
  pr["#define Cal(x,y,z) ((x)?(y):(z))\n\n"];

  pr["\n\n\n"];
];

functionname = "AHF_1d_1";


functionarguments = "tVarList *vars, tVarList *ahfvars";


functionloop = 
  "forinner25(level,oo2dx,oo2dy,oo2dz) {";
endfunctionloop = "} endforinner; /* loop i, j, k */";







(* custom variable declaration
   we have to translate between the tensor names and the C variables  
   we may need special variables, for example for parameters
*)
variabledeclarations[] := Module[{},

  pr["tL *level = ahfvars->level;\n"];
  pr["double *xp1  = Ptr(level, \"x\");\n"];
  pr["double *xp2  = Ptr(level, \"y\");\n"];
  pr["double *xp3  = Ptr(level, \"z\");\n"];
  pr["\n"];
  
 
  prdecvl[{ g[a,b],A[a,b],K,phi}, "vars"];
  prdecvl[{ s[a],H}, "ahfvars"];
 
  
  pr["\n"];
];    




(************************************************************************)
(* now we are ready to go *)

<< "../../math/MathToC/TensorEquationsToC.m"
