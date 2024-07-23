(* AHF_1d_2.m 
   mth 06/09 *)

(* compute right hand side of BSSN equations *)


(* variables *)
variables = { g[a,b],A[a,b],K,phi, s[a],H
            }
              


(* compute in this order *)
tocompute = {
    
    delg[c,a,b]         == del[c,g[a,b]],
    deldelg[a,b,c,d]    == deldel[a,b,g[c,d]],
    detginv             == 1/matrixdet[g],
    ginv[a,b]           == detginv matrixinvdet[g,a,b],
    gammado[c,a,b]      == 1/2 (delg[a,b,c] + delg[b,a,c] - delg[c,a,b]),
    gamma[c,a,b]        == ginv[c,d] gammado[d,a,b], 
    
    ds[a,b]        == del[a,s[b]],
    dphi[a]        == del[a,phi],
    ggamma[c,a,b]  == gamma[c,a,b] + 2*(delta[c,a]*dphi[b] + delta[c,b]*dphi[a] - g[a,b]*ginv[c,d]*dphi[d]),
   
    Ds[a,b]        == ds[a,b] + s[c]*ggamma[b,a,c],
    KK[a,b]        == exp[4*phi]*(A[a,b] + 1/3*g[a,b]*K),
    H              == Ds[a,a] - K + KK[a,b]*s[a]*s[b]
   
   
}



(* symmetries *)

g[a_,b_] :=  g[b,a] /; !OrderedQ[{a,b}]
A[a_,b_] :=  A[b,a] /; !OrderedQ[{a,b}]
 
delg[c_,a_,b_] := delg[c,b,a] /; !OrderedQ[{a,b}]
deldelg[a_,b_,c_,d_] := deldelg[b,a,c,d] /; !OrderedQ[{a,b}]
deldelg[a_,b_,c_,d_] := deldelg[a,b,d,c] /; !OrderedQ[{c,d}]
ginv[a_,b_] := ginv[b,a] /; !OrderedQ[{a,b}]

ds[a_,b_] :=  ds[b,a]/; !OrderedQ[{a,b}]
Ds[a_,b_] :=  Ds[b,a]/; !OrderedQ[{a,b}]
KK[a_,b_] :=  KK[b,a]/; !OrderedQ[{a,b}]




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

functionname = "AHF_1d_2";


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
  
  prdecvl[{ g[a,b],A[a,b],K,phi}, "vars"];
  prdecvl[{ s[a],H}, "ahfvars"];
  
  pr["\n"];
];    




(************************************************************************)
(* now we are ready to go *)

<< "../../math/MathToC/TensorEquationsToC.m"
