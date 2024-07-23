(* grhd_a2p.m
   mth 11/09 *)

   
(* compute right hand side of GRHD equations *)


(* variables *)
variables = { g[a,b],
              Madmrho,MadmS[a],MadmSS[a,b],MadmST,
              Mrho,Mepsl,Mp,Mv[a],Mvsqr
}


(* compute in this order *)
tocompute = {

  

  (* 3-velocity 
  vup[a]   == Mv[a],
  v[a]     == g[a,b] Mv[b],
  vsq      == vup[a] v[a],
  W        == 1/sqrt[1-vsq],
  Cinstruction == "if (vsq>GRHD.HRSC_VMAX) W = GRHD.HRSC_WlorMAX;",
  W2hrho   == W*W*(Mrho + Mrho*Mepsl + Mp),
  *)
  
  (* compute ADM vars 
  Madmrho          == (W2hrho - Mp),
  MadmS[a]         == (W2hrho*vup[a]),
  MadmSS[a,b]      == (W2hrho*v[a]*v[b]   + g[a,b]*Mp),
  MadmST           == (W2hrho*v[a]*vup[a] + 3     *Mp),
  *)

  Mv[a] == 0

}



(* symmetries *)

     g[a_,b_] :=  g[b,a] /; !OrderedQ[{a,b}]
     K[a_,b_] :=  K[b,a] /; !OrderedQ[{a,b}]
   Tij[a_,b_] :=  Tij[b,a] /; !OrderedQ[{a,b}]
   Kij[a_,b_] :=  Kij[b,a] /; !OrderedQ[{a,b}]
MadmSS[a_,b_] :=  MadmSS[b,a] /; !OrderedQ[{a,b}]

   delg[c_,a_,b_]    := delg[c,b,a] /; !OrderedQ[{a,b}]
   ginv[a_,b_]       := ginv[b,a] /; !OrderedQ[{a,b}]

MadmSSTF[a_,b_] :=  MadmSSTF[b,a] /; !OrderedQ[{a,b}]
      rA[a_,b_] := rA[b,a] /; !OrderedQ[{a,b}]
      tA[a_,b_] := tA[b,a] /; !OrderedQ[{a,b}]



(************************************************************************)
(* information for C output *)

(* information about the function 
   the C file will be functionname.c 
*)
includeANDdefine[] := Module[{},

  pr["#include \"bam.h\"\n"];
  pr["#include \"grhd.h\"\n"];
  pr["\n"];
  pr["#define Power(x,y) (pow((double) (x), (double) (y)))\n"];
  pr["#define Log(x)     log((double) (x)))\n"];
  pr["#define pow2(x)    ((x)*(x))\n"];
  pr["#define pow2inv(x) (1.0/((x)*(x)))\n\n"];
  pr["#define Cal(x,y,z) ((x)?(y):(z))\n\n"];

  pr["\n\n\n"];
];

functionname = "grhd_a2p";
functionarguments = "tL* level";

functionbegin     = "bampi_openmp_start\n";

functionloop      = "forinnerpoints_ijk_openmp(level) {";

endfunctionloop   = "} endfor_ijk_openmp; /* loop i, j, k */";

endfunction       = "bampi_openmp_stop\n";



(* custom variable declaration
   we have to translate between the tensor names and the C variables  
   we may need special variables, for example for parameters
*)
variabledeclarations[] := Module[{},
  
  prdecvarname[{g[a,b]}, "adm_gxx"];
  prdecvlf[{Madmrho,MadmS[a],MadmSS[a,b],MadmST}, "give_matterADMVars_vl"];
  prdecvlf[{Mrho,Mepsl,Mv[a],Mp,Mvsqr,detmetric}, "give_primVars_vl"];
  
  pr["\n"];
  
  
];




(************************************************************************)
(* now we are ready to go *)

(* assume that we are working on a box *)
BoxMode = True;


<< "../../math/MathToC/TensorEquationsToC.m"
