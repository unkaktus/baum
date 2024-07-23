(* Z4d_init.m 
   Wolfgang Tichy  4/2004       
   mth 06/2011 *)

(* initialize Z4d variables from ADM variables *)


(* variables *)
variables = {gb[a,b], Kb[a,b], 
             gt[a,b], At[a,b], G[a], K, chi, gtinv[a,b], Theta,
             rr,rp,xp,yp,zp }


(* compute in this order *)
tocompute = {

  (* conformal factors *)
  detgb       == matrixdet[gb],
  psim4       == detgb^(-1/3),
  psi         == detgb^(1/12),

  (* basic rescaling *)
  gt[a,b]     == psim4 gb[a,b],
  Kt[a,b]     == psim4 Kb[a,b],

  (* inverse Z4d metric *)
  detgt       == matrixdet[gt],
  detgtinv    == 1/detgt,
  gtinv[a,b]  == detgtinv matrixinvdet[gt,a,b],
  
  Cinstruction == "if (detgt<=0.00001) printf(\"detg<0 :  %e %e %e\\n\",detgb,detgt,chi[ijk]);",

  (* additional variables *)
  K           == gtinv[a,b] Kt[a,b],
  A[a,b]      == Kb[a,b] - gb[a,b] K / 3,
  At[a,b]     == psim4 A[a,b],
  chi         == Power[psi,chipsipower],
  Theta       == 0,




  (* now without boundary *)
  Cinstruction == "} endfor_ijk_openmp;",
  Cinstruction == "forinnerpoints_ijk_openmp(level) {",
  
  dgtinv[c,a,b] == del[c,gtinv[a,b]],
  dchi[a]       == del[a,chi],
  
  Cinstruction == "if (useShellsTransfo) {",
    
    (* derivatives for fisheye coords ... this is a speedup *)
    dRdr         == DRDr,
    ddRdr        == DDRDr,
    
    (* define jacobians for each box *)
    Cinstruction == "if (nbox<=1 ) {",
      Jac[a,b]    == JacobB01[a,b],
    Cinstruction == "} else if (nbox<=3) {",
      Jac[a,b]    == JacobB23[a,b],
    Cinstruction == "} else if (nbox<=5) {",
      Jac[a,b]    == JacobB45[a,b],
    Cinstruction == "}",

    dgtinvSST[a,b,c]  == Jac[d,a] dgtinv[d,b,c],
    dchiSST[a]        == Jac[b,a] dchi[b],
    
    dgtinv[a,b,c]    == dgtinvSST[a,b,c],
    dchi[a]           == dchiSST[a],
    
  Cinstruction == "}",
  
  
  
  (* G *)
  G[a] == - dgtinv[b,b,a]
  
}


(* symmetries *)
gb[a_,b_] := gb[b,a] /; !OrderedQ[{a,b}]
Kb[a_,b_] := Kb[b,a] /; !OrderedQ[{a,b}]

gt[a_,b_] := gt[b,a] /; !OrderedQ[{a,b}]
gtinv[a_,b_] := gtinv[b,a] /; !OrderedQ[{a,b}]
gtinvSST[a_,b_] := gtinvSST[b,a] /; !OrderedQ[{a,b}]
Kt[a_,b_] := Kt[b,a] /; !OrderedQ[{a,b}]
At[a_,b_] := At[b,a] /; !OrderedQ[{a,b}]

dgtinv[c_,a_,b_] := dgtinv[c,b,a] /; !OrderedQ[{a,b}]




(************************************************************************)
(* information for C output *)

(* information about the function 
   the C file will be functionname.c 
*)
includeANDdefine[] := Module[{},

  pr["#include \"bam.h\"\n"];
  pr["#include \"z4.h\"\n"];
  pr["\n"];
  pr["#define Power(x,y) (pow((double) (x), (double) (y)))\n"];
  pr["#define Log(x)     log((double) (x)))\n"];
  pr["#define pow2(x)    ((x)*(x))\n"];
  pr["#define pow2inv(x) (1.0/((x)*(x)))\n\n"];
  pr["#define Tanh(x)    tanh(x)\n"];
  pr["#define Sech(x)    (1/cosh(x))\n"];
  pr["#define Cal(x,y,z) ((x)?(y):(z))\n\n"];

  pr["\n\n\n"];
];

functionname = "z4_init";

functionarguments = "tVarList *ucur";

functionbegin     = "bampi_openmp_start\n";

functionloop      = "forallpoints_ijk_openmp(level) {";

endfunctionloop   = "} endfor_ijk_openmp; /* loop i, j, k */";

endfunction       = "bampi_openmp_stop\n";


(* custom variable declaration
   we have to translate between the tensor names and the C variables  
   we may need special variables, for example for parameters
*)
variabledeclarations[] := Module[{},

  pr["tL *level = ucur->level;\n"];
  prdecvlNR[{gb[a,b],Kb[a,b], gt[a,b],chi,At[a,b],K,G[a],gtinv[a,b],Theta}, "ucur", "METRIC_z4_INDX_VAR"];
  pr["\n"];
 
  pr["double chipsipower = Getd(\"z4_chi_psipower\");\n"];
  pr["\n"];
  
  (* where are we? *)
  pr["const double *xp = level->v[Ind(\"x\")];\n"];
  pr["const double *yp = level->v[Ind(\"y\")];\n"];
  pr["const double *zp = level->v[Ind(\"z\")];\n"];
  
  (* shells stuff *)
  pr["const int useShellsTransfo = level->shells;\n"];
  pr["const double *rp = level->v[IndLax(\"shells_R\")];\n"];
  pr["const double *rr = level->v[IndLax(\"shells_r\")];\n"];
  pr["const double shellsS = GetdLax(\"amr_shells_stretch\");\n"];
  pr["const double shellsR = GetdLax(\"amr_shells_r0\");\n"];
  pr["const double shellsE = GetdLax(\"amr_shells_eps\");\n"];

];    




(************************************************************************)
(* now we are ready to go *)

(* assume that we are working on a box *)
BoxMode = True;

<< "../../math/MathToC/TensorEquationsToC.m"



