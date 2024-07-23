(* bssn_init.m 
   Bernd Bruegmann 10/98, 10/02
   Wolfgang Tichy  4/2004       *)

(* initialize BSSN variables from ADM variables *)


(* variables *)
variables = {gb[a,b], Kb[a,b], 
             gt[a,b], At[a,b], G[a], K, chi, gtinv[a,b],
             shiftDr, rr,rp,xp,yp,zp}


(* compute in this order *)
tocompute = {

  (* conformal factors *)
  detgb       == matrixdet[gb],
  psim4       == detgb^(-1/3),
  psi         == detgb^(1/12),

  (* basic rescaling *)
  gt[a,b]     == psim4 gb[a,b],
  Kt[a,b]     == psim4 Kb[a,b], 

  (* inverse bssn metric *)
  detgt       == matrixdet[gt],
  detgtinv    == 1/detgt,
  gtinv[a,b]  == detgtinv matrixinvdet[gt,a,b],

  (* additional variables *)
  K           == gtinv[a,b] Kt[a,b],
  A[a,b]      == Kb[a,b] - gb[a,b] K / 3,
  At[a,b]     == psim4 A[a,b],
  chi         == Power[psi,chipsipower],




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
    
    dgtinv[a,b,c]     == dgtinvSST[a,b,c],
    dchi[a]           == dchiSST[a],
    
  Cinstruction == " }",
  
  
  (* G *)
  G[a] == - dgtinv[b,b,a],


  (* location-dependent shiftdriver *)
  Cinstruction == "if (use_eta) {",
    
    psim2     == psi^(-2),
    dpsim2[a] == -2*chipsipower*psim2*dchi[a]/chi,
    absdpsim2 == Sqrt[gtinv[a,b]dpsim2[a]dpsim2[b]],
    
    Cinstruction == "shiftDr[ijk] = bssn_eta_set(xp[ijk],yp[ijk],zp[ijk], absdpsim2,psim2);",
  
  Cinstruction == "}"
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
  pr["#include \"bssn.h\"\n"];
  pr["\n"];
  pr["#define Power(x,y) (pow((double) (x), (double) (y)))\n"];
  pr["#define Sqrt(x)    sqrt((double) (x))\n"];
  pr["#define Log(x)     log((double) (x))\n"];
  pr["#define pow2(x)    ((x)*(x))\n"];
  pr["#define pow2inv(x) (1.0/((x)*(x)))\n\n"];
  pr["#define Tanh(x)    tanh(x)\n"];
  pr["#define Sech(x)    (1/cosh(x))\n"];
  pr["#define Cal(x,y,z) ((x)?(y):(z))\n\n"];

  pr["\n\n\n"];
];

functionname      = "bssn_init";

functionbegin     = "bampi_openmp_start\n";

functionarguments = "tVarList *ucur";

functionloop      = "forallpoints_ijk_openmp(level) {";

endfunctionloop   = "} endfor_ijk_openmp; /* loop i, j, k */";

endfunction       = "bampi_openmp_stop\n";


(* custom variable declaration
   we have to translate between the tensor names and the C variables  
   we may need special variables, for example for parameters
*)
variabledeclarations[] := Module[{},

  pr["tL *level = ucur->level;\n"];
  prdecvl[{gb[a,b],Kb[a,b], gt[a,b],chi,At[a,b],K,G[a],gtinv[a,b],shiftDr}, "ucur"];
  pr["\n"];

  pr["double chipsipower  = Getd(\"bssn_chi_psipower\");\n"];
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
  
  
  pr["const int use_eta   = Getv(\"bssn_use_eta\", \"yes\");\n"];
  pr["bssn_eta_init(level->grid);\n"];

];



(************************************************************************)
(* now we are ready to go *)

<< "../../math/MathToC/TensorEquationsToC.m"
