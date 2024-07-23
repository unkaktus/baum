(* Z4d_rhs_add_source.m 
   Wolfgang Tichy 08/2009 
   mth 11/2010*)

(* compute sources of Z4d equations *)


(* variables *)
variables = { g[a,b],  A[a,b],  G[a],  Khat,  chi,  Theta,  alpha,  beta[a], B[a],
             ng[a,b], nA[a,b], nG[a], nKhat, nchi, nTheta, nalpha, nbeta[a],nB[a],
             Madmrho, MadmS[a], MadmSS[a,b],MadmST,
             psi }

(* compute in this order *)
tocompute = {

  (* full metric and chi *)
  detg          == matrixdet[g],
  (* Cinstruction  == "if ((detg<0.) || (!finite(detg))) detg = 1.;", *)
  detginv       == 1/detg,
  ginv[a,b]     == detginv matrixinvdet[g,a,b],
  chiguard      == chiDivFloor,
  chiguarded    == chi,
  Cinstruction  == "if (chiguarded < chiguard) chiguarded = chiguard;",
  f             == oochipsipower log[chiguarded],
  psi4          == exp[4 f],
  psim4         == 1/psi4,

  metric[a,b]   == psi4  g[a,b],
  metricinv[a,b]== psim4 ginv[a,b],

  Hhatmatter    == - 16 PI Madmrho,


  (* compute matter RHS for A_ij, Khat and G^i *)
  MadmSSTF[a,b] == MadmSS[a,b] - 1/3*metric[a,b]*MadmST,
  rA[a,b]       == - psim4 8*PI*alpha MadmSSTF[a,b],
  rKhat         == 4*PI*alpha*(Madmrho+MadmST),
  rG[a]         == 2*alpha*(-8*PI*psi4 MadmS[a]),
  rTheta        == alpha Hhatmatter/2,

  (* bssn uses R = AA-2/3KK insted of R = AA-2/3KK+16*pi*rho *)
  r2A[a,b]      == - 1/3* g[a,b]*  16*PI*alpha*Madmrho,

  (* beta uses rG[a] -> need extra term *)
  betaF         == shiftgammacoeff alpha^shiftalphapower,
  rB[a]         == withB * (gamma0factor + gamma2factor betaF) * rG[a],


  (* hack src if something went wrong *)
  Cinstruction == "if (CheckForNANandINF(19, rA11,rA12,rA13,rA22,rA23,rA33, 
    r2A11,r2A12,r2A13,r2A22,r2A23,r2A33, rG1,rG2,rG3, rK,rKhat,rTheta, rB1,rB2,rB3)) {",
        rA[a,b] == 0,
        r2A[a,b]== 0,
        rK      == 0,
        rKhat   == 0,
        rG[a]   == 0,
        rTheta  == 0,
        rB[a]   == 0,
    Cinstruction == "}",


  (* add matter RHS *)
  Cif == addlinear, 
    nKhat   == nKhat   + c rKhat,
    nA[a,b] == nA[a,b] + c rA[a,b],
    nG[a]   == nG[a]   + c rG[a],
    nTheta  == nTheta  + c rTheta,
    nB[a]   == nB[a]   + c rB[a],
  Cif == else,
    nKhat   == nKhat   +   rKhat,
    nA[a,b] == nA[a,b] +   rA[a,b],
    nG[a]   == nG[a]   +   rG[a],
    nTheta  == nTheta  +   rTheta,
    nB[a]   == nB[a]   +   rB[a],
  Cif == end
}


(* symmetries *)
        g[a_,b_] :=         g[b,a] /; !OrderedQ[{a,b}]
        A[a_,b_] :=         A[b,a] /; !OrderedQ[{a,b}]
       ng[a_,b_] :=        ng[b,a] /; !OrderedQ[{a,b}]
       nA[a_,b_] :=        nA[b,a] /; !OrderedQ[{a,b}]
       rA[a_,b_] :=        rA[b,a] /; !OrderedQ[{a,b}]
      r2A[a_,b_] :=       r2A[b,a] /; !OrderedQ[{a,b}]

   MadmSS[a_,b_] :=    MadmSS[b,a] /; !OrderedQ[{a,b}]
 MadmSSTF[a_,b_] :=  MadmSSTF[b,a] /; !OrderedQ[{a,b}]

     ginv[a_,b_] :=      ginv[b,a] /; !OrderedQ[{a,b}]
   metric[a_,b_] :=    metric[b,a] /; !OrderedQ[{a,b}]
metricinv[a_,b_] := metricinv[b,a] /; !OrderedQ[{a,b}]



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
  pr["#define Log(x)     log((double) (x))\n"];
  pr["#define pow2(x)    ((x)*(x))\n"];
  pr["#define pow2inv(x) (1.0/((x)*(x)))\n\n"];
  pr["#define Cal(x,y,z) ((x)?(y):(z))\n\n"];

  pr["\n\n\n"];
];

functionname = "z4_rhs_add_source";

functionbegin     = "bampi_openmp_start\n";

functionarguments = "tVarList *unew, tVarList *upre, double c, tVarList *ucur";

functionloop = "forinnerpoints_ijk_openmp(level) {";

endfunctionloop = "} endfor_ijk_openmp; /* loop i, j, k */"

endfunction       = "bampi_openmp_stop\n";

(* custom variable declaration
   we have to translate between the tensor names and the C variables  
   we may need special variables, for example for parameters
*)
variabledeclarations[] := Module[{},

  pr["tL *level = ucur->level;\n"];
  pr["int addlinear = (c != 0.0l);\n"];
  pr["const double oochipsipower = 1.0/Getd(\"z4_chi_psipower\");\n"];
  pr["const double chiDivFloor = Getd(\"z4_chi_div_floor\");\n"];
  pr["\n"];
  
  prdecvarname[{Madmrho},     "adm_rho"];
  prdecvarname[{MadmS[a]},    "adm_Sx"];
  prdecvarname[{MadmSS[a,b]}, "adm_SSxx"];
  prdecvarname[{MadmST},      "adm_ST"];
  prdecvlNR[{ g[a,b],  chi,  A[a,b],  Khat,  G[a],  Theta,  alpha,  beta[a], B[a]}, "ucur", "METRIC_z4_INDX_VAR"];
  prdecvlNR[{ng[a,b], nchi, nA[a,b], nKhat, nG[a], nTheta, nalpha, nbeta[a],nB[a]}, "unew", "METRIC_z4_INDX_VAR"];
  
  
  pr["const double gamma0factor    = Getv(\"z4_shift\", \"gamma0\");\n"];
  pr["const double gamma2factor    = Getv(\"z4_shift\", \"gamma2\");\n"];
  pr["const double withGadv        = Getv(\"z4_shift\", \"withGadv\");\n"];
  pr["const double withShiftadv    = Getv(\"z4_shift\", \"withShiftadv\");\n"];
  pr["const double withBadv        = Getv(\"z4_shift\", \"withBadv\");\n"];
  pr["const double withB           = Getv(\"z4_shift\", \"withB\");\n"];

  pr["const double lapseharmonicf  = Getd(\"z4_lapseharmonicf\");\n"];
  pr["const double shiftalphapower = Getd(\"z4_shiftalphapower\");\n"];
  pr["const double shiftgammacoeff = Getd(\"z4_shiftgammacoeff\");\n"];
  pr["const double shiftdriver     = Getd(\"z4_shiftdriver\");\n"];
 
  pr["\n"];
];


(************************************************************************)
(* now we are ready to go *)
BoxMode = True;

<< "../../math/MathToC/TensorEquationsToC.m"

