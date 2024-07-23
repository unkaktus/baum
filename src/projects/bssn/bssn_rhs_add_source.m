(* bssn_rhs_add_source.m
   mth 08/09 *)

(* compute sources of BSSN equations *)


(* variables *)
variables = { g[a,b],  A[a,b],  G[a],  K,  chi,  alpha,  beta[a],  B[a],
             ng[a,b], nA[a,b], nG[a], nK, nchi, nalpha, nbeta[a], nB[a],
             Madmrho,MadmS[a],MadmSS[a,b],MadmST, xp,yp,zp
            }


(* compute in this order *)
tocompute = {
    
    (* full metric *)
    detg            == matrixdet[g],
    detginv         == 1/detg,
    ginv[a,b]       == detginv matrixinvdet[g,a,b],
  
    chiguard        == chiDivFloor,
    chiguarded      == chi,
    Cinstruction    == "if (chiguarded < chiguard) chiguarded = chiguard;",
    ff              == chiguarded,
    oochipsipower   == 1/chipsipower,
    f               == oochipsipower log[ff],
    psim4           == exp[-4 f],
    psi4            == exp[ 4 f],
      
    metric[a,b]     == g[a,b]*psi4,
    metricinv[a,b]  == ginv[a,b]*psim4,

    
    Cinstruction    == "//if (detg<=0.) errorexit(\" detg <= 0.\"); ",
    Cinstruction    == "//if (psi4<=0.) errorexit(\" psi^4 <= 0.\"); ",
    
    
    (* compute RHSs for A K and G *)
    MadmSSTF[a,b]   == MadmSS[a,b] - 1/3*metric[a,b]*MadmST,
    rA[a,b]         == - psim4 8*PI*alpha MadmSSTF[a,b],
    rK              == 4*PI*alpha*(Madmrho+MadmST),
    rG[a]           == 2*alpha*(-8*PI*psi4 MadmS[a]),
    
    (* bssn uses R = AA-2/3KK insted of R = AA-2/3KK+16*pi*rho *)
    r2A[a,b]        == - 1/3* g[a,b]*  16*PI*alpha*Madmrho,
    
    (* beta uses rG[a] -> need extra term *)
    betaF           == shiftgammacoeff alpha^shiftalphapower,
    rB[a]           == withB * (gamma0factor + gamma2factor betaF) * rG[a],
    
    
    
    Cinstruction == "if (setRHSto0 || CheckForNANandINF(19, rA11,rA12,rA13,rA22,rA23,rA33, 
    r2A11,r2A12,r2A13,r2A22,r2A23,r2A33, rG1,rG2,rG3,rK, rB1,rB2,rB3)) {",
        rA[a,b] == 0,
        rK      == 0,
        rG[a]   == 0,
        r2A[a,b]== 0,
        rB[a]   == 0,
    Cinstruction == "}",
    
    
    (*
    Cinstruction == "if (CheckForNANandINF(19, rA11,rA12,rA13,rA22,rA23,rA33, 
    r2A11,r2A12,r2A13,r2A22,r2A23,r2A33, rG1,rG2,rG3,rK, rB1,rB2,rB3))
    printf(\"problem with nans inside bssn_rhs_add_source\\n\");",
    *)
    
    
    (* result *)
    Cif == addlinear, 
        nK      == nK      + c rK,
        nA[a,b] == nA[a,b] + c(rA[a,b]+r2A[a,b]),
        nG[a]   == nG[a]   + c rG[a],
        nB[a]   == nB[a]   + c rB[a],
    Cif == else,
        nK      == nK      +   rK,
        nA[a,b] == nA[a,b] +  (rA[a,b]+r2A[a,b]),
        nG[a]   == nG[a]   +   rG[a],
        nB[a]   == nB[a]   +   rB[a],
    Cif == end,
    
    Cinstruction == "if (CheckForNANandINF(13, nA11[ijk],nA12[ijk],nA13[ijk],nA22[ijk],nA23[ijk],nA33[ijk], nG1[ijk],nG2[ijk],nG3[ijk], nK[ijk], nB1[ijk],nB2[ijk],nB3[ijk]))
    printf(\"the ultimative test if there are no nans coming from the matter part ... failed (nans inside bssn_rhs_add_source)\\n\");"
    
}



(* symmetries *)
        g[a_,b_] :=         g[b,a] /; !OrderedQ[{a,b}]
        A[a_,b_] :=         A[b,a] /; !OrderedQ[{a,b}]
       ng[a_,b_] :=        ng[b,a] /; !OrderedQ[{a,b}]
       nA[a_,b_] :=        nA[b,a] /; !OrderedQ[{a,b}]
   MadmSS[a_,b_] :=    MadmSS[b,a] /; !OrderedQ[{a,b}]

     ginv[a_,b_] :=      ginv[b,a] /; !OrderedQ[{a,b}]
   metric[a_,b_] :=    metric[b,a] /; !OrderedQ[{a,b}]
metricinv[a_,b_] := metricinv[b,a] /; !OrderedQ[{a,b}]

 MadmSSTF[a_,b_] :=  MadmSSTF[b,a] /; !OrderedQ[{a,b}]
       rA[a_,b_] :=        rA[b,a] /; !OrderedQ[{a,b}]
      r2A[a_,b_] :=       r2A[b,a] /; !OrderedQ[{a,b}]




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
  pr["#define Log(x)     log((double) (x)))\n"];
  pr["#define pow2(x)    ((x)*(x))\n"];
  pr["#define pow2inv(x) (1.0/((x)*(x)))\n\n"];
  pr["#define Cal(x,y,z) ((x)?(y):(z))\n\n"];
  

  pr["\n\n\n"];
];

functionname      = "bssn_rhs_add_source";

functionarguments = "tVarList *unew, tVarList *upre, double c, tVarList *ucur";

functionbegin     = "bampi_openmp_start\n";

functionloop      = "forinnerpoints_ijk_openmp(level) {";

endfunctionloop   = "} endfor_ijk_openmp; /* loop i, j, k */";

endfunction       = "bampi_openmp_stop\n";










(* custom variable declaration
   we have to translate between the tensor names and the C variables  
   we may need special variables, for example for parameters
*)
variabledeclarations[] := Module[{},

  pr["tL *level = ucur->level;\n"];
  pr["int addlinear = (c != 0.0l);\n"];
  pr["\n"];
  
  pr["const double chiDivFloor = Getd(\"bssn_chi_div_floor\");\n"];
  pr["const double chipsipower = Getd(\"bssn_chi_psipower\");\n"];
  pr["\n"];
  
  prdecvarname[{Madmrho},     "adm_rho"];
  prdecvarname[{MadmS[a]},    "adm_Sx"];
  prdecvarname[{MadmSS[a,b]}, "adm_SSxx"];
  prdecvarname[{MadmST},      "adm_ST"];
  prdecvlNR[{ g[a,b],  chi,  A[a,b],  K,  G[a],  alpha,  beta[a],  B[a]}, "ucur", "METRIC_bssn_INDX_VAR"];
  prdecvlNR[{ng[a,b], nchi, nA[a,b], nK, nG[a], nalpha, nbeta[a], nB[a]}, "unew", "METRIC_bssn_INDX_VAR"];
  
  pr["const double *xp = level->v[Ind(\"x\")];\n"];
  pr["const double *yp = level->v[Ind(\"y\")];\n"];
  pr["const double *zp = level->v[Ind(\"z\")];\n"];
  
  
  pr["const double gamma0factor    = Getv(\"bssn_shift\", \"gamma0\");\n"];
  pr["const double gamma2factor    = Getv(\"bssn_shift\", \"gamma2\");\n"];
  pr["const double withGadv        = Getv(\"bssn_shift\", \"withGadv\");\n"];
  pr["const double withShiftadv    = Getv(\"bssn_shift\", \"withShiftadv\");\n"];
  pr["const double withBadv        = Getv(\"bssn_shift\", \"withBadv\");\n"];
  pr["const double withB           = !Getv(\"bssn_shift\", \"withoutB\");\n"];

  pr["const double lapseharmonicf  = Getd(\"bssn_lapseharmonicf\");\n"];
  pr["const double shiftalphapower = Getd(\"bssn_shiftalphapower\");\n"];
  pr["const double shiftgammacoeff = Getd(\"bssn_shiftgammacoeff\");\n"];
  pr["const double shiftdriver     = Getd(\"bssn_shiftdriver\");\n"];

  
  
  
  
  
  
  pr["const int setRHSto0 = Getv(\"bssnRHSto0\",\"yes\")+ Getv(\"bssnsourceRHSto0\",\"yes\");\n"];
  pr["\n"];
];




(************************************************************************)
(* now we are ready to go *)
BoxMode = True;


<< "../../math/MathToC/TensorEquationsToC.m"
