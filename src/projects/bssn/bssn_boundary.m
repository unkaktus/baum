
(************************************************************************)
(* Sommerfeld boundary conditions *)
(************************************************************************)

(* bssn_boundary.m mth 02/12 *)

(* variables *)
variables = { g[a,b],  chi,  A[a,b],  K,  G[a],  alpha,  beta[a],  B[a],  MD,  MS[a],  Mtau,
             ng[a,b], nchi, nA[a,b], nK, nG[a], nalpha, nbeta[a], nB[a], nMD, nMS[a], nMtau,
             pg[a,b], pchi, pA[a,b], pK, pG[a], palpha, pbeta[a], pB[a], pMD, pMS[a], pMtau,
             rp,rr,xp,yp,zp
            }

(* compute in this order *)
tocompute = {

  Cinstruction == "if (CheckForNANandINF(24, 
    g11[ijk],g12[ijk],g13[ijk],g22[ijk],g23[ijk],g33[ijk],
    A11[ijk],A12[ijk],A13[ijk],A22[ijk],A23[ijk],A33[ijk], 
    G1[ijk],G2[ijk],G3[ijk], K[ijk],chi[ijk],
    alpha[ijk],beta1[ijk],beta2[ijk],beta3[ijk],B1[ijk],B2[ijk],B3[ijk])) {
    printf(\"problem with vars in bssn_boundary.m\\n\");
    printf(\"x=%2.5e, y=%2.5e, z=%2.5e\\n\",xp[ijk],yp[ijk],zp[ijk]);}", 

(************************************************************************)
(* ---- Compute derivatives on the spherical grid ---- *)
(************************************************************************)

(* ---- derivatives of cartesian components in spherical directions (del_rpt V^xyz) ---- *)
    Cinstruction     == "if (order == 2 || boundaryNaway(1)) {",
        db[a,b]     == del[a,beta[b]],
        dB[a,b]     == del[a,B[b]],
        dg[c,a,b]   == del[c,g[a,b]],
        dK[a]       == del[a,K],
        dA[a,b,c]   == del[a,A[b,c]],
        dG[a,b]     == del[a,G[b]],
    Cinstruction == "} else if (order == 4 || boundaryNaway(2)) {",
        db[a,b]     == del4[a,beta[b]],
        dB[a,b]     == del4[a,B[b]],
        dg[c,a,b]   == del4[c,g[a,b]],
        dK[a]       == del4[a,K],
        dA[a,b,c]   == del4[a,A[b,c]],
        dG[a,b]     == del4[a,G[b]],
    Cinstruction == "} else if (order == 6 || boundaryNaway(3)) {",
        db[a,b]     == del6[a,beta[b]],
        dB[a,b]     == del6[a,B[b]],
        dg[c,a,b]   == del6[c,g[a,b]],
        dK[a]       == del6[a,K],
        dA[a,b,c]   == del6[a,A[b,c]],
        dG[a,b]     == del6[a,G[b]],
    Cinstruction == "} else if (order == 8 || boundaryNaway(4)) {",
        db[a,b]     == del8[a,beta[b]],
        dB[a,b]     == del8[a,B[b]],
        dg[c,a,b]   == del8[c,g[a,b]],
        dK[a]       == del8[a,K],
        dA[a,b,c]   == del8[a,A[b,c]],
        dG[a,b]     == del8[a,G[b]],
    Cinstruction == "} else if (order == 10 || boundaryNaway(5)) {",
        db[a,b]     == del10[a,beta[b]],
        dB[a,b]     == del10[a,B[b]],
        dg[c,a,b]   == del10[c,g[a,b]],
        dK[a]       == del10[a,K],
        dA[a,b,c]   == del10[a,A[b,c]],
        dG[a,b]     == del10[a,G[b]],
    Cinstruction == "} else errorexit(\"order not implemented\"); ",

(************************************************************************)
(* ---- Transform all derivatives to Cartesian ---- *)
(************************************************************************)

  (* ---- Define Jacobians for each box ---- *)
    
    (* derivatives for fisheye coords ... this is a speedup *)
    dRdr         == DRDr,
    ddRdr        == DDRDr,
    
    (* define jacobians for each box *)
    Cinstruction == "if (nbox<=1 ) {",
      Jac[a,b]    == JacobB01[a,b],
      DJac[a,b,c] == DJacobB01[a,b,c],
    Cinstruction == "} else if (nbox<=3) {",
      Jac[a,b]    == JacobB23[a,b],
      DJac[a,b,c] == DJacobB23[a,b,c],
    Cinstruction == "} else if (nbox<=5) {",
      Jac[a,b]    == JacobB45[a,b],
      DJac[a,b,c] == DJacobB45[a,b,c],
    Cinstruction == "}",

    (* ---- do the transformation ---- *)
    dbSST[a,b]      == Jac[c,a] db[c,b],
    dBSST[a,b]      == Jac[c,a] dB[c,b],
    dgSST[a,b,c]    == Jac[d,a] dg[d,b,c],
    dGSST[a,b]      == Jac[c,a] dG[c,b],
    dKSST[a]        == Jac[b,a] dK[b],
    dASST[a,b,c]    == Jac[d,a] dA[d,b,c],

    (* now give the tmp derivative the same name again *)
    db[a,b]       == dbSST[a,b],
    dB[a,b]       == dBSST[a,b],
    dg[a,b,c]     == dgSST[a,b,c],
    dG[a,b]       == dGSST[a,b],
    dK[a]         == dKSST[a],
    dA[a,b,c]     == dASST[a,b,c],

(************************************************************************)
(* ---- Define projection stuff in xyz and rpt; up and down differ! ---- *)
(************************************************************************)

  (* ---- Define normal vector ---- *)
    r    == 0,
    shat[a] == 0,
    sup[a] == 0,
    Cinstruction == "r=rp[ijk];",
    Cinstruction == "shat1=xp[ijk]/r;shat2=yp[ijk]/r;shat3=zp[ijk]/r;",
  
    detginv == 1/matrixdet[g],
  
    ginv[a,b]    == detginv matrixinvdet[g,a,b],
    ADMginv[a,b] == chi ginv[a,b],
    modshatARG   == ADMginv[c,d]shat[c]shat[d],
    Cinstruction == "if (modshatARG<0.00001) {
      printf(\"modshat is wrong (%e)\\n\",modshatARG);
      modshatARG = 0.00001;
    }",
    oomodshat    == 1/sqrt[modshatARG],
    sdown[a]     == oomodshat shat[a],
    sup[a]       == ADMginv[a,b] sdown[b],
    
  (* ---- Define 2+1 projection operator ---- *)
    qud[a,b] == delta[a,b] - sup[a]sdown[b],
    qdd[a,b] == g[a,b] / chi - sdown[a]sdown[b],
    quu[a,b] == ADMginv[a,b] - sup[a]sup[b],

  (* ---- Physical projection ---- *)
    qPhysuudd[c,d,a,b] == (qud[c,a]qud[d,b]+qud[c,b]qud[d,a]-quu[c,d]qdd[a,b])/2, 

  (* ---- Speeds for the boundary conditions ---- *)
    (* ---- lapse: muL ---- *)
    muL    ==  2 / alpha,
    (* ---- Shift: muS ---- *)
    (* ---- For compatibility with the standard BAM pars this should be changed! DH.  ---- *)
    muS    ==  shiftgammacoeff,
    (* ---- Longtitudinal shift speed ---- *)
    vbetas == 2 sqrt[ muS / (3 chi) ],
    (* ---- Transverse shift speed ---- *)
    vbetaA == sqrt[muS / chi],

(************************************************************************)
(* ---- Connections, curvatures, divs, aux vars ---- *)
(************************************************************************)

    (* tf part of conformal extrinsic curvature *)
    Ainv[a,b] == ginv[a,c] ginv[b,d] A[c,d],
    
    (* the shift terms *)
    divbeta == delta[a,b]db[a,b],
    totdivbeta == 2/3 divbeta,
    
    (* Lie derivative *)
    lieg[a,b] == beta[c]dg[c,a,b] - g[a,b] totdivbeta + db[a,c] g[b,c] + db[b,c] g[a,c],

(************************************************************************)
(* ---- Define some 2+1 split variables ---- *)
(************************************************************************)

    (* ---- Scalar sector  ---- *)
    ACss  == sup[a]sup[b]A[a,b],
    DACss == sup[c]sup[a]sup[b]dA[c,a,b],
    DK    == sup[a]dK[a],
    Bs    == sdown[a]B[a],
    DBs   == sup[a]sdown[b]dB[a,b],
    Gams  == sdown[a]G[a],
    DGams == sup[a]sdown[b]dG[a,b],
    
    (* ---- Vector sector  ---- *)
    BA[k] == qud[k,a]B[a],
    DBA[k] == sup[a]qud[k,b]dB[a,b],
    GamA[k] == qud[k,a]G[a],
    DGamA[k] == sup[a]qud[k,b]dG[a,b],
    ACsA[k] == sup[a]qud[b,k]A[a,b],
    DACsA[k] == sup[c]sup[a]qud[b,k]dA[c,a,b],
    
    (* ---- Tensor sector  ---- *)
    ACABTF[k,l]  == qPhysuudd[c,d,k,l]A[c,d],
    DACABTF[k,l] == sup[a]qPhysuudd[b,c,k,l]dA[a,b,c],

(************************************************************************)
(* ---- Construct RHSs of the 2+1 decomposed quantities ---- *)
(************************************************************************)

   (* ---- Scalar sector ---- *)
   (* Sommerfeld conditions everywhere *)
   (* ---- Gauge ---- *)
   rK      == beta[a] dK[a] - alpha sqrt[muL] ( DK + K / r),
   rGams   == sdown[l] beta[b] dG[b,l] - vbetas (DGams + Gams / r ),
   rBs     == sdown[l] beta[b] dB[b,l] - vbetas (DBs + Bs / r ),
   (* ---- Constraint violating conditions ---- *)
   rACss   == beta[a]sup[b]sup[c]dA[a,b,c] - ( DACss + ACss / r) ,
   (* ---- Algebraic constraint preserving condition ---- *)
   rACqq   == chi Ainv[a,b](-2 alpha A[a,b] + lieg[a,b]) - rACss,

   (* ---- Vector sector ---- *)
   (* Sommerfeld conditions everywhere *)
   (* ---- Gauge ---- *)
   rGamA[k] == qud[k,l] beta[b]dG[b,l] - vbetaA (DGamA[k] + GamA[k] / r ),
   rBA[k]   == qud[k,l] beta[b]dB[b,l] - vbetaA (DBA[k]   + BA[k] / r ),
   (* ---- Constraint violating conditions ---- *)
   rACsA[k]  == sup[c]qud[a,k]beta[b]dA[b,a,c] - alpha (DACsA[k] + ACsA[k] / r ) ,

   (* ---- Tensor sector ---- *)
   rACABTF[a,b] == qPhysuudd[c,d,a,b]beta[e]dA[e,c,d] 
                - alpha (DACABTF[a,b] + ACABTF[a,b] / r ),

(************************************************************************)
(* ---- Reconstruct full RHSs from decomposed ones ---- *)
(************************************************************************)

  (* ---- Extrinsic curvature ---- *)
  rA[a,b]   == rACABTF[a,b] + (qud[k,a]sdown[b]+qud[k,b]sdown[a])rACsA[k] 
              + sdown[a]sdown[b] rACss + qdd[a,b] / 2 rACqq,

  (* ---- G and B ---- *)
  rG[a]     == sup[a] rGams + qud[a,k] rGamA[k],
  rB[a]     == sup[a] rBs   + qud[a,k] rBA[k],

  (* check RHS *)
  Cinstruction == "if (CheckForNANandINF(12, rA11,rA12,rA13,rA22,rA23,rA33,
                       rB1,rB2,rB3,rG1,rG2,rK)) {
    printf(\"problem with RHS in bssn_boundary.m\\n\");
    printf(\"x=%2.5e, y=%2.5e, z=%2.5e\\n\",xp[ijk],yp[ijk],zp[ijk]);
  }",
  
  
(************************************************************************)
(* ---- Dissipation ---- *)
(************************************************************************)

  (* add dissipation terms we should test which terms are really needed *)
  Cinstruction == "if (order_dissipation == 4 && boundaryNormore(2)) {",
    rA[a,b] == rA[a,b] + dissfactor dissipation4[A[a,b]],
    rG[a]   == rG[a]   + dissfactor dissipation4[G[a]],
    rK      == rK      + dissfactor dissipation4[K],
    rB[a]   == rB[a]   + dissfactor dissipation4[B[a]],
  Cinstruction == "}",
  Cinstruction == "if (order_dissipation == 6 && boundaryNormore(3)) {",
    rA[a,b] == rA[a,b] + dissfactor dissipation6[A[a,b]],
    rG[a]   == rG[a]   + dissfactor dissipation6[G[a]],
    rK      == rK      + dissfactor dissipation6[K],
    rB[a]   == rB[a]   + dissfactor dissipation6[B[a]],
  Cinstruction == "}",
  Cinstruction == "if (order_dissipation == 8 && boundaryNormore(4)) {",
    rA[a,b] == rA[a,b] + dissfactor dissipation8[A[a,b]],
    rG[a]   == rG[a]   + dissfactor dissipation8[G[a]],
    rK      == rK      + dissfactor dissipation8[K],
    rB[a]   == rB[a]   + dissfactor dissipation8[B[a]],
  Cinstruction == "}",
  Cinstruction == "if (order_dissipation == 10 && boundaryNormore(5)) {",
    rA[a,b] == rA[a,b] + dissfactor dissipation10[A[a,b]],
    rG[a]   == rG[a]   + dissfactor dissipation10[G[a]],
    rK      == rK      + dissfactor dissipation10[K],
    rB[a]   == rB[a]   + dissfactor dissipation10[B[a]],
  Cinstruction == "}",
  Cinstruction == "if (order_dissipation == 12 && boundaryNormore(6)) {",
    rA[a,b] == rA[a,b] + dissfactor dissipation12[A[a,b]],
    rG[a]   == rG[a]   + dissfactor dissipation12[G[a]],
    rK      == rK      + dissfactor dissipation12[K],
    rB[a]   == rB[a]   + dissfactor dissipation12[B[a]],
  Cinstruction == "}",

(************************************************************************)
(* --- The rest of the script is standard BAM setup ... update the RHS---- *)
(************************************************************************)

  Cif == addlinear, 
      nA[a,b]     == pA[a,b]  + c rA[a,b],
      nG[a]       == pG[a]    + c rG[a],
      nK          == pK       + c rK,
      nB[a]       == pB[a]    + c rB[a], 
  Cif == else,
      nA[a,b]     == rA[a,b],
      nG[a]       == rG[a],
      nK          == rK,
      nB[a]       == rB[a],
  Cif == end

}

(************************************************************************)
(* symmetries *)
(************************************************************************)

  g[a_,b_]    :=  g[b,a] /; !OrderedQ[{a,b}]
  pg[a_,b_]   := pg[b,a] /; !OrderedQ[{a,b}]
  ng[a_,b_]   := ng[b,a] /; !OrderedQ[{a,b}]
   A[a_,b_]   :=  A[b,a] /; !OrderedQ[{a,b}]
  pA[a_,b_]   := pA[b,a] /; !OrderedQ[{a,b}]
  nA[a_,b_]   := nA[b,a] /; !OrderedQ[{a,b}]
  Ainv[a_,b_] := Ainv[b,a] /; !OrderedQ[{a,b}]

  lieg[a_,b_] := lieg[b,a] /; !OrderedQ[{a,b}]
  
  dg[c_,a_,b_]    := dg[c,b,a] /; !OrderedQ[{a,b}]
  dA[c_,a_,b_]    := dA[c,b,a] /; !OrderedQ[{a,b}]

  DJac[a_,b_,c_]  := DJac[a,c,b] /; !OrderedQ[{b,c}]

  dgSST[a_,b_,c_]         := dgSST[a,c,b] /; !OrderedQ[{b,c}]
  dASST[a_,b_,c_]         := dASST[a,c,b] /; !OrderedQ[{b,c}]        

  qdd[a_,b_]         := qdd[b,a] /; !OrderedQ[{a,b}]
  quu[a_,b_]         := quu[b,a] /; !OrderedQ[{a,b}]

  qPhysuudd[a_,b_,c_,d_] := qPhysuudd[b,a,c,d]  /; !OrderedQ[{a,b}]
  qPhysuudd[a_,b_,c_,d_] := qPhysuudd[a,b,d,c]  /; !OrderedQ[{c,d}]

  rACABTF[a_,b_]  := rACABTF[b,a] /; !OrderedQ[{a,b}]

  ADMginv[a_,b_]  := ADMginv[b,a] /; !OrderedQ[{a,b}]
  ginv[a_,b_]     := ginv[b,a] /; !OrderedQ[{a,b}]

  rA[a_,b_]       := rA[b,a] /; !OrderedQ[{a,b}]
 
(************************************************************************)
(* information for C output *)
(************************************************************************)

(* information about the function 
   the C file will be functionname.c 
*)
includeANDdefine[] := Module[{},

  pr["#include \"bam.h\"\n"];
  pr["#include \"bssn.h\"\n"];
  pr["\n"];
  pr["#define Power(x,y) (pow((double) (x), (double) (y)))\n"];
  pr["#define Sqrt(x)    sqrt(x)\n"];
  pr["#define Log(x)     log((double) (x))\n"];
  pr["#define pow2(x)    ((x)*(x))\n"];
  pr["#define pow3(x)    ((x)*(x)*(x))\n"];
  pr["#define pow4(x)    ((x)*(x)*(x)*(x))\n"];
  pr["#define pow2inv(x) (1.0/((x)*(x)))\n\n"];
  pr["#define Cal(x,y,z) ((x)?(y):(z))\n\n"];
  pr["#define Tan(x)     tan(x)\n"];
  pr["#define ArcTan(x)  atan(x)\n"];
  pr["#define Sin(x)     sin(x)\n"];
  pr["#define Cos(x)     cos(x)\n"];
  pr["#define Csc(x)     (1./sin(x))\n"];
  pr["#define Abs(x)     (fabs(x))\n"];
  pr["#define sqrt2      (sqrt(2))\n"];
  pr["#define Tanh(x)    tanh(x)\n"];
  pr["#define Sech(x)    (1/cosh(x))\n"];

  pr["\n\n\n"];
];

functionname      = "bssn_boundary";

functionarguments = "tVarList *unew, tVarList *upre, double c, tVarList *ucur";

functionbegin     = "bampi_openmp_start\n";

functionloop      = "forallpoints_ijk_openmp(level) {\n"<>
                    "if (!dequal(rp[ijk],level->bbox[1]-N*level->dx)) continue;\n" <>
                    "if (j==0 || j==box->n-1 || k==0 || k==box->o-1) continue;\n";"forinnerpoints_ijk_openmp(level) {";

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
  
  prdecvlNR[{ g[a,b],  chi,  A[a,b],  K,  G[a],  alpha,  beta[a], B[a]}, "ucur", "METRIC_bssn_INDX_VAR"];
  prdecvlNR[{ng[a,b], nchi, nA[a,b], nK, nG[a], nalpha, nbeta[a],nB[a]}, "unew", "METRIC_bssn_INDX_VAR"];
  prdecvlNR[{pg[a,b], pchi, pA[a,b], pK, pG[a], palpha, pbeta[a],pB[a]}, "upre", "METRIC_bssn_INDX_VAR"];
  pr["\n"];
   
  pr["const int order               = Geti(\"order_centered\");\n"];
  pr["const int N                   = Geti(\"boundary_N_extrapolate\");\n"];
  pr["const int order_dissipation   = Geti(\"order_dissipation\");\n"];
  pr["const double dissfactor       = get_dissipation_factor(level);\n"];
 
  pr["double shiftgammacoeff = Getd(\"bssn_shiftgammacoeff\");\n"];

  (* where are we? *)
  pr["const double *xp = level->v[Ind(\"x\")];\n"];
  pr["const double *yp = level->v[Ind(\"y\")];\n"];
  pr["const double *zp = level->v[Ind(\"z\")];\n"];
  
  (* shells stuff *)
  pr["const double *rp = level->v[IndLax(\"shells_R\")];\n"];
  pr["const double *rr = level->v[IndLax(\"shells_r\")];\n"];
  pr["const double shellsS = GetdLax(\"amr_shells_stretch\");\n"];
  pr["const double shellsR = GetdLax(\"amr_shells_r0\");\n"];
  pr["const double shellsE = GetdLax(\"amr_shells_eps\");\n"];

  pr["\n"];
];




(************************************************************************)
(* now we are ready to go *)


<< "../../math/MathToC/TensorEquationsToC.m"

