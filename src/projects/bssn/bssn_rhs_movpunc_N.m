(* ::Package:: *)

(* bssn_rhs_mov_punc_N.m 
   Bernd Bruegmann 10/98, 10/02
   Wolfgang Tichy  4/2004
   Jose Gonzalez 8/2006 
   Doreen Mueller 4/2010
   mth -> cleanup 12/2011
*)

(* compute right hand side of BSSN equations *)


(* variables *)
variables = { g[a,b],  A[a,b],  G[a],  K,  chi,  alpha,  beta[a],  B[a],
             ng[a,b], nA[a,b], nG[a], nK, nchi, nalpha, nbeta[a], nB[a],
             pg[a,b], pA[a,b], pG[a], pK, pchi, palpha, pbeta[a], pB[a],
             shiftDr,
             xp, yp, zp, rp,rr}



(* compute in this order *)
tocompute = {
    Cinstruction == " if (alpha[ijk]<=0.) { printf(\"ATTENTION:  alpha<=0 at %d\\n\",ijk);",
    Cinstruction == " printf(\"  (%e %e  %e %e)\\n\",alpha[ijk],palpha[ijk],K[ijk],pK[ijk]);",
    Cinstruction == " printf(\"  x=%e y=%e z=%e\\n\",Ptr(level,\"x\")[ijk],Ptr(level,\"y\")[ijk],Ptr(level,\"z\")[ijk]); ",
    Cinstruction == " if (ijkinsidefinerlevel(box,ijk)==0) printf(\"  point is NOT inside finer box NOR in some symmetry area\\n\"); ",
    Cinstruction == " else printf(\"  point is inside finer box/ in symmetry \\n\");} ",

(*---------------------------------------------------------------------------*)
(* partial derivatives, centered differences                                 *)
(*---------------------------------------------------------------------------*)
  Cinstruction == "if (order_centered == 2 || boundaryNaway(1)) {",
    da[a] == del[a,alpha],
    dda[a,b] == deldel[a,b,alpha],
    db[a,b] == del[a,beta[b]],
    ddb[a,b,c] == deldel[a,b,beta[c]],
    delg[c,a,b] == del[c,g[a,b]],
    deldelg[a,b,c,d] == deldel[a,b,g[c,d]],
    delG[a,b] == del[a, G[b]],
    dK[a] == del[a, K],
    dA[a,b,c] == del[a, A[b,c]],
    dchi[a] == del[a,chi],
    ddchi[a,b] == deldel[a,b,chi],
  Cinstruction == "} else if (order_centered == 4 || boundaryNaway(2)) {",
    Cinstruction == "#ifdef REDUCEORDERTO2\n  errorexit(\"CompilerFlag reduced order to 2\");\n #else\n",
    da[a] == del4[a,alpha],
    dda[a,b] == deldel4[a,b,alpha],
    db[a,b] == del4[a,beta[b]],
    ddb[a,b,c] == deldel4[a,b,beta[c]],
    delg[c,a,b] == del4[c,g[a,b]],
    deldelg[a,b,c,d] == deldel4[a,b,g[c,d]],
    delG[a,b] == del4[a, G[b]],
    dK[a] == del4[a, K],
    dA[a,b,c] == del4[a, A[b,c]],
    dchi[a] == del4[a,chi],
    ddchi[a,b] == deldel4[a,b,chi],
    Cinstruction == "#endif\n",
  Cinstruction == "} else if (order_centered == 6 || boundaryNaway(3)) {",
    Cinstruction == "#ifdef REDUCEORDERTO4\n  errorexit(\"CompilerFlag reduced order to 4\");\n #else\n",
    da[a] == del6[a,alpha],
    dda[a,b] == deldel6[a,b,alpha],
    db[a,b] == del6[a,beta[b]],
    ddb[a,b,c] == deldel6[a,b,beta[c]],
    delg[c,a,b] == del6[c,g[a,b]],
    deldelg[a,b,c,d] == deldel6[a,b,g[c,d]],
    delG[a,b] == del6[a, G[b]],
    dK[a] == del6[a, K],
    dA[a,b,c] == del6[a, A[b,c]],
    dchi[a] == del6[a,chi],
    ddchi[a,b] == deldel6[a,b,chi],
    Cinstruction == "#endif\n",
  Cinstruction == "} else if (order_centered == 8 || boundaryNaway(4)) {",
    Cinstruction == "#ifdef REDUCEORDERTO6\n  errorexit(\"CompilerFlag reduced order to 6\");\n #else\n",
    da[a] == del8[a,alpha],
    dda[a,b] == deldel8[a,b,alpha],
    db[a,b] == del8[a,beta[b]],
    ddb[a,b,c] == deldel8[a,b,beta[c]],
    delg[c,a,b] == del8[c,g[a,b]],
    deldelg[a,b,c,d] == deldel8[a,b,g[c,d]],
    delG[a,b] == del8[a, G[b]],
    dK[a] == del8[a, K],
    dA[a,b,c] == del8[a, A[b,c]],
    dchi[a] == del8[a,chi],
    ddchi[a,b] == deldel8[a,b,chi],
    Cinstruction == "#endif\n",
  Cinstruction == "} else {",
    Cinstruction == "#ifdef REDUCEORDERTO8\n  errorexit(\"CompilerFlag reduced order to 8\");\n #else\n",
    da[a] == del10[a,alpha],
    dda[a,b] == deldel10[a,b,alpha],
    db[a,b] == del10[a,beta[b]],
    ddb[a,b,c] == deldel10[a,b,beta[c]],
    delg[c,a,b] == del10[c,g[a,b]],
    deldelg[a,b,c,d] == deldel10[a,b,g[c,d]],
    delG[a,b] == del10[a, G[b]],
    dK[a] == del10[a, K],
    dA[a,b,c] == del10[a, A[b,c]],
    dchi[a] == del10[a,chi],
    ddchi[a,b] == deldel10[a,b,chi],
    Cinstruction == "#endif\n",
  Cinstruction == "}",



(*---------------------------------------------------------------------------*)
(* advection derivatives, treated differently than centered derivatives      *)
(* typically we store the advection derivative under lieSomething and        *)
(* add the missing terms later                                               *)
(*---------------------------------------------------------------------------*)
  w1 == beta1/dx, w2 == beta2/dy, w3 == beta3/dz,
  Cinstruction == "if (order_advection == 2 || boundaryNaway(1)) {",
    Cinstruction == "set_advection2(w1, w2, w3);",
    lieg[a,b] == adv[g[a,b]],
    lieA[a,b] == adv[A[a,b]],
    lieK == adv[K],
    liechi == adv[chi],
    liealpha == adv[alpha],
    advG[a] == adv[G[a]],
    advbeta[a] == adv[beta[a]],
    advB[a] == adv[B[a]],
  (* legacy: don't drop below order 4, should experiment *)
  Cinstruction == "} else if (order_advection == 4 || boundaryNaway(2)) {",
    Cinstruction == "#ifdef REDUCEORDERTO2\n  errorexit(\"CompilerFlag reduced order to 2\");\n #else\n",
    Cinstruction == "set_advection4(w1, w2, w3);",
    lieg[a,b] == adv4[g[a,b]],
    lieA[a,b] == adv4[A[a,b]],
    lieK == adv4[K],
    liechi == adv4[chi],
    liealpha == adv4[alpha],
    advG[a] == adv4[G[a]],
    advbeta[a] == adv4[beta[a]],
    advB[a] == adv4[B[a]],
    Cinstruction == "#endif\n",
  Cinstruction == "} else if (order_advection == 6 || boundaryNaway(3)) {",
    Cinstruction == "#ifdef REDUCEORDERTO4\n  errorexit(\"CompilerFlag reduced order to 4\");\n #else\n",
    Cinstruction == "set_advection6(w1, w2, w3);",
    lieg[a,b] == adv6[g[a,b]],
    lieA[a,b] == adv6[A[a,b]],
    lieK == adv6[K],
    liechi == adv6[chi],
    liealpha == adv6[alpha],
    advG[a] == adv6[G[a]],
    advbeta[a] == adv6[beta[a]],
    advB[a] == adv6[B[a]],
    Cinstruction == "#endif\n",
  Cinstruction == "} else if (order_advection == 8 || boundaryNaway(4)) {",
   Cinstruction == "#ifdef REDUCEORDERTO6\n  errorexit(\"CompilerFlag reduced order to 6\");\n #else\n",
    Cinstruction == "set_advection8(w1, w2, w3);",
    lieg[a,b] == adv8[g[a,b]],
    lieA[a,b] == adv8[A[a,b]],
    lieK == adv8[K],
    liechi == adv8[chi],
    liealpha == adv8[alpha],
    advG[a] == adv8[G[a]],
    advbeta[a] == adv8[beta[a]],
    advB[a] == adv8[B[a]],
    Cinstruction == "#endif\n",
  Cinstruction == "} else if (order_advection == 10 || boundaryNaway(5)) {",
    Cinstruction == "#ifdef REDUCEORDERTO8\n  errorexit(\"CompilerFlag reduced order to 8\");\n #else\n",
    Cinstruction == "set_advection10(w1, w2, w3);",
    lieg[a,b] == adv10[g[a,b]],
    lieA[a,b] == adv10[A[a,b]],
    lieK == adv10[K],
    liechi == adv10[chi],
    liealpha == adv10[alpha],
    advG[a] == adv10[G[a]],
    advbeta[a] == adv10[beta[a]],
    advB[a] == adv10[B[a]],
    Cinstruction == "#endif\n",
  Cinstruction == "} else errorexit(\"order not implemented\"); ",


(*---------------------------------------------------------------------------*)
(* do coordinatetransfo if we use shells                                     *)
(*---------------------------------------------------------------------------*)
  Cinstruction == "if (useShellsTransfo) {",
    
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
    
    Cinstruction == "if (order_advection == 2 || boundaryNaway(1)) {",
      dB[a,b] == del[a,B[b]],
    Cinstruction == "} else if (order_advection == 4 || boundaryNaway(2)) {",
      dB[a,b] == del4[a,B[b]],
    Cinstruction == "} else if (order_advection == 6 || boundaryNaway(3)) {",
      dB[a,b] == del6[a,B[b]],
    Cinstruction == "} else if (order_centered == 8 || boundaryNaway(4)) {",
      dB[a,b] == del8[a,B[b]],
    Cinstruction == "} else {",
      dB[a,b] == del10[a,B[b]],
    Cinstruction == "}",
    
    (* now do the transformation and store in a tmp name *)
    daSST[a]        == Jac[b,a] da[b],
    ddaSST[a,b]     == DJac[c,a,b] da[c] + Jac[c,a] Jac[d,b] dda[c,d],
    dbSST[a,b]      == Jac[c,a] db[c,b],
    ddbSST[a,b,c]   == DJac[d,a,b] db[d,c] + Jac[d,a] Jac[e,b] ddb[d,e,c],
    delgSST[a,b,c]  == Jac[d,a] delg[d,b,c],
    deldelgSST[a,b,c,d] == DJac[e,a,b] delg[e,c,d] + Jac[e,a] Jac[f,b] deldelg[e,f,c,d],
    delGSST[a,b]    == Jac[c,a] delG[c,b],
    dKSST[a]        == Jac[b,a] dK[b],
    dASST[a,b,c]    == Jac[d,a] dA[d,b,c],
    dchiSST[a]      == Jac[b,a] dchi[b],
    ddchiSST[a,b]   == DJac[c,a,b] dchi[c] + Jac[c,a] Jac[d,b] ddchi[c,d],
    dBSST[a,b]      == Jac[c,a] dB[c,b],
    
    (* now give the tmp derivative the same name again *)
    da[a]         == daSST[a],
    dda[a,b]      == ddaSST[a,b],
    db[a,b]       == dbSST[a,b],
    ddb[a,b,c]    == ddbSST[a,b,c],
    delg[a,b,c]   == delgSST[a,b,c],
    deldelg[a,b,c,d]  == deldelgSST[a,b,c,d],
    delG[a,b]     == delGSST[a,b],
    dK[a]         == dKSST[a],
    dA[a,b,c]     == dASST[a,b,c],
    dchi[a]       == dchiSST[a],
    ddchi[a,b]    == ddchiSST[a,b],
    dB[a,b]       == dBSST[a,b],
    
    (* now advection derivatives, FIXME: do it better or prove that this is enough at level 0 *)
    lieg[a,b]     == beta[c] delg[c,a,b],
    lieA[a,b]     == beta[c] dA[c,a,b],
    lieK          == beta[a] dK[a],
    liechi        == beta[a] dchi[a],
    liealpha      == beta[a] da[a],
    advG[a]       == beta[b] delG[b,a],
    advbeta[a]    == beta[b] db[b,a],
    advB[a]       == beta[b] dB[b,a],
    
  Cinstruction == "}",




  (* inverse conformal metric *)
  detg == matrixdet[g],
  Cinstruction == "if (detg<=0.01) { PROBLEM; }",
  detginv == 1/detg,
  ginv[a,b] == detginv matrixinvdet[g,a,b],

  (* derivatives of conformal metric *)
  (* dginv[a,b,c] == - ginv[b,d] ginv[c,e] delg[a,d,e], *)

  (* connection of conformal metric *) 
  gammado[c,a,b] == 1/2 (delg[a,b,c] + delg[b,a,c] - delg[c,a,b]),
  gamma[c,a,b] == ginv[c,d] gammado[d,a,b], 
  Gfromg[a] == ginv[b,c] gamma[a,b,c],

  (* curvature of conformal metric *)
  R[a,b] == ginv[c,d] ( -1/2 deldelg[c,d,a,b] +
                    gamma[e,c,a] gammado[b,e,d] +
                    gamma[e,c,b] gammado[a,e,d] +
                    gamma[e,a,d] gammado[e,c,b]) +
               1/2 (g[c,a] delG[b,c] + g[c,b] delG[a,c] +
               Gfromg[c] (gammado[a,b,c] + gammado[b,a,c])),


  (* conformal factor and its derivatives *)
  
  chiguard      == chiDivFloor,
  chiguarded    == chi,
  Cinstruction  == "if (chiguarded < chiguard) chiguarded = chiguard;",
  oochipsipower == 1/chipsipower,
  f             == oochipsipower log[chiguarded],
  df[a]         == oochipsipower dchi[a] / chiguarded,
  ddf[a,b]      == oochipsipower ddchi[a,b] / chiguarded  -  chipsipower df[a] df[b],
  
  psim2         == exp[-2 f],
  psim4         == exp[-4 f],

  cddf[a,b]     == ddf[a,b] - gamma[c,a,b] df[c],
  trcddf        == ginv[a,b] cddf[a,b],

  (* curvature contribution from conformal factor *)
  Rphi[a,b]     == 4 df[a] df[b] - 2 cddf[a,b] - 2 g[a,b] trcddf -
                   4 g[a,b] ginv[c,d] df[c] df[d],

  (* derivatives of the lapse *)
  cdda[a,b]     == dda[a,b] - gamma[c,a,b] da[c] -
                   2 (df[a] da[b] + df[b] da[a] - g[a,b] ginv[c,d] df[c] da[d]),
  trcdda        == psim4 (ginv[a,b] cdda[a,b]),

  (* set K zero *)
  K             == K forceKzerofactor,
  dK[a]         == dK[a] forceKzerofactor,

  (* tf part of conformal extrinsic curvature *)
  AA[a,b]       == ginv[c,d] A[a,c] A[d,b],
  AA            == ginv[a,b] AA[a,b],
  Ainv[a,b]     == ginv[a,c] ginv[b,d] A[c,d],
  divAinv[a]    == 2/3 ginv[a,b] dK[b] - 6 Ainv[a,b] df[b] -
                   gamma[a,b,c] Ainv[b,c],

  (* Ricci scalar *)
  R             == AA - 2 K K / 3,
  
  
  (* the shift terms *)
  divbeta       == delta[a,b] db[a,b],
  totdivbeta    == 2/3 divbeta,
  ootddivbeta[a] == 1/3 delta[b,c] ddb[a,b,c],
  
  lieg[a,b]     == lieg[a,b] - g[a,b] totdivbeta + db[a,c] g[b,c] + db[b,c] g[a,c],
  lieA[a,b]     == lieA[a,b] - A[a,b] totdivbeta + db[a,c] A[b,c] + db[b,c] A[a,c],
  liechi        == liechi + (chipsipower/6) chiguarded divbeta,

  pseudolieG[a] == ginv[b,c] ddb[b,c,a] + ginv[a,b] ootddivbeta[b] -
                   Gfromg[b] db[b,a] + Gfromg[a] totdivbeta + 
                   advG[a],



(* right hand sides *)
  rg[a,b]       == - 2 alpha A[a,b] + lieg[a,b],

  rA[a,b]       == psim4 ( - cdda[a,b] + alpha (R[a,b] + Rphi[a,b])) -
                   1/3 g[a,b] (- trcdda + alpha R) +
                   alpha (K A[a,b] - 2 AA[a,b]) + lieA[a,b],

  rG[a]         == - 2 Ainv[a,b] da[b] - 2 alpha divAinv[a] + pseudolieG[a], 

  rK            == - trcdda + alpha (AA + K K/3) + lieK, 

  rchi          == (-chipsipower/6) chiguarded alpha K + liechi,


(*---------------------------------------------------------------------------*)
(* gauge conditions, avoid unnecessary if statements, use flags *)
(*---------------------------------------------------------------------------*)
  falpha        == oploglapse lapseharmonicf + oploglapse2 8/3/(3-alpha) +
                   harmoniclapse alpha,
  ralpha        == oplogwithshift liealpha -
                   falpha alpha K,


  betaF         == shiftgammacoeff alpha^shiftalphapower,

  rbeta[a]      == (gamma2factor + gamma0factor betaF) withB B[a] +  withShiftadv advbeta[a] + 
                   (1-withB) (shiftgammacoeff alpha^shiftGalphapower G[a] - shiftdriver beta[a]),



(*---------------------------------------------------------------------------*)
(* location-dependent shiftdriver                                            *)
(*---------------------------------------------------------------------------*)
  Cinstruction == "if (use_eta) {",
    
    dpsim2[a]   == -2*psim2*df[a],
    absdpsim2   == Sqrt[ginv[a,b]dpsim2[a]dpsim2[b]],
    
    Cinstruction == "shiftDr[ijk] = bssn_eta_set(xp[ijk],yp[ijk],zp[ijk], absdpsim2,psim2);",

    rB[a]       == (gamma0factor + gamma2factor betaF) * (rG[a] - withGadv advG[a]) 
                   - shiftDr B[a] + withBadv advB[a],

  Cinstruction == "} else {", (* not useeta*)

    rB[a]       == withB * ((gamma0factor + gamma2factor betaF) *
                   (rG[a] - withGadv advG[a]) - shiftdriver B[a] + withBadv advB[a]),

  Cinstruction == "}",


(*---------------------------------------------------------------------------*)
(* add dissipation terms we should test which terms are really needed        *)
(*---------------------------------------------------------------------------*)
  Cinstruction == "if (order_dissipation == 4 && boundaryNormore(2)) {",
    rg[a,b] == rg[a,b] + dissfactor dissipation4[g[a,b]],
    rA[a,b] == rA[a,b] + dissfactor dissipation4[A[a,b]],
    rG[a]   == rG[a]   + dissfactor dissipation4[G[a]],
    rK      == rK      + dissfactor dissipation4[K],
    rchi    == rchi    + dissfactor dissipation4[chi],
  Cinstruction == "}",
  Cinstruction == "if (order_dissipation == 6 && boundaryNormore(3)) {",
    Cinstruction == "#ifdef REDUCEORDERTO2\n  errorexit(\"CompilerFlag reduced order to 2\");\n #else\n",
    rg[a,b] == rg[a,b] + dissfactor dissipation6[g[a,b]],
    rA[a,b] == rA[a,b] + dissfactor dissipation6[A[a,b]],
    rG[a]   == rG[a]   + dissfactor dissipation6[G[a]],
    rK      == rK      + dissfactor dissipation6[K],
    rchi    == rchi    + dissfactor dissipation6[chi],
    Cinstruction == "#endif\n",
  Cinstruction == "}",
  Cinstruction == "if (order_dissipation == 8 && boundaryNormore(4)) {",
    Cinstruction == "#ifdef REDUCEORDERTO4\n  errorexit(\"CompilerFlag reduced order to 4\");\n #else\n",
    rg[a,b] == rg[a,b] + dissfactor dissipation8[g[a,b]],
    rA[a,b] == rA[a,b] + dissfactor dissipation8[A[a,b]],
    rG[a]   == rG[a]   + dissfactor dissipation8[G[a]],
    rK      == rK      + dissfactor dissipation8[K],
    rchi    == rchi    + dissfactor dissipation8[chi],
    Cinstruction == "#endif\n",
  Cinstruction == "}",
  Cinstruction == "if (order_dissipation == 10 && boundaryNormore(5)) {",
    Cinstruction == "#ifdef REDUCEORDERTO6\n  errorexit(\"CompilerFlag reduced order to 6\");\n #else\n",
    rg[a,b] == rg[a,b] + dissfactor dissipation10[g[a,b]],
    rA[a,b] == rA[a,b] + dissfactor dissipation10[A[a,b]],
    rG[a]   == rG[a]   + dissfactor dissipation10[G[a]],
    rK      == rK      + dissfactor dissipation10[K],
    rchi    == rchi    + dissfactor dissipation10[chi],
    Cinstruction == "#endif\n",
  Cinstruction == "}",
  Cinstruction == "if (order_dissipation == 12 && boundaryNormore(6)) {",
    Cinstruction == "#ifdef REDUCEORDERTO8\n  errorexit(\"CompilerFlag reduced order to 8\");\n #else\n",
    rg[a,b] == rg[a,b] + dissfactor dissipation12[g[a,b]],
    rA[a,b] == rA[a,b] + dissfactor dissipation12[A[a,b]],
    rG[a]   == rG[a]   + dissfactor dissipation12[G[a]],
    rK      == rK      + dissfactor dissipation12[K],
    rchi    == rchi    + dissfactor dissipation12[chi],
    Cinstruction == "#endif\n",
  Cinstruction == "}",



(*---------------------------------------------------------------------------*)
(* test RHS                                                                  *)
(*---------------------------------------------------------------------------*)
  Cinstruction == "if (setRHSto0) {",
    rg[a,b] == 0,
    rA[a,b] == 0,
    rG[a]   == 0,
    rK      == 0,
    rchi    == 0,
    rbeta[a]== 0,
    rB[a]   == 0,
    ralpha  == 0,
  Cinstruction == "}",
  Cinstruction == "CheckForNANandINF(24, rA11,rA12,rA13,rA22,rA23,rA33, ralpha,rB1,rB2,rB3,rbeta1,rbeta2,rbeta3,rchi, rG1,rg11,rg12,rg13,rG2,rg22,rg23,rG3,g33,rK);",
  


(*---------------------------------------------------------------------------*)
(* result ... add RHS                                                        *)
(*---------------------------------------------------------------------------*)
  Cif == addlinear, 
    ng[a,b]   == pg[a,b] + c rg[a,b],
    nA[a,b]   == pA[a,b] + c rA[a,b],
    nG[a]     == pG[a] + c rG[a],
    nK        == forceKzerofactor (pK + c rK),
    nchi      == pchi + c rchi,
    nalpha    == palpha + c ralpha,
    nbeta[a]  == pbeta[a] + c rbeta[a],
    nB[a]     == pB[a] + c rB[a], 

    (* normalize detg and subtract trace from As *)
    detnginv == 1/matrixdet[ng],
    Cif == subtractA,
      traceA  == detnginv matrixinvdet[ng,a,b] nA[a,b],
      aux     == -traceA/3,
      nA[a,b] == nA[a,b] + aux ng[a,b],
    Cif == end,
    Cif == normalizedetg,
      aux     == detnginv^(1/3),
      ng[a,b] == aux ng[a,b], 
    Cif == end,

  Cif == else,
    ng[a,b]   == rg[a,b],
    nchi      == rchi,
    nA[a,b]   == rA[a,b],
    nG[a]     == rG[a],
    nK        == forceKzerofactor rK,
    nalpha    == ralpha,
    nbeta[a]  == rbeta[a],
    nB[a]     == rB[a],
  Cif == end

}


(* symmetries *)
 g[a_,b_] :=  g[b,a] /; !OrderedQ[{a,b}]
pg[a_,b_] := pg[b,a] /; !OrderedQ[{a,b}]
ng[a_,b_] := ng[b,a] /; !OrderedQ[{a,b}]
rg[a_,b_] := rg[b,a] /; !OrderedQ[{a,b}]

 A[a_,b_] :=  A[b,a] /; !OrderedQ[{a,b}]
pA[a_,b_] := pA[b,a] /; !OrderedQ[{a,b}]
nA[a_,b_] := nA[b,a] /; !OrderedQ[{a,b}]
rA[a_,b_] := rA[b,a] /; !OrderedQ[{a,b}]

ginv[a_,b_] := ginv[b,a] /; !OrderedQ[{a,b}]
Ainv[a_,b_] := Ainv[b,a] /; !OrderedQ[{a,b}]
Kinv[a_,b_] := Kinv[b,a] /; !OrderedQ[{a,b}]
K[a_,b_] := K[b,a] /; !OrderedQ[{a,b}]
AA[a_,b_] := AA[b,a] /; !OrderedQ[{a,b}]
R[a_,b_] := R[b,a] /; !OrderedQ[{a,b}]
Rphi[a_,b_] := Rphi[b,a] /; !OrderedQ[{a,b}]
lieg[a_,b_] := lieg[b,a] /; !OrderedQ[{a,b}]
lieA[a_,b_] := lieA[b,a] /; !OrderedQ[{a,b}]
result[a_,b_] := result[b,a] /; !OrderedQ[{a,b}]

ddchi[a_,b_] := ddchi[b,a] /; !OrderedQ[{a,b}]
ddf[a_,b_] := ddf[b,a] /; !OrderedQ[{a,b}]
cddf[a_,b_] := cddf[b,a] /; !OrderedQ[{a,b}]
dda[a_,b_] := dda[b,a] /; !OrderedQ[{a,b}]
cdda[a_,b_] := cdda[b,a] /; !OrderedQ[{a,b}]
ddpop[a_,b_] := ddpop[b,a] /; !OrderedQ[{a,b}]
ddb[a_,b_,c_] := ddb[b,a,c] /; !OrderedQ[{a,b}]

deldel[a_,b_,c_] := deldel[b,a,c] /; !OrderedQ[{a,b}]
deldel4[a_,b_,c_] := deldel4[b,a,c] /; !OrderedQ[{a,b}]
deldel6[a_,b_,c_] := deldel6[b,a,c] /; !OrderedQ[{a,b}]

delg[c_,a_,b_] := delg[c,b,a] /; !OrderedQ[{a,b}]
delgb[c_,a_,b_] := delgb[c,b,a] /; !OrderedQ[{a,b}]
dginv[c_,a_,b_] := dginv[c,b,a] /; !OrderedQ[{a,b}]
dA[c_,a_,b_] := dA[c,b,a] /; !OrderedQ[{a,b}]

deldelg[a_,b_,c_,d_] := deldelg[b,a,c,d] /; !OrderedQ[{a,b}]
deldelg[a_,b_,c_,d_] := deldelg[a,b,d,c] /; !OrderedQ[{c,d}]

codel[a_, (x_ /; NumberQ[x])] = 0
codelK[c_,a_,b_] := codelK[c,b,a] /; !OrderedQ[{a,b}]

ddfSST[a_,b_]     := ddfSST[b,a] /; !OrderedQ[{a,b}]
ddchiSST[a_,b_]   := ddchiSST[b,a] /; !OrderedQ[{a,b}]
ddaSST[a_,b_]     := ddaSST[b,a] /; !OrderedQ[{a,b}]
ddbSST[a_,b_,c_]  := ddbSST[b,a,c] /; !OrderedQ[{a,b}]
delgSST[c_,a_,b_] := delgSST[c,b,a] /; !OrderedQ[{a,b}]
dASST[c_,a_,b_]   := dASST[c,b,a] /; !OrderedQ[{a,b}]
deldelgSST[a_,b_,c_,d_] := deldelgSST[b,a,c,d] /; !OrderedQ[{a,b}]
deldelgSST[a_,b_,c_,d_] := deldelgSST[a,b,d,c] /; !OrderedQ[{c,d}]

DJac[a_,b_,c_]        := DJac[a,c,b] /; !OrderedQ[{b,c}]




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
  pr["#define Sqrt(x)    sqrt(x)\n"];
  pr["#define Log(x)     log((double) (x))\n"];
  pr["#define pow2(x)    ((x)*(x))\n"];
  pr["#define pow3(x)    ((x)*(x)*(x))\n"];
  pr["#define pow4(x)    ((x)*(x)*(x)*(x))\n"];
  pr["#define pow2inv(x) (1.0/((x)*(x)))\n"];
  pr["#define Tanh(x)    tanh(x)\n"];
  pr["#define Sech(x)    (1/cosh(x))\n"];
  pr["#define Cal(x,y,z) ((x)?(y):(z))\n"];
  pr["#define PROBLEM    printf(\"  %d pts away from boundary\\n\",boundaryaway(6)); \\\n"];
  pr["                   printf(\"    %e %e %e \\n    detg=%e \\n\", \\\n"];
  pr["                   Ptr(level,\"x\")[ijk],Ptr(level,\"y\")[ijk],Ptr(level,\"z\")[ijk],detg); \\\n"];
  pr["                   if (ijkinsidefinerlevel(box,ijk)==0) printf(\"    point is NOT inside finer box NOR in some symmetry area\\n\"); \\\n"];
  pr["                   else printf(\"    point is inside finer box/ in symmetry \\n\"); \\\n"];
  pr["                   printf(\"    gtilde = %e %e %e %e %e %e \\n\"  ,g11[ijk],g12[ijk],g13[ijk],g22[ijk], g23[ijk],g33[ijk]);\\\n"];
  pr["\n\n\n"];
];

functionname      = "bssn_rhs_movpunc_N";

functionarguments = "tVarList *unew, tVarList *upre, double c, tVarList *ucur";

functionbegin     = "int_advectionN_stencil;\n" <>
                    "bampi_openmp_start\n" <>
                    "int_advectionN_vars;\n";

functionloop      = "forinnerpoints_ijk_openmp(level) {";

endfunctionloop   = "} endfor_ijk_openmp; /* loop i, j, k */";

endfunction       = "bampi_openmp_stop\n";




(* custom variable declaration
   we have to translate between the tensor names and the C variables  
   we may need special variables, for example for parameters
*)
variabledeclarations[] := Module[{},

  pr["tL *level = ucur->level;\n"];
  pr["\n"];
  
  prdecvlNR[{ g[a,b],  chi,  A[a,b],  K,  G[a],  alpha,  beta[a],  B[a]}, "ucur", "METRIC_bssn_INDX_VAR"];
  prdecvlNR[{ng[a,b], nchi, nA[a,b], nK, nG[a], nalpha, nbeta[a], nB[a]}, "unew", "METRIC_bssn_INDX_VAR"];
  prdecvlNR[{pg[a,b], pchi, pA[a,b], pK, pG[a], palpha, pbeta[a], pB[a]}, "upre", "METRIC_bssn_INDX_VAR"];
  pr["\n"];
  
  prdecvarname[{shiftDr}, "bssn_eta"];
  pr["\n"];
 
  pr["const int addlinear           = (c != 0.0l);\n"];
  
  pr["const int order_centered      = Geti(\"order_centered\");\n"];
  pr["const int order_advection     = Geti(\"order_advection\");\n"];
  pr["const int advectionlopsided   = Geti(\"advection_lopsided\");\n"];
  pr["const int advectionlopsided6  = Geti(\"advection_lopsided6\");\n"];
  pr["const int advectionlopsided8  = Geti(\"advection_lopsided8\");\n"];
  pr["const int advectionlopsided10 = Geti(\"advection_lopsided10\");\n"];

  pr["const int order_dissipation = Geti(\"order_dissipation\");\n"];
  pr["double dissfactor = get_dissipation_factor(level);\n"];

  pr["const double forceKzerofactor = Getv(\"bssn_forceKzero\", \"no\");\n"];
  pr["const int subtractA      = Getv(\"bssn_subtractA\", \"yes\");\n"];
  pr["const int normalizedetg  = Getv(\"bssn_normalizedetg\", \"yes\");\n"];
  pr["const int constantlapse  = Getv(\"bssn_lapse\", \"constant\");\n"];
  pr["const int oploglapse     = Getv(\"bssn_lapse\", \"1+log\");\n"];
  pr["const int oploglapse2    = Getv(\"bssn_lapse\", \"1+log2\");\n"];
  pr["const int oplogwithshift = Getv(\"bssn_lapse\", \"withshift\");\n"];
  pr["const int harmoniclapse  = Getv(\"bssn_lapse\", \"harmonic\");\n"];

  pr["const double gamma0factor    = Getv(\"bssn_shift\", \"gamma0\");\n"];
  pr["const double gamma2factor    = Getv(\"bssn_shift\", \"gamma2\");\n"];
  pr["const double withGadv        = Getv(\"bssn_shift\", \"withGadv\");\n"];
  pr["const double withShiftadv    = Getv(\"bssn_shift\", \"withShiftadv\");\n"];
  pr["const double withBadv        = Getv(\"bssn_shift\", \"withBadv\");\n"];
  pr["const double withB           = !Getv(\"bssn_shift\", \"withoutB\");\n"];
 
  pr["const double lapseharmonicf  = Getd(\"bssn_lapseharmonicf\");\n"];
  pr["const double shiftalphapower = Getd(\"bssn_shiftalphapower\");\n"];
  pr["const double shiftgammacoeff = Getd(\"bssn_shiftgammacoeff\");\n"];
  pr["const double shiftGalphapower= Getd(\"bssn_shiftgammaalphapower\");\n"];
  pr["const double shiftdriver     = Getd(\"bssn_shiftdriver\");\n"];

  pr["const double chiDivFloor     = Getd(\"bssn_chi_div_floor\");\n"];
  pr["const double chipsipower     = Getd(\"bssn_chi_psipower\");\n"];
  
  
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
  
  
  pr["const int setRHSto0 = Getv(\"bssnRHSto0\",\"yes\");\n"];
  pr["\n"];
  
  
];




(************************************************************************)
(* now we are ready to go *)



<< "../../math/MathToC/TensorEquationsToC.m"

