(* Z4d_rhs_mov_punc_6.m 
   Wolfgang Tichy 04/2004 
   SB,DH 12/2015 add linear extrapolation to PHYSBOUND for certain BCs *)

(* compute right hand side of Z4d equations *)


(* variables *)
variables = {g[a,b],  A[a,b],  G[a],  Khat,  chi,   Theta,  alpha,  beta[a],  B[a], 
             ng[a,b], nA[a,b], nG[a], nKhat, nchi, nTheta, nalpha, nbeta[a], nB[a],
             pg[a,b], pA[a,b], pG[a], pKhat, pchi, pTheta, palpha, pbeta[a], pB[a],
             Zinv[a],
             xp, yp, zp, rr,rp}


(* compute in this order *)
tocompute = {

  Cinstruction == "CheckForNANandINF(25, 
    g11[ijk],g12[ijk],g13[ijk],g22[ijk],g23[ijk],g33[ijk],
    A11[ijk],A12[ijk],A13[ijk],A22[ijk],A23[ijk],A33[ijk], 
    G1[ijk],G2[ijk],G3[ijk], Khat[ijk],chi[ijk], Theta[ijk],
    alpha[ijk],beta1[ijk],beta2[ijk],beta3[ijk],B1[ijk],B2[ijk],B3[ijk]);", 

  (* partial derivatives, centered differences *)
  Cinstruction == "if (order_centered == 2 || boundary1away) {",
    da[a] == del[a,alpha],
    dda[a,b] == deldel[a,b,alpha],
    db[a,b] == del[a,beta[b]],
    ddb[a,b,c] == deldel[a,b,beta[c]],
    dB[a,b] == del[a,B[b]],
    delg[c,a,b] == del[c,g[a,b]],
    deldelg[a,b,c,d] == deldel[a,b,g[c,d]],
    delG[a,b] == del[a, G[b]],
    dKhat[a] == del[a, Khat],
    dA[a,b,c] == del[a, A[b,c]],
    dchi[a] == del[a,chi],
    ddchi[a,b] == deldel[a,b,chi],
    dTheta[a] == del[a,Theta],
  Cinstruction == "} else if (order_centered == 4 || boundaryNaway(2))  {",
    da[a] == del4[a,alpha],
    dda[a,b] == deldel4[a,b,alpha],
    db[a,b] == del4[a,beta[b]],
    ddb[a,b,c] == deldel4[a,b,beta[c]],
    dB[a,b] == del4[a,B[b]],
    delg[c,a,b] == del4[c,g[a,b]],
    deldelg[a,b,c,d] == deldel4[a,b,g[c,d]],
    delG[a,b] == del4[a, G[b]],
    dKhat[a] == del4[a, Khat],
    dA[a,b,c] == del4[a, A[b,c]],
    dchi[a] == del4[a,chi],
    ddchi[a,b] == deldel4[a,b,chi],
    dTheta[a] == del4[a,Theta],
  Cinstruction == "} else if (order_centered == 6 || boundaryNaway(3))  {",
    da[a] == del6[a,alpha],
    dda[a,b] == deldel6[a,b,alpha],
    db[a,b] == del6[a,beta[b]],
    ddb[a,b,c] == deldel6[a,b,beta[c]],
    dB[a,b] == del6[a,B[b]],
    delg[c,a,b] == del6[c,g[a,b]],
    deldelg[a,b,c,d] == deldel6[a,b,g[c,d]],
    delG[a,b] == del6[a, G[b]],
    dKhat[a] == del6[a, Khat],
    dA[a,b,c] == del6[a, A[b,c]],
    dchi[a] == del6[a,chi],
    ddchi[a,b] == deldel6[a,b,chi],
    dTheta[a] == del6[a,Theta],
  Cinstruction == "} else errorexit(\"this order is not yet implemented\"); ",

  (* for faster counting: 1 10 19 37 55 91 100 103 121 130 *)
  Cinstruction == "CheckForNANandINF(132, 
    da1,da2,da3, dda11,dda12,dda13,dda22,dda23,dda33,
    db11,db12,db13,db21,db22,db23,db31,db32,db33,
    ddb111,ddb121,ddb131,ddb221,ddb231,ddb331,
    ddb112,ddb122,ddb132,ddb222,ddb232,ddb332,
    ddb113,ddb123,ddb133,ddb223,ddb233,ddb333,
    delg111,delg112,delg113,delg122,delg123,delg133,
    delg211,delg212,delg213,delg222,delg223,delg233,
    delg311,delg312,delg313,delg322,delg323,delg333,
    deldelg1111,deldelg1112,deldelg1113,deldelg1122,deldelg1123,deldelg1133,
    deldelg1211,deldelg1212,deldelg1213,deldelg1222,deldelg1223,deldelg1233,
    deldelg1311,deldelg1312,deldelg1313,deldelg1322,deldelg1323,deldelg1333,
    deldelg2211,deldelg2212,deldelg2213,deldelg2222,deldelg2223,deldelg2233,
    deldelg2311,deldelg2312,deldelg2313,deldelg2322,deldelg2323,deldelg2333,
    deldelg3311,deldelg3312,deldelg3313,deldelg3322,deldelg3323,deldelg3333,
    delG11,delG12,delG13,delG21,delG22,delG23,delG31,delG32,delG33,
    dKhat1,dKhat2,dKhat3,
    dA111,dA112,dA113,dA122,dA123,dA133,
    dA211,dA212,dA213,dA222,dA223,dA233,
    dA311,dA312,dA313,dA322,dA323,dA333,
    dchi1,dchi2,dchi3,ddchi11,ddchi12,ddchi13,ddchi22,ddchi23,ddchi33,
    dTheta1,dTheta2,dTheta3);",





  (* do coordinatetransfo if we use shells + compute adv derivs in a simple way *)
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
    
    (* now do the transformation and store in a tmp name *)
    daSST[a]        == Jac[b,a] da[b],
    ddaSST[a,b]     == DJac[c,a,b] da[c] + Jac[c,a] Jac[d,b] dda[c,d],
    dbSST[a,b]      == Jac[c,a] db[c,b],
    ddbSST[a,b,c]   == DJac[d,a,b] db[d,c] + Jac[d,a] Jac[e,b] ddb[d,e,c],
    delgSST[a,b,c]  == Jac[d,a] delg[d,b,c],
    deldelgSST[a,b,c,d] == DJac[e,a,b] delg[e,c,d] + Jac[e,a] Jac[f,b] deldelg[e,f,c,d],
    delGSST[a,b]    == Jac[c,a] delG[c,b],
    dKhatSST[a]     == Jac[b,a] dKhat[b],
    dASST[a,b,c]    == Jac[d,a] dA[d,b,c],
    dchiSST[a]      == Jac[b,a] dchi[b],
    ddchiSST[a,b]   == DJac[c,a,b] dchi[c] + Jac[c,a] Jac[d,b] ddchi[c,d],
    dThetaSST[a]    == Jac[b,a] dTheta[b],
    
    (* now give the tmp derivative the same name again *)
    da[a]         == daSST[a],
    dda[a,b]      == ddaSST[a,b],
    db[a,b]       == dbSST[a,b],
    ddb[a,b,c]    == ddbSST[a,b,c],
    delg[a,b,c]   == delgSST[a,b,c],
    deldelg[a,b,c,d]  == deldelgSST[a,b,c,d],
    delG[a,b]     == delGSST[a,b],
    dKhat[a]      == dKhatSST[a],
    dA[a,b,c]     == dASST[a,b,c],
    dchi[a]       == dchiSST[a],
    ddchi[a,b]    == ddchiSST[a,b],
    dTheta[a]     == dThetaSST[a],
    
    (* now advection derivatives, FIXME: do it better or prove that this is enough at level 0 *)
    lieg[a,b]     == beta[c] delg[c,a,b],
    lieA[a,b]     == beta[c] dA[c,a,b],
    lieKhat       == beta[a] dKhat[a],
    liechi        == beta[a] dchi[a],
    liealpha      == beta[a] da[a],
    advG[a]       == beta[b] delG[b,a],
    advbeta[a]    == beta[b] db[b,a],
    advB[a]       == beta[b] dB[b,a],
    advTheta      == beta[a] dTheta[a],
    
  Cinstruction == "} else {",

    (* advection derivatives, treated differently than centered derivatives
       typically we store the advection derivative under lieSomething and
       add the missing terms later
    *)
    w1 == beta1/dx, w2 == beta2/dy, w3 == beta3/dz,
    Cinstruction == "if (order_advection == 2 || boundary1away) {",
      Cinstruction == "set_advection2(w1, w2, w3);",
      lieg[a,b] == adv[g[a,b]],
      lieA[a,b] == adv[A[a,b]],
      lieKhat == adv[Khat],
      liechi == adv[chi],
      liealpha == adv[alpha],
      advG[a] == adv[G[a]],
      advbeta[a] == adv[beta[a]],
      advB[a] == adv[B[a]],
      advTheta == adv[Theta],
    Cinstruction == "} else if (order_advection == 4 || boundaryNaway(2)) {",
      Cinstruction == "set_advection4(w1, w2, w3);",
      lieg[a,b] == adv4[g[a,b]],
      lieA[a,b] == adv4[A[a,b]],
      lieKhat == adv4[Khat],
      liechi == adv4[chi],
      liealpha == adv4[alpha],
      advG[a] == adv4[G[a]],
      advbeta[a] == adv4[beta[a]],
      advB[a] == adv4[B[a]],
      advTheta == adv4[Theta],
    Cinstruction == "} else if (order_advection == 6 || boundaryNaway(3)) {",
      Cinstruction == "set_advection6(w1, w2, w3);",
      lieg[a,b] == adv6[g[a,b]],
      lieA[a,b] == adv6[A[a,b]],
      lieKhat == adv6[Khat],
      liechi == adv6[chi],
      liealpha == adv6[alpha],
      advG[a] == adv6[G[a]],
      advbeta[a] == adv6[beta[a]],
      advB[a] == adv6[B[a]],
      advTheta == adv6[Theta],
    Cinstruction == "} else if (order_advection == 8 || boundaryNaway(4)) {",
      Cinstruction == "set_advection6(w1, w2, w3);",
      lieg[a,b] == adv8[g[a,b]],
      lieA[a,b] == adv8[A[a,b]],
      lieKhat == adv8[Khat],
      liechi == adv8[chi],
      liealpha == adv8[alpha],
      advG[a] == adv8[G[a]],
      advbeta[a] == adv8[beta[a]],
      advB[a] == adv8[B[a]],
      advTheta == adv8[Theta],
    Cinstruction == "} else errorexit(\"this order is not yet implemented\"); ",
    
  Cinstruction == "}",







  (* get K from Khat *)
  K == Khat + 2 Theta, 
  dK[a] == dKhat[a] + 2 dTheta[a],

  (* inverse conformal metric *)
  detginv == 1/matrixdet[g],
  ginv[a,b] == detginv matrixinvdet[g,a,b],

  (* derivatives of conformal metric *)
  dginv[a,b,c] == - ginv[b,d] ginv[c,e] delg[a,d,e],

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

  chiguard   == chiDivFloor,
  chiguarded == chi,
  Cinstruction == "if (chiguarded < chiguard) chiguarded = chiguard;",
  ff == chiguarded,
  oochipsipower == 1/chipsipower,
  f == oochipsipower log[ff],
  psim4 == exp[-4 f],
  (* phi == f, *)
  df[a] == oochipsipower dchi[a] / chiguarded,
  ddf[a,b] == oochipsipower ddchi[a,b] / chiguarded -chipsipower df[a] df[b],

  cddf[a,b] == ddf[a,b] - gamma[c,a,b] df[c],
  trcddf == ginv[a,b] cddf[a,b],


  (* curvature contribution from conformal factor *)
  Rphi[a,b] == 4 df[a] df[b] - 2 cddf[a,b] - 2 g[a,b] trcddf -
               4 g[a,b] ginv[c,d] df[c] df[d],

  (* derivatives of the lapse *)
  cdda[a,b] == dda[a,b] - gamma[c,a,b] da[c] -
               2 (df[a] da[b] + df[b] da[a] - g[a,b] ginv[c,d] df[c] da[d]),
  trcdda == psim4 (ginv[a,b] cdda[a,b]),

  (* tf part of conformal extrinsic curvature *)
  AA[a,b] == ginv[c,d] A[a,c] A[d,b],
  cAA == ginv[a,b] AA[a,b],
  Ainv[a,b] == ginv[a,c] ginv[b,d] A[c,d],
  divAinv[a] == gamma[a,b,c] Ainv[b,c] - 3/2 Ainv[a,b] dchi[b] / chiguarded
                - 1/3 ginv[a,b] (2 dKhat[b] + dTheta[b]),

  (* Ricci scalar *)
  (* In bams bssn:  R == cAA - 2 K K / 3,  <-- from Hamiltonian=0 
     here Rhat comes from the Ricci, but with constr. addition in R[a,b],
     since undifferentiated G[a] was replaced by its definition in R[a,b] *)
  Rhat == psim4 ginv[a,b] (R[a,b] + Rphi[a,b]),

  (* actual Hamiltonian constraint, but with constr. addition in R[a,b] *)
  Hhat == Rhat + 2 K K / 3 - cAA,

  (* the shift terms *)
  divbeta == delta[a,b] db[a,b],
  totdivbeta == 2/3 divbeta,
  ootddivbeta[a] == 1/3 delta[b,c] ddb[a,b,c],
  
  lieg[a,b] == lieg[a,b] - g[a,b] totdivbeta + db[a,c] g[b,c] + db[b,c] g[a,c],
  lieA[a,b] == lieA[a,b] - A[a,b] totdivbeta + db[a,c] A[b,c] + db[b,c] A[a,c],

  liechi == liechi + (chipsipower/6) chiguarded divbeta,

  pseudolieG[a] == ginv[b,c] ddb[b,c,a] + ginv[a,b] ootddivbeta[b] -
                   Gfromg[b] db[b,a] + Gfromg[a] totdivbeta + 
                   advG[a],






  (* right hand sides *)
  rg[a,b]   == - 2 alpha A[a,b] + lieg[a,b],

  rA[a,b]   == psim4 ( - cdda[a,b] + alpha (R[a,b] + Rphi[a,b])) -
               1/3 g[a,b] (- trcdda + alpha Rhat) +
               alpha (K A[a,b] - 2 AA[a,b]) + lieA[a,b],

  rG[a]     == - 2 Ainv[a,b] da[b] + 2 alpha ( divAinv[a] ) -
               2 alpha kappa1 ( G[a] - Gfromg[a] ) + pseudolieG[a],

  rKhat     == - trcdda + alpha (cAA + K K/3) + lieKhat + 
               kappa1 (1-kappa2) alpha Theta, 

  rchi      == (-chipsipower/6) chiguarded alpha K + liechi,

  rTheta    == ( alpha ( Hhat/2 - (2+kappa2) kappa1 Theta ) + advTheta ) * 
               zeroThetaRHS,

  Cinstruction == "CheckForNANandINF(18, rA11,rA12,rA13,rA22,rA23,rA33, rchi, rG1,rg11,rg12,rg13,rG2,rg22,rg23,rG3,rg33,rKhat,rTheta);",
  
  
  



  (* non principle terms *)
  Cinstruction == "if (useNPTs) { ",
    
    dginv[a, b,c] == - ginv[b,d] ginv[c,e] delg[a, d,e],
    
    (* Transformation for Christoffel from Alc S.85 *)
    (* chi = e^(-4phi), here phi is not the evolved var form Alic et al *)
    dphi[a]       == -1/(4*chiguarded) dchi[a],
    (* we have to use the full metric, but conformal factor cancels out *)
    gammaF[a,b,c] == gamma[a,b,c] + 2 (delta[a,b] dphi[c] + delta[a,c] dphi[b] - g[b,c] ginv[a,d] dphi[d]),
    
    
    Gd[b]       == ginv[b,c] ginv[d,e] delg[e, c,d],
    dGd[a,b]    == dginv[a, b,c] ginv[d,e] delg[e, c,d] +
                   ginv[b,c] dginv[a, d,e] delg[e, c,d] +
                   ginv[b,c] ginv[d,e] deldelg[a,e, c,d],
      
    Zinv[b]     == 1/2 ( G[b] - Gd[b] ),
    dZinv[a,b]  == 1/2 ( delG[a,b] - dGd[a,b] ),
    Z[b]        == g[b,c] Zinv[c],
    dZ[a,b]     == delg[a,b,c] Zinv[c] + g[b,c] dZinv[a,c],

    DZinv[a,b]  == dZinv[a,b] + gammaF[b,a,c] Zinv[c] ,
    DZ[a,b]     == dZ[a,b]    - gammaF[c,a,b] Z[c],
    DZsym[a,b]  == DZ[a,b] + DZ[b,a],
    trDZsym     == psim4 (ginv[a,b] DZsym[a,b]),
   
    rA[a,b] == rA[a,b] + chi alpha ( DZsym[a,b] - 1/3 g[a,b] trDZsym )
                       - 2 alpha A[a,b] Theta,
  
    rTheta  == rTheta  + alpha DZinv[a,a] - Zinv[a] da[a],
    
    rG[a]   == rG[a]   - 2 ginv[b,a] ( Theta da[b] + 2/3 alpha K Z[b] ) 
                       - 2 alpha kappa1 ginv[a,b] Z[b],

    Cinstruction == "if (CheckForNANandINF(10, rA11,rA12,rA13,rA22,rA23,rA33, rG1,rG2,rG3,rTheta)) errorexit(\"better stop\");",

  Cinstruction == "} else {",
    Zinv[a] == 0,
  Cinstruction == "}",
  
  
  
  


  (* gauge conditions, avoid unnecessary if statements, use flags *)
  falpha      == oploglapse lapseharmonicf + oploglapse2 8/3/(3-alpha) +
                 harmoniclapse alpha,
  ralpha      == oplogwithshift liealpha - falpha alpha Khat ,


  Cinstruction == "if (withB) { ",
  
    betaF     == shiftgammacoeff alpha^shiftalphapower,

    rbeta[a]  == (gamma2factor + gamma0factor betaF) B[a] + withShiftadv advbeta[a],

    rB[a]     == (gamma0factor + gamma2factor betaF) * (rG[a] - withGadv advG[a]) 
                 - shiftdriver B[a] + withBadv advB[a],

  Cinstruction == "} else { ",

    rB[a]     == 0,
    rbeta[a]  == G[a] + withShiftadv advbeta[a] - shiftdriver beta[a] ,

  Cinstruction == "} ",



  Cinstruction == "CheckForNANandINF(25, rA11,rA12,rA13,rA22,rA23,rA33, ralpha,rB1,rB2,rB3,rbeta1,rbeta2,rbeta3,rchi, rG1,rg11,rg12,rg13,rG2,rg22,rg23,rG3,rg33,rKhat,rTheta);",

  (* add dissipation terms
     we should test which terms are really needed
  *)
  Cinstruction == "if (order_dissipation == 4 && boundary2ormore) {",
    rg[a,b] == rg[a,b] + dissfactor dissipation4[g[a,b]],
    rA[a,b] == rA[a,b] + dissfactor dissipation4[A[a,b]],
    rG[a]   == rG[a]   + dissfactor dissipation4[G[a]],
    rKhat   == rKhat   + dissfactor dissipation4[Khat],
    rchi    == rchi    + dissfactor dissipation4[chi],
    rTheta  == rTheta  + dissfactor dissipation4[Theta],
    ralpha  == ralpha  + dissfactor dissipation4[alpha],
    rbeta[a]== rbeta[a]+ dissfactor dissipation4[beta[a]],
    rB[a]   == rB[a]   + dissfactor dissipation4[B[a]],
  Cinstruction == "}",
  Cinstruction == "if (order_dissipation == 6 && boundary3ormore) {",
    rg[a,b] == rg[a,b] + dissfactor dissipation6[g[a,b]],
    rA[a,b] == rA[a,b] + dissfactor dissipation6[A[a,b]],
    rG[a]   == rG[a]   + dissfactor dissipation6[G[a]],
    rKhat   == rKhat   + dissfactor dissipation6[Khat],
    rchi    == rchi    + dissfactor dissipation6[chi],
    rTheta  == rTheta  + dissfactor dissipation6[Theta],
    ralpha  == ralpha  + dissfactor dissipation6[alpha],
    rbeta[a]== rbeta[a]+ dissfactor dissipation6[beta[a]],
    rB[a]   == rB[a]   + dissfactor dissipation6[B[a]],
  Cinstruction == "}",
  Cinstruction == "if (order_dissipation == 8 && boundary4ormore) {",
    rg[a,b] == rg[a,b] + dissfactor dissipation8[g[a,b]],
    rA[a,b] == rA[a,b] + dissfactor dissipation8[A[a,b]],
    rG[a]   == rG[a]   + dissfactor dissipation8[G[a]],
    rKhat   == rKhat   + dissfactor dissipation8[Khat],
    rchi    == rchi    + dissfactor dissipation8[chi],
    rTheta  == rTheta  + dissfactor dissipation8[Theta],
    ralpha  == ralpha  + dissfactor dissipation8[alpha],
    rbeta[a]== rbeta[a]+ dissfactor dissipation8[beta[a]],
    rB[a]   == rB[a]   + dissfactor dissipation8[B[a]],
  Cinstruction == "}",
  

  Cinstruction == "if (CheckForNANandINF(25, rA11,rA12,rA13,rA22,rA23,rA33, ralpha,rB1,rB2,rB3,rbeta1,rbeta2,rbeta3,rchi, rG1,rg11,rg12,rg13,rG2,rg22,rg23,rG3,rg33,rKhat,rTheta)) {
    printf(\"x=%2.5e, y=%2.5e, z=%2.5e\\n\",xp[ijk],yp[ijk],zp[ijk]);
    printf(\"  (%e %e  %e %e)\\n\",alpha[ijk],palpha[ijk],Khat[ijk],pKhat[ijk]);
    printf(\"  %e %e %e %e\\n\",detginv,oochipsipower,chiguarded,psim4);
    printf(\"  %e %e %e   %e %e %e %e %e %e\\n\",A33[ijk],lieg33,rg33, rA11,rA12,rA13,rA22,rA23,rA33);
    errorexit(\"better stop\");}",


  (* result *) 
  Cif == addlinear, 
    ng[a,b]     == pg[a,b]  + c rg[a,b],
    nA[a,b]     == pA[a,b]  + c rA[a,b],
    nG[a]       == pG[a]    + c rG[a],
    nKhat       == (pKhat   + c rKhat),
    nchi        == pchi     + c rchi,
    nTheta      == pTheta   + c rTheta,
    nalpha      == palpha    + c ralpha,
    nbeta[a]    == pbeta[a] + c rbeta[a],
    nB[a]       == pB[a]    + c rB[a], 

    (* normalize detg and subtract trace from As *)
    detnginv    == 1/matrixdet[ng],
    Cif == subtractA,
      traceA    == detnginv matrixinvdet[ng,a,b] nA[a,b],
      aux       == -traceA/3,
      nA[a,b]   == nA[a,b] + aux ng[a,b],
    Cif == end,
    Cif == normalizedetg,
      aux       == detnginv^(1/3),
      ng[a,b]   == aux ng[a,b],
    Cif == end,

  Cif == else,
    ng[a,b]     == rg[a,b],
    nA[a,b]     == rA[a,b],
    nG[a]       == rG[a],
    nKhat       == rKhat,
    nchi        == rchi,
    nTheta      == rTheta,
    nalpha      == ralpha,
    nbeta[a]    == rbeta[a],
    nB[a]       == rB[a],
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
Khat[a_,b_] := Khat[b,a] /; !OrderedQ[{a,b}]
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



ddchiSST[a_,b_]   := ddchiSST[b,a] /; !OrderedQ[{a,b}]
ddaSST[a_,b_]     := ddaSST[b,a] /; !OrderedQ[{a,b}]
ddbSST[a_,b_,c_]  := ddbSST[b,a,c] /; !OrderedQ[{a,b}]
delgSST[c_,a_,b_] := delgSST[c,b,a] /; !OrderedQ[{a,b}]
dASST[c_,a_,b_] := dASST[c,b,a] /; !OrderedQ[{a,b}]
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
  pr["#include \"z4.h\"\n"];
  pr["\n"];
  pr["#define Power(x,y) (pow((double) (x), (double) (y)))\n"];
  pr["#define Sqrt(x)    sqrt(x)\n"];
  pr["#define Log(x)     log((double) (x))\n"];
  pr["#define pow2(x)    ((x)*(x))\n"];
  pr["#define pow3(x)    ((x)*(x)*(x))\n"];
  pr["#define pow4(x)    ((x)*(x)*(x)*(x))\n"];
  pr["#define pow2inv(x) (1.0/((x)*(x)))\n\n"];
  pr["#define Tanh(x)    tanh(x)\n"];
  pr["#define Sech(x)    (1/cosh(x))\n"];
  pr["#define Cal(x,y,z) ((x)?(y):(z))\n\n"];

  pr["\n\n\n"];
];

functionname      = "z4_rhs_movpunc_N";

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

  prdecvlNR[{ g[a,b],  chi,  A[a,b],  Khat,  G[a],  Theta,  alpha,  beta[a], B[a]}, "ucur", "METRIC_z4_INDX_VAR"];
  prdecvlNR[{ng[a,b], nchi, nA[a,b], nKhat, nG[a], nTheta, nalpha, nbeta[a],nB[a]}, "unew", "METRIC_z4_INDX_VAR"];
  prdecvlNR[{pg[a,b], pchi, pA[a,b], pKhat, pG[a], pTheta, palpha, pbeta[a],pB[a]}, "upre", "METRIC_z4_INDX_VAR"];
  prdecvarname[{Zinv[a]}, "z4_Zx"];
  pr["\n"];
 

  pr["const int addlinear = (c != 0.0l);\n"];
  
  pr["const int order_centered      = Geti(\"order_centered\");\n"];
  pr["const int order_advection     = Geti(\"order_advection\");\n"];
  pr["const int advectionlopsided   = Geti(\"advection_lopsided\");\n"];
  pr["const int advectionlopsided6  = Geti(\"advection_lopsided6\");\n"];
  pr["const int order_dissipation   = Geti(\"order_dissipation\");\n"];
  pr["const double dissfactor       = get_dissipation_factor(level);\n"];

  pr["const double kappa1       = Getd(\"z4_kappa1\");\n"];
  pr["const double kappa2       = Getd(\"z4_kappa2\");\n"];
  pr["const double zeroThetaRHS = Getv(\"z4_zeroThetaRHS\", \"no\");\n"];
  pr["const double useNPTs      = Getv(\"z4_useNPTs\",\"yes\");\n"];

  pr["const double forceKzerofactor = Getv(\"z4_forceKzero\", \"no\");\n"];
  pr["const int subtractA      = Getv(\"z4_subtractA\", \"yes\");\n"];
  pr["const int normalizedetg  = Getv(\"z4_normalizedetg\", \"yes\");\n"];
  pr["const int constantlapse  = Getv(\"z4_lapse\", \"constant\");\n"];
  pr["const int oploglapse     = Getv(\"z4_lapse\", \"1+log\");\n"];
  pr["const int oploglapse2    = Getv(\"z4_lapse\", \"1+log2\");\n"];
  pr["const int oplogwithshift = Getv(\"z4_lapse\", \"withshift\");\n"];
  pr["const int harmoniclapse  = Getv(\"z4_lapse\", \"harmonic\");\n"];

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

  pr["const double chiDivFloor = Getd(\"z4_chi_div_floor\");\n"];
  pr["const double chipsipower = Getd(\"z4_chi_psipower\");\n"];

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

  (* extrapolate fields for background BCs on boxes  *)
  pr["if (!Getv(\"grid\",\"shells\") && Getv(\"boundary\", \"background\") && (level->l==0) ) {\n/* extrapolate field values */\n"];
  pr["int vindex;\n"];
  pr["for (vindex = 0; vindex < ucur->n; vindex++) {\n"];
  pr["	set_boundary_extrapolate(level, ucur->index[vindex]);}}\n"];

];    




(************************************************************************)
(* now we are ready to go *)

(* assume that we are working on a box *)
BoxMode = True;


<< "../../math/MathToC/TensorEquationsToC.m"
