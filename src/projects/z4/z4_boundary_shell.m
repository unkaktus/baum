
(************************************************************************)
(* Second order constraint preserving boundary conditions *)
(************************************************************************)

(* background.m mth 02/10 *)
(* background.m dh  08/10 *)
(* background.m mth 05/11   -> now only shells *)
(* Z4d_boundary.m mth 06/11 -> moved from boundary to Z4d *)
(* Z4d_boundary.m dh 10/11  -> Conditions without linearization *)
(* z4_boundary_shell.m sb 01/15  -> just changed name *)

(* variables *)
variables = { g[a,b],  chi,  A[a,b],  Khat,  G[a],  Theta,  alpha,  beta[a],  B[a], 
             ng[a,b], nchi, nA[a,b], nKhat, nG[a], nTheta, nalpha, nbeta[a], nB[a],
             pg[a,b], pchi, pA[a,b], pKhat, pG[a], pTheta, palpha, pbeta[a], pB[a],
             rp,rr,xp,yp,zp
            }

(* compute in this order *)
tocompute = {

  Cinstruction == "if (CheckForNANandINF(25, 
    g11[ijk],g12[ijk],g13[ijk],g22[ijk],g23[ijk],g33[ijk],
    A11[ijk],A12[ijk],A13[ijk],A22[ijk],A23[ijk],A33[ijk], 
    G1[ijk],G2[ijk],G3[ijk], Khat[ijk],chi[ijk], Theta[ijk],
    alpha[ijk],beta1[ijk],beta2[ijk],beta3[ijk],B1[ijk],B2[ijk],B3[ijk])) {
    printf(\"problem with vars in Z4d_boundary.m\\n\");
    printf(\"x=%2.5e, y=%2.5e, z=%2.5e\\n\",xp[ijk],yp[ijk],zp[ijk]);}", 



(************************************************************************)
(* ---- Compute derivatives on the spherical grid ---- *)
(************************************************************************)

(* ---- derivatives of cartesian components in spherical directions (del_rpt V^xyz) ---- *)
    Cinstruction     == "if (order == 2 || boundaryNaway(1)) {",
        da[a]       == del[a,alpha],
        dda[a,b]    == deldel[a,b,alpha],
        db[a,b]     == del[a,beta[b]],
        ddb[a,b,c]  == deldel[a,b,beta[c]],
        dchi[a]     == del[a,chi],
        ddchi[a,b]  == deldel[a,b,chi],
        dg[c,a,b]   == del[c,g[a,b]],
        ddg[a,b,c,d]== deldel[a,b,g[c,d]],
        dKhat[a]    == del[a,Khat],
        dA[a,b,c]   == del[a,A[b,c]],
        dG[a,b]     == del[a,G[b]],
        dTheta[a]   == del[a,Theta],
    Cinstruction == "} else if (order == 4 || boundaryNaway(2)) {",
        da[a]       == del4[a,alpha],
        dda[a,b]    == deldel4[a,b,alpha],
        db[a,b]     == del4[a,beta[b]],
        ddb[a,b,c]  == deldel4[a,b,beta[c]],
        dchi[a]     == del4[a,chi],
        ddchi[a,b]  == deldel4[a,b,chi],
        dg[c,a,b]   == del4[c,g[a,b]],
        ddg[a,b,c,d]== deldel4[a,b,g[c,d]],
        dKhat[a]    == del4[a,Khat],
        dA[a,b,c]   == del4[a,A[b,c]],
        dG[a,b]     == del4[a,G[b]],
        dTheta[a]   == del4[a,Theta],
    Cinstruction == "} else if (order == 6 || boundaryNaway(3)) {",
        da[a]       == del6[a,alpha],
        dda[a,b]    == deldel6[a,b,alpha],
        db[a,b]     == del6[a,beta[b]],
        ddb[a,b,c]  == deldel6[a,b,beta[c]],
        dchi[a]     == del6[a,chi],
        ddchi[a,b]  == deldel6[a,b,chi],
        dg[c,a,b]   == del6[c,g[a,b]],
        ddg[a,b,c,d]== deldel6[a,b,g[c,d]],
        dKhat[a]    == del6[a,Khat],
        dA[a,b,c]   == del6[a,A[b,c]],
        dG[a,b]     == del6[a,G[b]],
        dTheta[a]   == del6[a,Theta],
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
    daSST[a]        == Jac[b,a] da[b],
    ddaSST[a,b]     == DJac[c,a,b] da[c] + Jac[c,a] Jac[d,b] dda[c,d],
    dbSST[a,b]      == Jac[c,a] db[c,b],
    ddbSST[a,b,c]   == DJac[d,a,b] db[d,c] + Jac[d,a] Jac[e,b] ddb[d,e,c],
    dgSST[a,b,c]  == Jac[d,a] dg[d,b,c],
    ddgSST[a,b,c,d] == DJac[e,a,b] dg[e,c,d] + Jac[e,a] Jac[f,b] ddg[e,f,c,d],
    dGSST[a,b]    == Jac[c,a] dG[c,b],
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
    dg[a,b,c]   == dgSST[a,b,c],
    ddg[a,b,c,d]  == ddgSST[a,b,c,d],
    dG[a,b]     == dGSST[a,b],
    dKhat[a]      == dKhatSST[a],
    dA[a,b,c]     == dASST[a,b,c],
    dchi[a]       == dchiSST[a],
    ddchi[a,b]    == ddchiSST[a,b],
    dTheta[a]     == dThetaSST[a],

(************************************************************************)
(* ---- Define projection stuff in xyz and rpt; up and down differ! ---- *)
(************************************************************************)

    chiguarded == chi,
    chiguard   == chiDivFloor,
    Cinstruction == "if (chiguarded<chiguard) {
      printf(\"chi is wrong (%e)\\n\",chi[ijk]);
      chiguarded = chiguard;
    }",


    (* ---- Define normal vector ---- *)
    r    == 0,
    shat[a] == 0,
    sup[a] == 0,    
    Cinstruction == "r=rp[ijk];",
    Cinstruction == "shat1=xp[ijk]/r;shat2=yp[ijk]/r;shat3=zp[ijk]/r;",

    detginv == 1/matrixdet[g],
    Cinstruction == "if (detginv<0.00001) {
      printf(\"detginv is wrong (%e)\\n\",detginv);
      detginv = 1.0;
    }",

    ginv[a,b]    == detginv matrixinvdet[g,a,b],
    ADMginv[a,b] == chiguarded ginv[a,b],
    modshatARG   == ADMginv[c,d]shat[c]shat[d],
    Cinstruction == "if (modshatARG<0.00001) {
      printf(\"modshat is wrong (%e)\\n\",modshatARG);
      // modshatARG = 0.00001;
      modshatARG = pow2(shat1) + pow2(shat2) + pow2(shat3);
    }",
    oomodshat    == 1/sqrt[modshatARG],
    sdown[a]     == oomodshat shat[a],
    sup[a]       == ADMginv[a,b] sdown[b],
        
    (* ---- Define 2+1 projection operator ---- *)
    qud[a,b] == delta[a,b] - sup[a]sdown[b],
    qdd[a,b] == g[a,b] / chiguarded - sdown[a]sdown[b],
    quu[a,b] == ADMginv[a,b] - sup[a]sup[b],

    (* ---- Physical projection ---- *)
    qPhysuudd[c,d,a,b] == (qud[c,a]qud[d,b]+qud[c,b]qud[d,a]
                          - quu[c,d]qdd[a,b])/2, 

(* ---- Speeds for the boundary conditions ---- *)
    (* ---- lapse: muL ---- *)
    muL    ==  2 / alpha,
    (* ---- Shift: muS ---- *)
    muStilde ==  1/chiguarded,
    (* ---- Longtitudinal shift speed ---- *)
    vbetas == 2 sqrt[ muStilde / 3 ],
    (* ---- Transverse shift speed ---- *)
    vbetaA == sqrt[muStilde],

(************************************************************************)
(* ---- Connections, curvatures, divs, aux vars ---- *)
(************************************************************************)

   (* get K from Khat *)
   K == Khat + 2 Theta, 
   dK[a] == dKhat[a] + 2 dTheta[a],

   (* derivatives of conformal metric *)
   dginv[a,b,c] == - ginv[b,d] ginv[c,e] dg[a,d,e],

   (* connection of conformal metric *)
   gammado[c,a,b] == 1/2 (dg[a,b,c] + dg[b,a,c] - dg[c,a,b]),
   gamma[c,a,b] == ginv[c,d] gammado[d,a,b],  
   Gfromg[a] == ginv[b,c] gamma[a,b,c],
   dGfromgdu[a,b] == ginv[b,c]ginv[d,e]ddg[a,d,e,c]
                     -ginv[b,c]ginv[d,f]ginv[e,g]dg[g,e,f]dg[a,c,d]
                     -ginv[b,e]ginv[g,c]ginv[d,f]dg[g,e,f]dg[a,c,d],

   (* curvature of conformal metric *)
   R[a,b] == ginv[c,d] ( -1/2 ddg[c,d,a,b] +
                     gamma[e,c,a] gammado[b,e,d] +
                     gamma[e,c,b] gammado[a,e,d] +
                     gamma[e,a,d] gammado[e,c,b]) +
                     1/2 (g[c,a] dG[b,c] + g[c,b] dG[a,c] + 
                     Gfromg[c] (gammado[a,b,c] + gammado[b,a,c])),

   ff == chiguarded,
   oochipsipower == 1/chipsipower,
   f == oochipsipower log[ff],
   psim4 == exp[-4 f],
   (* phi == f, *)
   df[a] == oochipsipower dchi[a] / chiguarded,
   ddf[a,b] == oochipsipower ddchi[a,b] / chiguarded - chipsipower df[a] df[b],

   cddf[a,b] == ddf[a,b] - gamma[c,a,b] df[c],
   trcddf == ginv[a,b] cddf[a,b],

   (* curvature contribution from conformal factor *)
   Rphi[a,b] == 4 df[a] df[b] - 2 cddf[a,b] - 2 g[a,b] trcddf -
                4 g[a,b] ginv[c,d] df[c] df[d],

   (* full curvature *)
   Rf[a,b] == R[a,b] + Rphi[a,b],
   Rhat == psim4 ginv[a,b] Rf[a,b],

   (* derivatives of the lapse *)
   cdda[a,b] == dda[a,b] - gamma[c,a,b] da[c] -
                2 (df[a] da[b] + df[b] da[a] - g[a,b] ginv[c,d] df[c] da[d]),
   trcdda == psim4 (ginv[a,b] cdda[a,b]),

   (* tf part of conformal extrinsic curvature *)
   AA[a,b] == ginv[c,d] A[a,c] A[d,b],
   Ainv[a,b] == ginv[a,c] ginv[b,d] A[c,d],

   (* cov div. A *)
   cdA[c,a,b] == dA[c,a,b] - gamma[d,c,a]A[b,d] - gamma[d,c,b]A[a,d], 

   (* the shift terms *)
   divbeta == delta[a,b]db[a,b],
   totdivbeta == 2/3 divbeta,
  
   (* Lie derivatives *)
   lieg[a,b] == beta[c]dg[c,a,b] - g[a,b] totdivbeta + db[a,c] g[b,c] 
                + db[b,c] g[a,c],
   lieA[a,b] == beta[c]dA[c,a,b] - A[a,b] totdivbeta + db[a,c] A[b,c] 
                + db[b,c] A[a,c],

(************************************************************************)
(* ---- Define some 2+1 split variables ---- *)
(************************************************************************)

   (* ---- Scalar sector  ---- *)
   betas  == sdown[a]beta[a],
   Dbetas == sup[a]sdown[b]db[a,b],
   Dalpha == sup[a]da[a],
   DKhat  == sup[a]dKhat[a],
   DK     == sup[a]dK[a],
   DTheta == sup[a]dTheta[a],
   Gams   == sdown[a]G[a],
   DGams  == sup[a]sdown[b]dG[a,b],


   (* ---- Vector sector  ---- *)
   GamA[k]    == qud[k,a]G[a],
   DGamA[k]   == sup[a]qud[k,b]dG[a,b],
   betaA[k]   == qud[k,a]beta[a],
   DbetaA[k]  == sup[a]qud[k,b]db[a,b],

   (* ---- normal der. for conditions required in others  ---- *)

  lienKhat == - sqrt[muL] DKhat - sqrt[muL] Khat / r, 
  
  (*
   lienKhat == - sqrt[muL] DKhat - sqrt[muL] Khat / r
   + Dalpha / r / alpha - betas Khat / r / alpha 
   + betas Dalpha / r / sqrt[muL] / alpha / alpha
   - quu[a,b]cdda[a,b] / 2 / alpha
   - 1/2 (ginv[i,j]AA[i,j] + K K / 3 
	  + kappa1 (1-kappa2) Theta), 
   *)

   lienTheta == - (DTheta + Theta / r) - kappa1 (2 + kappa2) Theta,

   lienK == lienKhat + 2 lienTheta, 

(************************************************************************)
(* ---- Construct RHSs of the 2+1 decomposed quantities ---- *)
(************************************************************************)

(* ---- Scalar sector ---- *)
    (* ---- Gauge ---- *)
    rKhat     == beta[a] dKhat[a] + alpha lienKhat,

    rGams     == - vbetas dG[a,a] + sdown[b]beta[a]dG[a,b]
                + 2 shiftdriver/sqrt[3]/vbetaA db[a,a]
                + quu[a,b]sdown[c]ddb[a,b,c]/chiguarded
                - sup[a]qud[b,c]ddb[a,b,c]/chiguarded
                + 4alpha sqrt[muL]/(3(vbetas+sqrt[muL]))sup[a]dKhat[a]/chiguarded
                + 2alpha/(3(vbetas+1))sup[a]dTheta[a]/chiguarded,
  (*  
  rGams     == shiftdriver beta[i]sdown[j]dG[i,j]
               + shiftdriver beta[i]beta[j]sdown[k]ddb[i,j,k]
               + shiftdriver sdown[k]beta[i]db[i,j]db[j,k]
               - shiftdriver shiftdriver sdown[j]beta[i]db[i,j],
   *)

    (* ---- Constraints ---- *)
    rTheta    == alpha lienTheta + beta[a]dTheta[a],
    rACss     == - alpha chiguarded (2sup[b]ginv[c,a]cdA[c,a,b]
                            - 4/3sup[a]dK[a]+2DTheta 
                            - 3sup[a]ginv[b,c]dchi[c]A[a,b]/chiguarded
                            - kappa1 chiguarded sdown[a](G[a]-Gfromg[a])
                              + 2/3chiguarded sup[a]sdown[b](dG[a,b]-dGfromgdu[a,b])
                              - 1/3chiguarded qud[a,b](dG[a,b]-dGfromgdu[a,b])
                            )
                 + sup[a]sup[b](psim4(-cdda[a,b]+alpha(Rf[a,b]))-
                 1/3 g[a,b](-trcdda+alpha Rhat)+
                 alpha (K A[a,b]-2AA[a,b])+lieA[a,b] ),
    rACqq     == chiguarded Ainv[a,b](-2alpha A[a,b]+lieg[a,b])-rACss,

(* ---- Vector sector ---- *)
    (* ---- Gauge ---- *)
    (* General condition under construction. *)
    rGamA[k]     == -vbetaA(sup[a]qud[k,b]dG[a,b]-quu[k,a]sdown[b]dG[a,b]/chiguarded)
                    +quu[a,b]qud[k,c]ddb[a,b,c]/chiguarded
                    +4sup[b]sdown[c]quu[k,a]ddb[a,b,c]/3/chiguarded
                    +qud[b,c]quu[k,a]ddb[a,b,c]/3/chiguarded
                    -4alpha quu[k,a]dKhat[a]/3/chiguarded-2alpha quu[k,a]dTheta[a]/3/chiguarded
                    +shiftdriver/(chiguarded vbetaA)(sup[a]qud[k,b]db[a,b]
                    -sdown[a]quu[k,b]db[b,a])
                    +qud[k,b]beta[a]dG[a,b], 

(* ---- Constraints ---- *)
    rACsA[a]     == - alpha chiguarded (ginv[b,c]qud[d,a]cdA[b,c,d]
                          -2/3qud[d,a]dK[d]+qud[d,a]dTheta[d]  
                          - 3qud[d,a]ginv[b,c]A[b,d]dchi[c]/(2chi)
                          - kappa1 chiguarded qdd[d,a](G[d]-Gfromg[d])/2 
                          - qud[d,a]sup[b]Rf[b,d]
                          + chiguarded sup[b]qdd[d,a](dG[b,d]-dGfromgdu[b,d])/2 
                          )
                    - chiguarded qud[d,a]sup[b]cdda[b,d]
                    + alpha sup[b]qud[d,a] (A[b,d] K - 2 AA[b,d])
                    + sup[b]qud[d,a]lieA[b,d],

(* ---- Tensor sector ---- *)
    rACABTF[a,b] == - alpha qPhysuudd[c,d,a,b](sup[e]cdA[e,c,d]
                                          -1/2 sup[e]cdA[c,d,e]
                                          -1/2 sup[e]cdA[d,c,e]
                                          +1/(4chiguarded)sup[e]A[e,c]dchi[d]
                                          +1/(4chiguarded)sup[e]A[e,d]dchi[c]
                                          -1/(2chiguarded)A[c,d]sup[e]dchi[e]
                                          )
                    -chiguarded qPhysuudd[c,d,a,b]cdda[c,d] 
                    -alpha qPhysuudd[c,d,a,b]AA[c,d]
                    +2/3alpha K qPhysuudd[c,d,a,b]A[c,d]
                    +qPhysuudd[c,d,a,b]lieA[c,d],

(* ---- add incomming radiation? ---- *)
     Cinstruction == "if (givehPsi0) {" ,

        (* specify incoming radiation by adding terms hPsi0 *)

        (* construct angular normal vectors *)             
        (* Gram-Schmidt taken from Invariants/curvature_invariants_N.m *)
        (* This way we control exactly the thing computed in post-processing. *)

        gADM[a,b] == g[a,b] / chiguarded,

        vu[a] == -yp delta[a,1] + xp delta[a,2],
        wu[a] == chiguarded^(-3/2) ADMginv[a,d] epsmatrix[d,b,c] sup[b] vu[c],

        sdotv ==  gADM[a,b] sup[a] vu[b],
        vu[a] == vu[a] - sup[a] sdotv,
        vdotv == gADM[a,b] vu[a] vu[b],
        vu[a] == vu[a]/(vdotv)^(1/2),
        sdotw == gADM[a,b] sup[a] wu[b] ,
        vdotw == gADM[a,b] vu[a] wu[b],
        wu[a] == wu[a] - sup[a] sdotw - vu[a] vdotw,
        wdotw == gADM[a,b] wu[a] wu[b],
        wu[a] == wu[a]/(wdotw)^(1/2),

        vd[a] == gADM[a,b]vu[b],
        wd[a] == gADM[a,b]wu[b],

        (* free data *)
        (* for the moment we can specify only Gaussian data *)
        RehPsi0 == hPsi0para Exp[ - (hPsi0parb (time-hPsi0parc ))^2 ],
        ImhPsi0 == 0,

        (* add Psi0 terms to curvature *) 
        rACABTF[a,b] == rACABTF[a,b] 
                           + alpha chiguarded RehPsi0 (vd[a]vd[b]-wd[a]wd[b]) 
                           + alpha chiguarded ImhPsi0 (vd[a]wd[b]+vd[b]wd[a]),

     Cinstruction == " } ",



(* ---- Reconstruct full RHSs from decomposed ones ---- *)

    (* ---- Extrinsic curvature ---- *)
    rA[a,b]   == rACABTF[a,b]+(qud[k,a]sdown[b]+qud[k,b]sdown[a])rACsA[k] 
                 + sdown[a]sdown[b]rACss+qdd[a,b]/2 rACqq,
    (* ---- G and Theta ---- *)
    rG[a]     == sup[a]rGams+qud[a,k]rGamA[k],

(* check for nans *)
  (* all devisions: 
    df,ddf,rACABTF is equal to 1/chiguarded, tested in qdd
    1/r should be ok
  *)
  Cinstruction == "if (CheckForNANandINF(15, 
    detginv, sdown1,sdown2,sdown3, sup1,sup2,sup3,
    qdd11,qdd12,qdd13,qdd22,qdd23,qdd33,
    muL,lienKhat)) {
    printf(\"problem with devisions in Z4_boundary.m\\n\");
    printf(\"x=%2.5e, y=%2.5e, z=%2.5e, r=%2.5e\\n\",xp[ijk],yp[ijk],zp[ijk],r);
  }",

  (* check RHS *)
  Cinstruction == "if (CheckForNANandINF(11, 
    rA11,rA12,rA13,rA22,rA23,rA33, 
    rG1,rG2,rG3,rKhat,rTheta)) {
    printf(\"nans in RHS in Z4_boundary.m\\n\");
    printf(\"x=%2.5e, y=%2.5e, z=%2.5e\\n\",xp[ijk],yp[ijk],zp[ijk]);
  }",

  (* Cinstruction == "printf(\" %03d %03d   %+2.8e   %+2.8e   %+2.8\\n\",j,k,rTheta, nTheta[ijk-di], nTheta[ijk+di]);", *)

(************************************************************************)
(* ---- Dissipation ---- *)
(************************************************************************)

  (* add dissipation terms we should test which terms are really needed *)
  Cinstruction == "if (order_dissipation == 4 && boundary2ormore) {",
    rA[a,b] == rA[a,b] + dissfactor dissipation4[A[a,b]],
    rG[a]   == rG[a]   + dissfactor dissipation4[G[a]],
    rKhat   == rKhat   + dissfactor dissipation4[Khat],
    rTheta  == rTheta  + dissfactor dissipation4[Theta],
  Cinstruction == "}",
  Cinstruction == "if (order_dissipation == 6 && boundary3ormore) {",
    rA[a,b] == rA[a,b] + dissfactor dissipation6[A[a,b]],
    rG[a]   == rG[a]   + dissfactor dissipation6[G[a]],
    rKhat   == rKhat   + dissfactor dissipation6[Khat],
    rTheta  == rTheta  + dissfactor dissipation6[Theta],
  Cinstruction == "}",
  Cinstruction == "if (order_dissipation == 8 && boundary4ormore) {",
    rA[a,b] == rA[a,b] + dissfactor dissipation8[A[a,b]],
    rG[a]   == rG[a]   + dissfactor dissipation8[G[a]],
    rKhat   == rKhat   + dissfactor dissipation8[Khat],
    rTheta  == rTheta  + dissfactor dissipation8[Theta],
  Cinstruction == "}",



(************************************************************************)
(* --- The rest of the script is standard BAM setup ... update the RHS---- *)
(************************************************************************)

   Cif == addlinear, 
       nKhat       == pKhat   + c rKhat,
       nA[i,j]     == pA[i,j] + c rA[i,j],
       nG[i]       == pG[i]   + c rG[i],
       nTheta      == pTheta  + c rTheta,
   Cif == else,
       nKhat       == rKhat,
       nA[i,j]     == rA[i,j],
       nG[i]       == rG[i],
       nTheta      == rTheta,
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
  lieA[a_,b_] := lieA[b,a] /; !OrderedQ[{a,b}]
  
  dda[a_,b_]      := dda[b,a] /; !OrderedQ[{a,b}]
  cdda[a_,b_]     := dda[b,a] /; !OrderedQ[{a,b}]
  ddb[a_,b_,c_]   := ddb[b,a,c] /; !OrderedQ[{a,b}]
  ddchi[a_,b_]    := ddchi[b,a] /; !OrderedQ[{a,b}]
  ddf[a_,b_]      := ddf[b,a] /; !OrderedQ[{a,b}]
  cddf[a_,b_]     := cddf[b,a] /; !OrderedQ[{a,b}]
  dg[c_,a_,b_]    := dg[c,b,a] /; !OrderedQ[{a,b}]
  ddg[a_,b_,c_,d_]:= ddg[b,a,c,d] /; !OrderedQ[{a,b}]
  ddg[a_,b_,c_,d_]:= ddg[a,b,d,c] /; !OrderedQ[{c,d}]
  dA[c_,a_,b_]    := dA[c,b,a] /; !OrderedQ[{a,b}]
  cdA[c_,a_,b_]   := cdA[c,b,a] /; !OrderedQ[{a,b}]
  dginv[c_,a_,b_] := dginv[c,b,a] /; !OrderedQ[{a,b}]

  R[a_,b_]      := R[b,a] /; !OrderedQ[{a,b}]
  Rphi[a_,b_]   := Rphi[b,a] /; !OrderedQ[{a,b}]
  Rf[a_,b_]     := Rf[b,a] /; !OrderedQ[{a,b}]

  DJac[a_,b_,c_]  := DJac[a,c,b] /; !OrderedQ[{b,c}]

  ddchiSST[a_,b_]         := ddchiSST[b,a] /; !OrderedQ[{a,b}] 
  ddaSST[a_,b_]           := ddaSST[b,a] /; !OrderedQ[{a,b}]
  ddbSST[a_,b_,c_]        := ddbSST[b,a,c] /; !OrderedQ[{a,b}]
  dgSST[a_,b_,c_]         := dgSST[a,c,b] /; !OrderedQ[{b,c}]
  dASST[a_,b_,c_]         := dASST[a,c,b] /; !OrderedQ[{b,c}]        
  ddgSST[a_,b_,c_,d_]     := ddgSST[b,a,c,d]  /; !OrderedQ[{a,b}]
  ddgSST[a_,b_,c_,d_]     := ddgSST[a,b,d,c]  /; !OrderedQ[{c,d}]

  qdd[a_,b_]         := qdd[b,a] /; !OrderedQ[{a,b}]
  quu[a_,b_]         := quu[b,a] /; !OrderedQ[{a,b}]
  qPhysuudd[a_,b_,c_,d_] := qPhysuudd[b,a,c,d]  /; !OrderedQ[{a,b}]
  qPhysuudd[a_,b_,c_,d_] := qPhysuudd[a,b,d,c]  /; !OrderedQ[{c,d}]

  rACABTF[a_,b_]  := rACABTF[b,a] /; !OrderedQ[{a,b}]

  ADMginv[a_,b_]  := ADMginv[b,a] /; !OrderedQ[{a,b}]
  ginv[a_,b_]     := ginv[b,a] /; !OrderedQ[{a,b}]

  rg[a_,b_]       := rg[b,a] /; !OrderedQ[{a,b}]
  rA[a_,b_]       := rA[b,a] /; !OrderedQ[{a,b}]
 

  (* Taken from Invariants/curvature_invariants_N.m *)
  (* The totally antisymmetric matrix *)
  epsmatrix[a_Integer,b_Integer,c_Integer] := (b-a)(c-b)(c-a)/2
  epsmatrix[a_,b_,c_] := -epsmatrix[b,a,c] /; !OrderedQ[{a,b}]
  epsmatrix[a_,b_,c_] := -epsmatrix[a,c,b] /; !OrderedQ[{b,c}] 
  
(************************************************************************)
(* information for C output *)
(************************************************************************)

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

functionname      = "z4_boundary_shell";

functionarguments = "tVarList *unew, tVarList *upre, double c, tVarList *ucur";

functionbegin     = "bampi_openmp_start\n"
;
functionloop      = 
  "forallpoints_ijk_openmp(level) {\n"<>
  "if (!dequal(rp[ijk],level->bbox[1]-N*level->dx)) continue;\n" <>
  "if (j==0 || j==box->n-1 || k==0 || k==box->o-1) continue;\n";

endfunctionloop   = "} endfor_ijk_openmp; /* loop i, j, k */"

endfunction       = "bampi_openmp_stop\n";




(* custom variable declaration
   we have to translate between the tensor names and the C variables  
   we may need special variables, for example for parameters
*)
variabledeclarations[] := Module[{},

  pr["tL *level = ucur->level;\n"];
  pr["int addlinear = (c != 0.0l);\n"];
  pr["double time = level->time;\n"];                    
  pr["\n"];
  
  prdecvlNR[{ g[a,b],  chi,  A[a,b],  Khat,  G[a],  Theta,  alpha,  beta[a],  B[a]}, "ucur", "METRIC_z4_INDX_VAR"];
  prdecvlNR[{ng[a,b], nchi, nA[a,b], nKhat, nG[a], nTheta, nalpha, nbeta[a], nB[a]}, "unew", "METRIC_z4_INDX_VAR"];
  prdecvlNR[{pg[a,b], pchi, pA[a,b], pKhat, pG[a], pTheta, palpha, pbeta[a], pB[a]}, "upre", "METRIC_z4_INDX_VAR"];
  pr["\n"];
  
  pr["const int order               = Geti(\"order_centered\");\n"];
  pr["const int N                   = Geti(\"boundary_N_extrapolate\");\n"];
  pr["const int order_dissipation   = Geti(\"order_dissipation\");\n"];
  pr["const double dissfactor       = get_dissipation_factor(level);\n"];

  pr["const double kappa1      = Getd(\"z4_kappa1\");\n"];
  pr["const double kappa2      = Getd(\"z4_kappa2\");\n"];
  pr["const double chiDivFloor = Getd(\"z4_chi_div_floor\");\n"];
  pr["const double chipsipower = Getd(\"z4_chi_psipower\");\n"];
  pr["const double shiftdriver = Getd(\"z4_shiftdriver\") * Getv(\"z4_bc_use_eta\",\"yes\");\n"];
  
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
                    
  pr["const int givehPsi0 = Getv(\"z4_bc_psi0\",\"yes\");\n"];
  pr["const double hPsi0para = Getd(\"z4_bc_psi0_a\");\n"]; 
  pr["const double hPsi0parb = Getd(\"z4_bc_psi0_b\");\n"]; 
  pr["const double hPsi0parc = Getd(\"z4_bc_psi0_c\");\n"]; 

  pr["\n"];
];




(************************************************************************)
(* now we are ready to go *)


<< "../../math/MathToC/TensorEquationsToC.m"

