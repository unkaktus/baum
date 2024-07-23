(* grhd_sources.m
   mth 11/09 *)

   
(* compute right hand side of GRHD equations *)


(* variables *)
variables = {  MD,  Mtau,  MS[a],
              nMD, nMtau, nMS[a], 
              pMD, pMtau, pMS[a], 
              g[a,b], K[a,b], alpha, beta[a],
              Madmrho,MadmS[a],MadmSS[a,b],MadmST, detmetric, 
              Mrho,Mepsl,Mp,Mv[a],Mvsqr,
              xp,yp,zp
            }


(* compute in this order *)
tocompute = {

    matrixdetg        == matrixdet[g],
    Cinstruction == " if (((MATTER.USEMASK) && (mask[ijk]>.99)) || matrixdetg<=0. || alpha[ijk]<=0.) {" ,
        
                    
        sqD          == 0,
        sqS[a]       == 0,
        sqT          == 0,
        Madmrho      == 0,
        MadmS[a]     == 0,
        MadmSS[a,b]  == 0,
        MadmST       == 0,
        
                    
    Cinstruction == " } else {",
        
                    
        (* check if alpha<=0 which can cause to troubles *)
        Cinstruction == " if (alpha[ijk]<=ALPHAMIN){ ",
        Cinstruction == " printf(\"alpha = %e < alpha_min = %e\\n\",alpha[ijk],ALPHAMIN);",
        Cinstruction == " printf(\"  x=%e y=%e z=%e\\n\",Ptr(level,\"x\")[ijk],Ptr(level,\"y\")[ijk],Ptr(level,\"z\")[ijk]); ",
        Cinstruction == " if (ijkinsidefinerlevel(box,ijk)==0) printf(\"  point is NOT inside finer box NOR in some symmetry area\\n\"); ",
        Cinstruction == " else printf(\"  point is inside finer box/ in symmetry \\n\");} ",
                
        Cinstruction == " 
        CheckForNANandINF(11,MD[ijk],MS1[ijk],MS2[ijk],MS3[ijk],Mtau[ijk], Mrho[ijk],Mepsl[ijk],Mp[ijk],Mv1[ijk],Mv2[ijk],Mv3[ijk]); ",
        
        
        (* compute metric derivatives *)
        Cinstruction == "if (order_centered == 2 || boundaryNaway(1)) {",
            delg[c,a,b]   == del[c,g[a,b]],
            dalpha[a]     == del[a,alpha],
            dbeta[a,b]    == del[a,beta[b]],
        Cinstruction == "} else if (order_centered == 4 || boundaryNaway(2)) {",
            delg[c,a,b]   == del4[c,g[a,b]],
            dalpha[a]     == del4[a,alpha],
            dbeta[a,b]    == del4[a,beta[b]],
        Cinstruction == "} else if (order_centered == 6 || boundaryNaway(3)) {",
            delg[c,a,b]   == del6[c,g[a,b]],
            dalpha[a]     == del6[a,alpha],
            dbeta[a,b]    == del6[a,beta[b]],
        Cinstruction == "} else errorexit(\"order is not implemented yet\");",
                
          
        (* determinant *)
        matrixdetg        == matrixdet[g],
        detginv           == 1/matrixdetg,
        ginv[a,b]         == detginv matrixinvdet[g,a,b],
        
        sqrtdetgamma      == sqrt[matrixdetg],
        sqrtmdetg         == alpha*sqrtdetgamma,
        

        (* 3-velocity *)
        vup[a]   == Mv[a],
        v[a]     == g[a,b] Mv[b],
        vsq      == vup[a] v[a],
        W        == 1/sqrt[1-vsq],        
    Cinstruction == "if (vsq>=1.) W = GRHD.HRSC_WlorMAX;", 
        W2hrho   == W*W*(Mrho + Mrho*Mepsl + Mp),
        
        
        (* cooling stuff*)         
        Cif == cooling,
        uupt     == W/ alpha,
        uup[a]   == W (vup[a] - beta[a]/alpha),
        betadown[b] == g[a,b] beta[a],
        udown[a] == g[a,b] uup[b] + betadown[a] uupt,
        EOSKoGmo    == (EOSK/EOSGmo), 
        epsth    == Mepsl-EOSKoGmo*Power[Mrho,EOSGmo] ,
        Lambda   == epsth*Mrho/tcool,
        Cif == else,       
        Lambda   == cooling,
        Cif == end,       
        
        
        (* compute Tmunu *)
        Cif == useprim,
            (* compute Tmunu by using primitive variables *)
            vshift[a]== vup[a]-beta[a]/alpha,
            T00      == (W2hrho - Mp)/(alpha*alpha),
            T0i[a]   == W2hrho*vshift[a]/alpha     + Mp*beta[a]/(alpha*alpha),
            Tij[a,b] == W2hrho*vshift[a]*vshift[b] + Mp*(ginv[a,b]-beta[a]*beta[b]/(alpha*alpha)),
            T0j[a]   == W2hrho*v[a]/alpha,
        Cif == else,
            (* compute Tmunu by using conservative variables *)
            conD     == MD/sqrtdetgamma,
            conS[a]  == MS[a]/sqrtdetgamma,
            conSi[a] == ginv[a,b] conS[b],
            conT     == Mtau/sqrtdetgamma,
            T00      == (conD+conT)/(alpha*alpha),
            T0i[a]   == conSi[a]/alpha - beta[a]*T00,
            Tij[a,b] == beta[a]*beta[b]*T00 + conSi[a]*vup[b] - 
                        (conSi[a]*beta[b]+conSi[b]*beta[a])/alpha +
                        Mp*ginv[a,b],
            T0j[a]   == conS[a]/alpha,
        Cif == end,
        
        
        (* compute source terms *)
        sqD              == 0.,
        sqT              == sqrtmdetg*( 
                            T00*(beta[a]*beta[b]*K[a,b] - beta[a]*dalpha[a]) +
                            T0i[a]*(2*beta[b]*K[a,b]-dalpha[a]) + 
                            Tij[a,b]*K[a,b]-alpha*alpha*uupt*Lambda),
        sqS[c]          == sqrtmdetg*(
                            T00*(0.5*beta[a]*beta[b]*delg[c,a,b]-alpha*dalpha[c]) +
                            T0i[a]*beta[b]*delg[c,a,b] +
                            T0j[a]*dbeta[c,a] +
                            Tij[a,b]*0.5*delg[c,a,b]-alpha*udown[c]*Lambda),
                            
                            
        (* compute ADM vars *)
        (* see Alcubierre book p.248 (7.3.28-30) *)
        Madmrho          == (W2hrho - Mp),
        MadmS[a]         == (W2hrho*vup[a]),
        MadmSS[a,b]      == (W2hrho*v[a]*v[b]   + g[a,b]*Mp),
        MadmST           == (W2hrho*v[a]*vup[a] + 3     *Mp),
        (* Alt. computation, dookie and bam show different behaviour for collapse *)
          (* 
          Madmrho          == T00*(alpha*alpha),
          MadmS[a]         == alpha*g[a,b]*T0i[b],
          MadmSS[a,b]      == g[a,c]*g[b,d]*Tij[c,d], 
          *)

          
        (* set source to 0 if needed *)
        Cinstruction == "if (SetSourceToZero) {",
                              sqD         == 0,
                              sqT         == 0,
                              sqS[a]      == 0,
                              Cinstruction == "}",

	      (* additional hack for ATM *) (* necessary? *) 
        Cinstruction == "if (Mrho[ijk]<GRHD.ATM_RHOATM) {",
            sqS[a]      == 0,
        Cinstruction == "}",

        (* set Tau source = 0, if cold EoS is used Tau rhs = 0 *)
        Cinstruction == "if (EOS.COLD) {",
            sqT          == 0, 
        Cinstruction == "}", 
            

        Cinstruction == "if (CheckForNANandINF(16, Madmrho[ijk], MadmS1[ijk],MadmS2[ijk],MadmS3[ijk],   MadmSS11[ijk],MadmSS12[ijk],MadmSS13[ijk],MadmSS22[ijk],MadmSS23[ijk],MadmSS33[ijk], MadmST[ijk], sqD,sqT,sqS1,sqS2,sqS3) ) {",
            sqD          == 0,
            sqT          == 0,
            sqS[a]       == 0,
            Madmrho      == 0,
            MadmS[a]     == 0,
            MadmSS[a,b]  == 0,
            MadmST       == 0,
        Cinstruction == "}",

          
    Cinstruction == "}",
    
        
    (* add to rhs, source term are always added first to rhs, fluxes later  *)
    Cif == addlinear, 
        nMD          == nMD    + c sqD,
        nMtau        == nMtau  + c sqT,
        nMS[a]       == nMS[a] + c sqS[a],
    Cif == else,
        nMD          == nMD    +   sqD,
        nMtau        == nMtau  +   sqT,
        nMS[a]       == nMS[a] +   sqS[a],
    Cif == end,
    
    
    Cinstruction == "if (CheckForNANandINF(16, Madmrho[ijk], MadmS1[ijk],MadmS2[ijk],MadmS3[ijk], MadmSS11[ijk],MadmSS12[ijk],MadmSS13[ijk],MadmSS22[ijk],MadmSS23[ijk],MadmSS33[ijk], MadmST[ijk], sqD,sqT,sqS1,sqS2,sqS3)) { ",
    Cinstruction == " printf(\"problem with nan's inside con2source\\n\"); ",
    Cinstruction == " printf(\"  x=%e y=%e z=%e\\n\",Ptr(level,\"x\")[ijk],Ptr(level,\"y\")[ijk],Ptr(level,\"z\")[ijk]); ",
    Cinstruction == " printf(\"  sD=%e sT=%e sSx=%e sSy=%e sSz=%e\\n\",sqD,sqT,sqS1,sqS2,sqS3); ",
    Cinstruction == " printf(\"  %e %e %e %e   %e %e %e  %e\\n\",v1,v2,v3,vsq,  W2hrho,g11[ijk],Mp[ijk],  MadmSS12[ijk]); ",
    Cinstruction == " if (ijkinsidefinerlevel(box,ijk)==0) printf(\"  point is NOT inside finer box NOR in some symmetry area\\n\"); ",
    Cinstruction == " else printf(\"  point is inside finer box/ in symmetry \\n\"); ",
    Cinstruction == " } "
    
}



(* symmetries *)
yyy[a_,b_] :=  yyy[b,a] /; !OrderedQ[{a,b}]

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
  pr["#define SetSourceToZero 0\n\n"];

  pr["\n\n\n"];
];

functionname = "grhd_sources";
functionarguments = "tVarList *unew, tVarList *upre, double c, tVarList *ucur, tVarList *primVars, tVarList *otherVars, tVarList *matterADMVars";

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
  pr["const int order_centered = Geti(\"order_centered\");\n"];
  pr["const int useprim = Getv(\"grhd_source_computation\",\"prim\");\n"];
  pr["double cooling = Getd(\"grhd_cooling\");\n"];
  pr["double tcool = Getd(\"grhd_cooling_time\");\n"];
  pr["double EOSK = EOS.K;\n"];
  pr["double EOSGmo = EOS.GAMMAMO;\n"]; 
  
  pr["\n"];
  
  prdecvlNR[{ MD,  Mtau , MS[a]}, "ucur", "MATTER.INDX_VAR_q"];
  prdecvlNR[{nMD, nMtau, nMS[a]}, "unew", "MATTER.INDX_VAR_q"];
  prdecvlNR[{pMD, pMtau, pMS[a]}, "upre", "MATTER.INDX_VAR_q"];
  
  prdecvlNR[{alpha},   "ucur", "MATTER.INDX_VAR_alp"];
  prdecvlNR[{beta[a]}, "ucur", "MATTER.INDX_VAR_bet"];
  
  prdecvarname[{g[a,b]}, "adm_gxx"];
  prdecvarname[{K[a,b]}, "adm_Kxx"];
  
  prdecvarname[{mask}, "matter_mask"];
    
  pr["const double *xp = level->v[Ind(\"x\")];\n"];
  pr["const double *yp = level->v[Ind(\"y\")];\n"];
  pr["const double *zp = level->v[Ind(\"z\")];\n"];
  
  prdecvl[{Mrho,Mepsl,Mv[a],Mp,Mvsqr,detmetric}, "primVars"];
  prdecvl[{Madmrho,MadmS[a],MadmSS[a,b],MadmST}, "matterADMVars"];
  
  pr["\n"];
  
  
];




(************************************************************************)
(* now we are ready to go *)

(* assume that we are working on a box *)
BoxMode = True;


<< "../../math/MathToC/TensorEquationsToC.m"
