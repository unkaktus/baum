(* BACKUPFILE
 radiative on boxes using centered stencils option 7 *)

(************************************************************************)
(* Second order Sommerfeld boundary conditions on boxes *)
(************************************************************************)

(* background.m mth 02/10 *)
(* background.m dh  08/10 *)
(* background.m mth 05/11   -> now only shells *)
(* Z4d_boundary.m mth 06/11 -> moved from boundary to Z4d *)
(* Z4d_boundary.m dh 10/11  -> Conditions without linearization *)
(* Z4c_boundary_box_rad.m sb/dh 11/15  -> radiative on boxes using centered stencils *)

(* variables *)
variables = { g[a,b],  chi,  A[a,b],  Khat,  G[a],  Theta,  alpha,  beta[a],  B[a], 
	      ng[a,b], nchi, nA[a,b], nKhat, nG[a], nTheta, nalpha, nbeta[a], nB[a],
	      pg[a,b], pchi, pA[a,b], pKhat, pG[a], pTheta, palpha, pbeta[a], pB[a],
	      xp,yp,zp
            }
  
(* compute in this order *)
tocompute = {

  (* Loop on all innerpoints, but want only those 'N_point_away' from boundary *)
  (* Is there a better way to do this? *)
  Cinstruction == "if (!(boundaryNaway(N_points_away))) continue;",

  Cinstruction == "if (CheckForNANandINF(25, 
    g11[ijk],g12[ijk],g13[ijk],g22[ijk],g23[ijk],g33[ijk],
    A11[ijk],A12[ijk],A13[ijk],A22[ijk],A23[ijk],A33[ijk], 
    G1[ijk],G2[ijk],G3[ijk], Khat[ijk],chi[ijk], Theta[ijk],
    alpha[ijk],beta1[ijk],beta2[ijk],beta3[ijk],B1[ijk],B2[ijk],B3[ijk])) {
    printf(\"problem with vars in Z4d_boundary.m\\n\");
    printf(\"x=%2.5e, y=%2.5e, z=%2.5e\\n\",xp[ijk],yp[ijk],zp[ijk]);
  }\n",

(************************************************************************)
(* ---- Compute derivatives, 2nd order centered, assume one point always available  ---- *)
(************************************************************************)
  
  (*
        da[a]       == del[a,alpha],
        db[a,b]     == del[a,beta[b]],
        dchi[a]     == del[a,chi],
        dg[c,a,b]   == del[c,g[a,b]],
  *)
        dKhat[a]    == del[a,Khat],
        dA[a,b,c]   == del[a,A[b,c]],
        dG[a,b]     == del[a,G[b]],
        dTheta[a]   == del[a,Theta],

(************************************************************************)
(* ----  ---- *)
(************************************************************************)

    (* ---- Define normal vector ---- *)
    r == 0,
    x == 0,
    y == 0,
    z == 0,
    sup[a] == 0,    

    Cinstruction == "x=xp[ijk]; y=yp[ijk]; z=zp[ijk];",
    Cinstruction == "r=sqrt(x*x + y*y + z*z);",
    Cinstruction == "sup1 = x/r;",
    Cinstruction == "sup2 = y/r;",
    Cinstruction == "sup3 = z/r;",

    rTheta == - sup[a] dTheta[a] - Theta/r,
    rKhat == - sqrt[2] ( sup[a] dKhat[a] + (Khat)/r ),
    rG[a] == - sup[b] dG[b,a] - (G[a])/r,
    rA[a,b] == - sup[c] dA[c,a,b] - (A[a,b])/r,

  (*
    rchi == - sqrt[2] ( sup[a] dchi[a] + (chi-1)/r ),
    ralpha == - sqrt[2] ( sup[a] da[a] + (alpha-1)/r ),
    rg[a,b] == - sup[c] dg[c,a,b] - (g[a,b] - delta[a,b])/r,
    rbeta[a] == - sup[b] db[b,a] - (beta[a])/r,
  *)

  (* check RHS *)
  Cinstruction == "if (CheckForNANandINF(11, 
    rA11,rA12,rA13,rA22,rA23,rA33, 
    rG1,rG2,rG3,rKhat,rTheta)) {
    printf(\"nans in RHS in z4_boundary_box_rad.m\\n\");
    printf(\"x=%2.5e, y=%2.5e, z=%2.5e\\n\",xp[ijk],yp[ijk],zp[ijk]);
  }",

(************************************************************************)
(* ---- Dissipation ---- *)
(************************************************************************)

  (* add dissipation terms we should test which terms are really needed *)
  (*
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
   *)

(************************************************************************)
(* ---  update the RHS---- *)
(************************************************************************)

   Cif == addlinear, 
       nKhat       == pKhat   + c rKhat,
       nA[i,j]     == pA[i,j] + c rA[i,j],
       nG[i]       == pG[i]   + c rG[i],
       nTheta      == pTheta  + c rTheta,
       (* 
       nchi        == pchi     + c rchi,
       nalpha      == palpha   + c ralpha, 
       ng[i,j]     == pg[i,j]  + c rg[i,j],
       nbeta[i]    == pbeta[i] + c rbeta[i], 
	*) 
  Cif == else,
       nKhat       == rKhat,
       nA[i,j]     == rA[i,j],
       nG[i]       == rG[i],
       nTheta      == rTheta,
       (* 
       nchi        == rchi,
       nalpha      == ralpha,
       ng[i,j]     == rg[i,j],
       nbeta[i]    == rbeta[i],
	*) 
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
  
  dg[c_,a_,b_]    := dg[c,b,a] /; !OrderedQ[{a,b}]
  dA[c_,a_,b_]    := dA[c,b,a] /; !OrderedQ[{a,b}]

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

functionname      = "z4_boundary_box_rad";

functionarguments = "tVarList *unew, tVarList *upre, double c, tVarList *ucur";

(* functionbegin     = "bampi_openmp_start\n"; *)

functionloop      = "forinnerpoints_ijk_openmp(level) {";

endfunctionloop   = "} endfor_ijk_openmp; /* loop i, j, k */";

(* endfunction       = "bampi_openmp_stop\n"; *)




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
  
  (*
  pr["const int order               = Geti(\"order_centered\");\n"];
  pr["const int N                   = Geti(\"boundary_N_extrapolate\");\n"];
  *)

  (* where to apply boundaries*)
  pr["const int N_points_away = 1;\n"];

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
  
  pr["\n"];
];




(************************************************************************)
(* now we are ready to go *)


<< "../../math/MathToC/TensorEquationsToC.m"

