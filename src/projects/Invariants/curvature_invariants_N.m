(* curvature_invariants_N.m *)
(* Jose Gonzalez 01.06 *)
(* mth 01/2011 *)



(* VARIABLES *)

variables = {g[a,b], K[a,b], alpha, beta[a], 
             xp, yp, zp, rp,rr,
             rpsi0, ipsi0, rpsi1, ipsi1, 
             rpsi2, ipsi2, rpsi3, ipsi3, rpsi4, ipsi4,
             rI, iI, rJ, iJ, rS, iS, Csqr}



(* COMPUTE IN THIS ORDER *)

tocompute = {

  (* PARTIAL DERIVATIVES *)

  Cinstruction == "if (order_centered == 2 || boundary1away) {",
    dg[c,a,b]      == del[c,g[a,b]],
    ddg[a,b,c,d]   == deldel[a,b,g[c,d]],
    dK[a,b,c]      == del[a, K[b,c]],
  Cinstruction == "} else if (order_centered == 4 || boundaryNaway(2)) {",
    dg[c,a,b]      == del4[c,g[a,b]],
    ddg[a,b,c,d]   == deldel4[a,b,g[c,d]],
    dK[a,b,c]      == del4[a, K[b,c]],
  Cinstruction == "} else if (order_centered == 6 || boundaryNaway(3)) {",
    dg[c,a,b]      == del6[c,g[a,b]],
    ddg[a,b,c,d]   == deldel6[a,b,g[c,d]],
    dK[a,b,c]       == del6[a, K[b,c]],
  Cinstruction == "} else errorexit(\"order is not implemented yet\");",

  (* do coordinatetransfo if we use shells *)
  
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
    dgSST[a,b,c]     == Jac[d,a] dg[d,b,c],
    ddgSST[a,b,c,d]  == DJac[e,a,b] dg[e,c,d] + Jac[e,a] Jac[f,b] ddg[e,f,c,d],
    dKSST[a,b,c]     == Jac[d,a] dK[d,b,c],
    
    (* now give the tmp derivative the same name again *)
    dg[a,b,c]    == dgSST[a,b,c],
    ddg[a,b,c,d] == ddgSST[a,b,c,d],
    dK[a,b,c]    == dKSST[a,b,c],
    
  Cinstruction == "}",
  


  (* PHYSICAL METRIC AND ITS INVERSE *)
  
   ginvdet[a,b] == matrixinvdet[g,a,b],
   detg == g[1,c] ginvdet[1,c],
   ginv[a,b] == ginvdet[a,b] / detg, 


  (* CHRISTOFFEL *) 

   gammado[c,a,b] == 1/2 (dg[b,a,c] + dg[a,b,c] - dg[c,a,b]),
   gamma[c,a,b] == ginv[c,d] gammado[d,a,b],


  (* RIEMANN AND RICCI: Using MTW conventions *)

   R[a,b] == ginv[c,d] (1/2(
             - ddg[c,d,a,b] - ddg[a,b,c,d] + ddg[a,c,b,d] + ddg[b,c,a,d])
             + gamma[e,a,c] gammado[e,b,d] - gamma[e,a,b] gammado[e,c,d]),
   R == ginv[a,b] R[a,b],
   R[a,b,c,d] == g[a,c] R[b,d] + g[b,d] R[a,c] - 1/2 R g[a,c] g[b,d] 
               - g[a,d] R[b,c] - g[b,c] R[a,d] + 1/2 R g[a,d] g[b,c],


  (* TRANSLATE EXTRINSIC CURVATURE TO LOCAL VARIABLES *)
   trK == ginv[c,d] K[c,d],


  (* 4-RIEMANN FROM THE GAUSS-CODAZZI EQUATIONS *)

   Riemm[a,b,c,d] == R[a,b,c,d] + K[a,c] K[b,d] - K[a,d] K[b,c],
   (* Possible BUG: I do not think, that agrees with my definition *)
   Riemm[a,b,c] == - (dK[c,a,b] - dK[b,a,c] + gamma[d,a,b] K[c,d] - gamma[d,a,c] K[b,d] ),
   Riemm[a,b] == R[a,b] + trK K[a,b] - ginv[c,d] K[a,c] K[d,b],

   tfactor == 1/2^(1/2),


  (* ELECTRIC AND MAGNETIC PARTS *)
  (* Added Kretschmann invariant, Mark Hannam, August 2006 *)

   Cinstruction == "if (Csqr) {",

     Sqrtdetg == (detg)^(1/2),

     E[a,b] == Riemm[a,b],
     B[a,b] == -1/2/Sqrtdetg g[a,c] epsmatrix[c,e,f] Riemm[b,e,f],

     Csqr == (ginv[a,c] ginv[b,d] E[a,b] E[c,d] - ginv[a,c] ginv[b,d] B[a,b] B[c,d]),

   Cinstruction == "}",

  (* coordinates *)

   x == xp - x0,
   y == yp - y0,
   z == zp - z0,

   r == (x x + y y + z z)^(1/2),

   costheta == z/r,
   sintheta == (1-(z/r)^2)^(1/2),
   cosphi == x/(r sintheta),
   sinphi == y/(r sintheta),

   Cif == kinnersley,

     (* KINNERSLEY TETRAD *)
 
     n0 == 1/2,
     n[a] == -1/2 (delta[a,1] sintheta cosphi + delta[a,2] sintheta sinphi + delta[a,3] costheta),
     rm[a] == tfactor (delta[a,1] costheta cosphi + delta[a,2] costheta sinphi - delta[a,3] sintheta),
     im[a] == tfactor (-delta[a,1] sinphi + delta[a,2] cosphi),

   Cif == else,

     (* NEW TETRAD USING GRAM_SCHMIDT PROCEDURE *)

     v1[a] == -y delta[a,1] + x delta[a,2],
     v2[a] ==  x delta[a,1] + y delta[a,2] + z delta[a,3], 
     v3[a] == detg^(1/2) ginv[a,d] epsmatrix[d,b,c] v1[b] v2[c],

     (* BUG: in v2[a]=...: first subtract the v1-part, then take the norm! *)
     (*      same for v3[a]=... *)
     w11 == v1[a] v1[b] g[a,b],
     v1[a] == v1[a]/(w11)^(1/2),
     w12 == v1[a] v2[b] g[a,b],
     v2[a] == (v2[a] - v1[a] w12),
     w22 == v2[a] v2[b] g[a,b],
     v2[a] == v2[a]/(w22)^(1/2),
     w13 == v1[a] v3[b] g[a,b],
     w23 == v2[a] v3[b] g[a,b],
     v3[a] == (v3[a] - v1[a] w13 - v2[a] w23),
     w33 == v3[a] v3[b] g[a,b],
     v3[a] == v3[a]/(w33)^(1/2),

     vecr[a] == v2[a],
     vect[a] == v3[a],
     vecp[a] == v1[a],
     vecu0 == 1/alpha,
     vecu[a] == -beta[a]/alpha,   

     n0 == tfactor vecu0,
     n[a] == tfactor (vecu[a] - vecr[a]),
     rm[a] == tfactor vect[a],
     im[a] == tfactor vecp[a],

   Cif == end,

   rmb[a] == rm[a],
   imb[a] == -im[a],

  (* PSI4 *)
  (* Possible BUG: here we split off projections of Riemann onto
                   the time component. We should split off the
                   projections onto n, however *)

   Cif == gausscodacciflat,

     rpsi4 == -(Riemm[a,b,c,d] n[a] n[c] + 2 Riemm[b,c,d] n0 n[c]
              + Riemm[b,d] n0 n0) (rmb[b]rmb[d] - imb[b]imb[d]),
     ipsi4 == -(Riemm[a,b,c,d] n[a] n[c] + 2 Riemm[b,c,d] n0 n[c]
              + Riemm[b,d] n0 n0) (rmb[b]imb[d] + imb[b]rmb[d]),

   Cif == else,

     rpsi4 == -(Riemm[a,b,c,d] vecr[a] vecr[c] - 2 Riemm[b,c,d] vecr[c]
              + Riemm[b,d]) (rmb[b]rmb[d] - imb[b]imb[d]) / 2,
     ipsi4 == -(Riemm[a,b,c,d] vecr[a] vecr[c] - 2 Riemm[b,c,d] vecr[c]
              + Riemm[b,d]) (rmb[b]imb[d] + imb[b]rmb[d]) / 2,

   Cif == end


} 


(**********************************************************************)
(* symmetries *)


  (* NONLOCAL VARIABLES *)

   g[a_,b_] := g[b,a] /; !OrderedQ[{a,b}]
   K[a_,b_] := K[b,a] /; !OrderedQ[{a,b}]


  (* LOCAL VARIABLES *)
 
   deldel[a_,b_,c_] := deldel[b,a,c] /; !OrderedQ[{a,b}]

  
   ginvdet[a_,b_] := ginvdet[b,a] /; !OrderedQ[{a,b}]


   dg[c_,a_,b_] := dg[c,b,a] /; !OrderedQ[{a,b}]

   ddg[a_,b_,c_,d_] := ddg[b,a,c,d] /; !OrderedQ[{a,b}]
   ddg[a_,b_,c_,d_] := ddg[a,b,d,c] /; !OrderedQ[{c,d}]
   
   gamma[c_,a_,b_] := gamma[c,b,a] /; !OrderedQ[{a,b}]
   gammado[c_,a_,b_] := gammado[c,b,a] /; !OrderedQ[{a,b}]
   (* hidden symmetries for contraction! careful, not general *)
   gamma/: gamma[a_, b_, c_] gammado[a_, d_, e_] :=
   gammado[a, b, c] gamma[a, d, e] /; !OrderedQ[{{b,c},{d,e}}]
  
   dK[c_,a_,b_] := dK[c,b,a] /; !OrderedQ[{a,b}]
   
   R[a_,a_,c_,d_] := 0
   R[a_,b_,c_,c_] := 0
   R[a_,b_,c_,d_] := - R[b,a,c,d] /; !OrderedQ[{a,b}]
   R[a_,b_,c_,d_] := - R[a,b,d,c] /; !OrderedQ[{c,d}]
   R[a_,b_,c_,d_] := R[c,d,a,b] /; !OrderedQ[{{a,b},{c,d}}] 

   R[a_,b_] := R[b,a] /; !OrderedQ[{a,b}]

   Riemm[a_,a_,c_,d_] := 0
   Riemm[a_,b_,c_,c_] := 0
   Riemm[a_,b_,c_,d_] := - Riemm[b,a,c,d] /; !OrderedQ[{a,b}]
   Riemm[a_,b_,c_,d_] := - Riemm[a,b,d,c] /; !OrderedQ[{c,d}]
   Riemm[a_,b_,c_,d_] := Riemm[c,d,a,b] /; !OrderedQ[{{a,b},{c,d}}]   
   
   Riemm[a_,b_,b_] := 0
   Riemm[a_,b_,c_] := - Riemm[a,c,b] /; !OrderedQ[{b,c}]

   Riemm[a_,b_] := Riemm[b,a] /; !OrderedQ[{a,b}]

   (* The totally antisymmetric matrix *)
   (* Check: do eps*eps terms work? *)
   epsmatrix[a_Integer,b_Integer,c_Integer] := (b-a)(c-b)(c-a)/2
   epsmatrix[a_,b_,c_] := -epsmatrix[b,a,c] /; !OrderedQ[{a,b}]
   epsmatrix[a_,b_,c_] := -epsmatrix[a,c,b] /; !OrderedQ[{b,c}]


  dgbSST[c_,a_,b_]     := dgbSST[c,b,a] /; !OrderedQ[{a,b}]
  ddgbSST[a_,b_,c_,d_] := ddgbSST[b,a,c,d] /; !OrderedQ[{a,b}]
  ddgbSST[a_,b_,c_,d_] := ddgbSST[a,b,d,c] /; !OrderedQ[{c,d}]
  dKbSST[c_,a_,b_]     := dKbSST[c,b,a] /; !OrderedQ[{a,b}]

  DJac[a_,b_,c_]       := DJac[a,c,b] /; !OrderedQ[{b,c}]


(************************************************************************)
(************************************************************************)
(* information for C output *)

(* information about the function 
   the C file will be functionname.c 
*)
includeANDdefine[] := Module[{},
  pr["#include \"bam.h\"\n"];
  pr["#include \"Invariants.h\"\n"];
  pr["\n"];
  pr["#define Power(x,y) pow((double) (x), (double) (y))\n"];
  pr["#define Sqrt(x)    sqrt((double) (x))\n"];
  pr["#define Log(x)     log((double) (x))\n"];
  pr["#define pow2(x)    ((x)*(x))\n"];
  pr["#define pow4(x)    ((x)*(x)*(x)*(x))\n"];
  pr["#define pow2inv(x) (1.0/((x)*(x)))\n\n"];
  pr["#define Tanh(x)    tanh(x)\n"];
  pr["#define Sech(x)    (1/cosh(x))\n"];
  pr["#define Cal(x,y,z) ((x)?(y):(z))\n\n"];
  pr["\n\n\n"];
];

functionname = "curvature_invariants_N";

functionarguments = "tVarList *u";

functionbegin     = "bampi_openmp_start\n";

functionloop      = "forinnerpoints_ijk_openmp(level) {";

endfunctionloop   = "} endfor_ijk_openmp; /* loop i, j, k */";

endfunction       = "bampi_openmp_stop\n";




(* custom variable declaration
   we have to translate between the tensor names and the C variables  
   we may need special variables, for example for parameters
*)
variabledeclarations[] := Module[{},

  prdecvl[{g[a,b], K[a,b], alpha, beta[a], xp, yp, zp, 
        (* rpsi0, ipsi0,
           rpsi1, ipsi1, 
           rpsi2, ipsi2,
           rpsi3, ipsi3, *)
           rpsi4, ipsi4,
        (* rI, iI, rJ, iJ, rS, iS *)
           Csqr
          }, "u"];



  pr["tL *level = u->level;\n"];
  pr["\n"];
  pr["double oosqrt2 = 1.0/sqrt(2);\n"];
  pr["double sqrt2 = sqrt(2);\n"];
  pr["\n"];
  pr["const int kinnersley = Getv(\"invariants_tetrad\",\"kinnersley\");\n"];
  pr["const int gausscodacciflat = Getv(\"gauss_codacci_mainardi\", \"flat\");\n"];
  
  pr["const int order_centered = Geti(\"order_centered\");\n"];
  
  pr["const double x0 = Getd(\"sphere_x0\");\n"];
  pr["const double y0 = Getd(\"sphere_y0\");\n"];
  pr["const double z0 = Getd(\"sphere_z0\");\n"];
  pr["\n"];
  pr["const int useShellsTransfo = (Getv(\"grid\",\"shells\") && level->l==0);\n"];
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

