(* derivative of one time integral *)
(* mth 10/2012 *)
(* sb  11/2012 *)



(* VARIABLES *)

variables = {I2rpsi4,I2ipsi4, dphiI2rpsi4,dphiI2ipsi4,
             xp, yp, zp, rp,rr }



(* COMPUTE IN THIS ORDER *)

tocompute = {

  (* RADIUS *)

  x == xp - x0,
  y == yp - y0,
  z == zp - z0,

  r == (x x + y y + z z)^(1/2),

  
  (* PARTIAL DERIVATIVES *)

  Cinstruction == "if (order_centered == 2 || boundary1away) {",
    dI2rpsi4[a]      == del[a,I2rpsi4],
    dI2ipsi4[a]      == del[a,I2ipsi4],
  Cinstruction == "} else if (order_centered == 4 || boundaryNaway(2)) {",
    dI2rpsi4[a]      == del4[a,I2rpsi4],
    dI2ipsi4[a]      == del4[a,I2ipsi4],
  Cinstruction == "} else if (order_centered == 6 || boundaryNaway(3)) {",
    dI2rpsi4[a]      == del6[a,I2rpsi4],
    dI2ipsi4[a]      == del6[a,I2ipsi4],
  Cinstruction == "} else errorexit(\"order is not implemented yet\");",

  (* do coordinatetransfo if we use shells *)
  
  Cinstruction == "if (useShellsTransfo) {",


    (* RADIUS ON SHELL *)

     r == rp,
    
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
    dI2rpsi4SST[a]     == Jac[b,a] dI2rpsi4[b],
    dI2ipsi4SST[a]     == Jac[b,a] dI2ipsi4[b],
    
    (* now give the tmp derivative the same name again *)
    dI2rpsi4[a]    == dI2rpsi4SST[a],
    dI2ipsi4[a]    == dI2ipsi4SST[a],
    
  Cinstruction == "}",

  (* COORDINATES *)

  costheta == z/r,
  sintheta == (1-(z/r)^2)^(1/2),
  cosphi == x/(r sintheta),
  sinphi == y/(r sintheta),
  
  (* TAKE INT DRVT *)

  dxyzdphi[a] == 0.,
  Cinstruction == "dxyzdphi1 = -r*sinphi*sintheta;",
  Cinstruction == "dxyzdphi2 =  r*cosphi*sintheta;",
  
  dphiI2rpsi4 == dxyzdphi[a] dI2rpsi4[a],
  dphiI2ipsi4 == dxyzdphi[a] dI2ipsi4[a]
  

} 


(**********************************************************************)
(* symmetries *)


  (* NONLOCAL VARIABLES *)



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

functionname = "dphi_int2_psi4_N";

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

  prdecvl[{I2rpsi4, I2ipsi4, dphiI2rpsi4, dphiI2ipsi4}, "u"];


  pr["tL *level = u->level;\n"];
  pr["\n"];
  pr["double oosqrt2 = 1.0/sqrt(2);\n"];
  pr["double sqrt2 = sqrt(2);\n"];
  pr["\n"];
  
  pr["const int order_centered = Geti(\"order_centered\");\n"];
  
  pr["const double x0 = Getd(\"sphere_x0\");\n"];
  pr["const double y0 = Getd(\"sphere_y0\");\n"];
  pr["const double z0 = Getd(\"sphere_z0\");\n"];
  pr["if (!(x0==y0==z0==0.)) errorexit(\"need sphere center at 0\");\n"];
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
  pr["\n"];
];




(************************************************************************)
(* now we are ready to go *)



<< "../../math/MathToC/TensorEquationsToC.m"

