(* adm_mass_integrand_csi.m 
   mth 2011 *)

(* compute integrand of the ADM mass *)
(* using simple formula that involves only the conformal factor *)
(* See notes in Invariants/doc                                  *)

(* variables *)
variables = {x[a], chi, integrand,
             xp,yp,zp,rp,rr}



(* compute in this order *)
tocompute = {

  delchi[a] == del[a,chi],


  (* do coordinatetransfo if we use shells *)
  Cinstruction == "if (useShellsTransfo) {",
  
    (* derivatives for fisheye coords ... this is a speedup *)
    dRdr         == DRDr,
    ddRdr        == DDRDr,
    
    (* define jacobians for each box *)
    Cinstruction == "if (nbox<=1 ) {",
      Jac[a,b]    == JacobB01[a,b],
    Cinstruction == "} else if (nbox<=3) {",
      Jac[a,b]    == JacobB23[a,b],
    Cinstruction == "} else if (nbox<=5) {",
      Jac[a,b]    == JacobB45[a,b],
    Cinstruction == "}",
    
    (* now do the transformation and store in a tmp name *)
    delchiSST[a]      == Jac[b,a] delchi[b],
    
    (* now give the tmp derivative the same name again *)
    delchi[a]         == delchiSST[a],
    
  Cinstruction == "}",




  (* R *)
  R == sqrt [delta[a,b] x[a] x[b]],
    
  (* Normal *)
  NN[a] == x[a]/R,

  (* Integrand of the ADM mass *)
  integrand ==  - 1.0/2.0/PI (-1/4*chi^(-5/4)) delchi[c] NN[c] R*R

}

(* symmetries *)
g[a_,b_]       := g[b,a] /; !OrderedQ[{a,b}]
ginv[a_,b_]    := ginv[b,a] /; !OrderedQ[{a,b}]
delg[c_,a_,b_] := delg[c,b,a] /; !OrderedQ[{a,b}]
delgSST[c_,a_,b_] := delgSST[c,b,a] /; !OrderedQ[{a,b}]


(************************************************************************)
(* information for C output *)

(* information about the function 
   the C file will be functionname.c 
*)
includeANDdefine[] := Module[{},

  pr["#include \"bam.h\"\n"];
  pr["#include \"ADM_mass.h\"\n"];
  pr["\n"];
  pr["#define Power(x,y) (pow((double) (x), (double) (y)))\n"];
  pr["#define pow2(x)    ((x)*(x))\n"];
  pr["#define pow3(x)    ((x)*(x)*(x))\n"];
  pr["#define Sqrt(x)    sqrt(x)\n"];  
  pr["#define Tanh(x)    tanh(x)\n"];
  pr["#define Sech(x)    (1/cosh(x))\n"];
  pr["#define PI  3.1415926535897932\n\n"];
  pr["\n\n\n"];

];

functionname = "adm_mass_integrand_chi";

functionarguments = "tL *level, int i_x, int i_chi, int i_integrand";

functionloop = "forinnerpoints_ijk(level) {";
endfunctionloop = "} endfor_ijk; /* loop i, j, k */";


(* custom variable declaration
   we have to translate between the tensor names and the C variables  
   we may need special variables, for example for parameters
*)
variabledeclarations[] := Module[{},

  prdecvar[{x[a]},      "i_x"];
  prdecvar[{chi},       "i_chi"];
  prdecvar[{integrand}, "i_integrand"];
  
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
  
];




(************************************************************************)
(* now we are ready to go *)

(* assume that we are working on a box *)
BoxMode = True;


<< "../../math/MathToC/TensorEquationsToC.m"
