(* adm_mass_integrand.m 
   Jose Gonzalez 10/04 *)

(* compute integrand of the ADM mass *)


(* variables *)
variables = {x[a], g[a,b], psi, dpop[a], integrand,
             xp, yp, zp,rp,rr }



(* compute in this order *)
tocompute = {

  (* partial derivatives *)
  delg[c,a,b] == del[c,g[a,b]],
  
  
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
    delgSST[a,b,c]      == Jac[d,a] delg[d,b,c],
    
    (* now give the tmp derivative the same name again *)
    delg[a,b,c]         == delgSST[a,b,c],
    
  Cinstruction == "}",

  (* Determinant of 3-metric *)
  det == matrixdet[g],

  (* inverse physical metric *)
  detginvf == 1/(matrixdet[g]),
  ginv[a,b] == detginvf matrixinvdet[g,a,b],

  (* R *)
  R == sqrt [delta[a,b] x[a] x[b]],

  (* Normal *)
  NN[a] == x[a]/R,

  (* Integrand of the ADM mass *)
  integrand ==  1.0/16.0/PI sqrt[det] ginv[a,b] (delg[b,a,c] - delg[c,a,b]) NN[c] R*R

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
  pr["#define Sqrt(x)    sqrt(x)\n"];
  pr["#define Log(x)     log((double) (x)))\n"];
  pr["#define pow2(x)    ((x)*(x))\n"];
  pr["#define pow3(x)    ((x)*(x)*(x))\n"];
  pr["#define pow4(x)    ((x)*(x)*(x)*(x))\n"];
  pr["#define pow2inv(x) (1.0/((x)*(x)))\n\n"]; 
  pr["#define Tanh(x)    tanh(x)\n"];
  pr["#define Sech(x)    (1/cosh(x))\n"];
  pr["#define Cal(x,y,z) ((x)?(y):(z))\n\n"];
  pr["#define PI  3.1415926535897932\n\n"];
  pr["\n\n\n"];

];

functionname = "adm_mass_integrand";

functionarguments = "tL *level, int i_x, int i_g, int i_integrand";

functionloop = "forinnerpoints_ijk(level) {";
endfunctionloop = "} endfor_ijk; /* loop i, j, k */";


(* custom variable declaration
   we have to translate between the tensor names and the C variables  
   we may need special variables, for example for parameters
*)
variabledeclarations[] := Module[{},

  prdecvar[{x[a]},      "i_x"];
  prdecvar[{g[a,b]},    "i_g"];
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

<< "../../math/MathToC/TensorEquationsToC.m"
