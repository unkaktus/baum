(* Z4d_algcon.m 
   Wolfgang Tichy 8/2009
*)

(* compute trA and detg, enforce algebraic constraints *)


(* variables *)
variables = {g[a,b], A[a,b],  trAstore, detgstore}


(* compute in this order *)
tocompute = {

  (* determinant of metric *)
  detg == matrixdet[g],
  Cinstruction == "if (detg<=0.) { PROBLEM; /*errorexit(\" detg <= 0.\");*/ detg=1.;} ",
  detginv == 1/detg,
  
  (* enforce detg = 1 *)
  Cif == normalizedetg,
    epsilon == detg-1.,
    Cinstruction == "if (fabs(epsilon)<1e-12) aux = 1.-oot*epsilon; else ",
    aux == detg^(-oot),
    g[a,b] == aux g[a,b],
  Cif == end,

  (* trace of extrinsic curvature using new metric if rescaled *)
  detg == matrixdet[g],
  Cinstruction == "if (detg<=0.) { PROBLEM; /*errorexit(\" detg <= 0.\");*/ detg=1.;} ",
  detginv == 1/detg,
  trA == detginv matrixinvdet[g,a,b] A[a,b],

  (* enforce trA = 0 *)
  Cif == subtractA,
    aux == -trA oot,
    A[a,b] == A[a,b] + aux g[a,b],
  Cif == end,

  trA == detginv matrixinvdet[g,a,b] A[a,b],

  (* store result if storage is enabled *)
  Cinstruction == "if (store) { ",
    detgstore  == detg,
    trAstore   == trA,
  Cinstruction == "} "
}




(* symmetries *)
 g[a_,b_] :=  g[b,a] /; !OrderedQ[{a,b}]
 A[a_,b_] :=  A[b,a] /; !OrderedQ[{a,b}]
ginv[a_,b_] := ginv[b,a] /; !OrderedQ[{a,b}]
Ainv[a_,b_] := Ainv[b,a] /; !OrderedQ[{a,b}]




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
  pr["#define Log(x)     log((double) (x)))\n"];
  pr["#define pow2(x)    ((x)*(x))\n"];
  pr["#define pow2inv(x) (1.0/((x)*(x)))\n\n"];
  pr["#define Cal(x,y,z) ((x)?(y):(z))\n\n"];
  pr["#define PROBLEM    printf(\"  %d pts away from boundary\\n\",boundaryaway(6)); \\\n"];
  pr["                   printf(\"  %e %e %e \\n  detg=%e,trA=%e \\n\", \\\n"];
  pr["                   Ptr(level,\"x\")[ijk],Ptr(level,\"y\")[ijk],Ptr(level,\"z\")[ijk], \\\n"];
  pr["                   detg, trA); stop++ \n"];
  
  pr["\n\n\n"];
];

functionname = "z4_algcon";

functionbegin     = "bampi_openmp_start\n";

functionarguments = "tVarList *u";

functionloop = "int stop=0; forallpoints_ijk_openmp(level) {";

endfunctionloop = "} endfor_ijk_openmp;/* if (stop) errorexit(\"\"); */ /* loop ccc */";

endfunction       = "bampi_openmp_stop\n";


(* custom variable declaration
   we have to translate between the tensor names and the C variables  
   we may need special variables, for example for parameters
*)
variabledeclarations[] := Module[{},

  prdecvl[{g[a,b], A[a,b], trAstore, detgstore}, "u"];

  pr["tL *level = u->level;\n"];
  pr["\n"];

  pr["const int subtractA      = Getv(\"z4_subtractA\", \"yes\");\n"];
  pr["const int normalizedetg  = Getv(\"z4_normalizedetg\", \"yes\");\n"];
  pr["const double oot = 1.0/3.0;\n"];
  pr["const int store          = Getv(\"z4_register_algcon\", \"store\");\n"];
];    




(************************************************************************)
(* now we are ready to go *)

<< "../../math/MathToC/TensorEquationsToC.m"
