(* ::Package:: *)

(* compute_hydrovars.m *)
(* S Bernuzzi 10.2012 *)


(* VARIABLES *)

variables = {conD, pres, rho, epsl, v[a],
             g[a,b], alpha, beta[a], 
             xp, yp, zp,             
             Dh, Db, Du, ut,
             etot, uesc,
             vor[a], ivor[a], 
	     P[a], Pu[a]}


(* COMPUTE IN THIS ORDER *)

tocompute = {



    (* skip atm points *)

    Cinstruction == " if ((MATTER.USEMASK) && (mask[ijk]>.99)) {" ,

    Db == 0,
    Du == 0,  
    Dh == 0,
    ivor[a] == 0,
    
    Cinstruction == " continue;",
    Cinstruction == " }",
    


    (* 4-velocity component u_t *)
    
    vup[a] == v[a],    
    vsq == vup[a] vup[b] g[a,b],
    
    ooalpha == 1/alpha,
    Cinstruction == "if (alpha[ijk]<=ALPHAMIN){",
    Cinstruction == " printf(\" alpha<%e in hydroanalysis, reset\\n\",ALPHAMIN); ooalpha=1./ALPHAMIN;",    
    Cinstruction == "}",

    W == (1-vsq)^(-1/2),
    Cinstruction == "if (vsq>=1.) W = Wmax;", 

    uup[a] == W (vup[a] - beta[a] ooalpha),
    
    uupt == W ooalpha,

    betadown[b] == g[a,b] beta[a],

    udown[a] == g[a,b] uup[b] + betadown[a] uupt,
    
    h == 1+epsl+pres/rho,
    
    betasq ==  beta[b] betadown[b],
    betaduu == betadown[a] uup[a],

    udownt == (-alpha^2 + betasq ) uupt + betaduu,

    ut == udownt,

    ur     == xp vup[1] + yp vup[2] + zp vup[3],
    
    P[a] == conD h udown[a],    
    
    (* vorticity *)
    (* almost: it requires a derivative, done later *)

    ivor[a] == h udown[a],
    


    (* set various rest-masses *)
    
    Db   == conD,
    Du   == conD,
    Dh   == conD,

    (* compute bound/unbound rest-mass *)
    
    Cinstruction == " if ( (- udownt > 1)  && (ur > 0) ){" ,
    
    Db == 0,

    auupt == alpha uupt,    
    
    etot == conD*(h auupt - pres /(rho auupt)),
   
    uesc == conD epsl, 
   
    Pu[a] == P[a],    
 
    Cinstruction == " } else {",
    
    Du == 0,
    uesc == 0,
    etot == 0,
    Pu[a] == 0,

    Cinstruction == " } "

}


(**********************************************************************)
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
  pr["#include \"hydroanalysis.h\"\n"];
  pr["\n"];
  pr["#define Power(x,y) pow((double) (x), (double) (y))\n"];
  pr["#define Sqrt(x)    sqrt((double) (x))\n"];
  pr["#define Log(x)     log((double) (x))\n"];
  pr["#define pow2(x)    ((x)*(x))\n"];
  pr["#define pow4(x)    ((x)*(x)*(x)*(x))\n"];
  pr["#define pow2inv(x) (1.0/((x)*(x)))\n\n"];
  pr["#define Tanh(x)    tanh(x)\n"];
  pr["#define Sech(x)    (1/cosh(x))\n"];
  pr["#define Cal(x,y,z) ((x)?(y):(z))\n\n\n"];

  pr["#define ALPHAMIN 1e-6\n"]; 
  pr["\n\n\n"];
];

functionname = "compute_hydrovars";

functionarguments = "tVarList *u";

(* with OMP: *)
(* *)
functionbegin     = "bampi_openmp_start\n"; 
functionloop      = "forinnerpoints_ijk_openmp(level) {";
endfunctionloop   = "} endfor_ijk_openmp; /* loop i, j, k */";
endfunction       = "bampi_openmp_stop\n"; 
(* *)

(* without OMP: *)
(*
 functionloop      = "forinnerpoints_ijk(level) {";
 endfunctionloop   = "} endfor_ijk; /* loop i, j, k */";
*)



(* custom variable declaration
   we have to translate between the tensor names and the C variables  
   we may need special variables, for example for parameters
*)
variabledeclarations[] := Module[{},

  prdecvl[{conD, pres, rho, epsl, v[a],
           g[a,b], alpha, beta[a], 
           xp, yp, zp,             
           Dh, Db, Du,
           etot, uesc,
           vor[a], ivor[a],
	   P[a], Pu[a], ut}, 
          "u"];
                                        
  pr["tL *level = u->level;\n"];
  
  prdecvarname[{mask}, "matter_mask"];
  pr["\n"];
                    
  pr["const double Wmax = Getd(\"hydroanalysis_wMAX\");\n"];                    
                    
                    
];


(************************************************************************)
(* now we are ready to go *)


<< "../../math/MathToC/TensorEquationsToC.m"

