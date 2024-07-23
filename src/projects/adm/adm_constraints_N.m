(* admconstraints.m 
   Bernd Bruegmann 11/02 *)

(* compute ADM constraints from ADM variables *)


(* variables *)
variables = {g[a,b], K[a,b], psi, dpop[a], ddpop[a,b], ham, mom[a],
             normham, normmom[a],rho,S[a],
             xp, yp, zp, rp,rr }



(* compute in this order *)
tocompute = {

  (* partial derivatives *)
  Cinstruction == "if (order_centered == 2 || boundary1away) {",
    delg[c,a,b]       == del[c,g[a,b]],
    deldelg[a,b,c,d]  == deldel[a,b,g[c,d]],
    delK[a,b,c]       == del[a, K[b,c]],
  Cinstruction == "} else if (order_centered == 4 || boundaryNaway(2)) {",
    delg[c,a,b]       == del4[c,g[a,b]],
    deldelg[a,b,c,d]  == deldel4[a,b,g[c,d]],
    delK[a,b,c]       == del4[a, K[b,c]],
  Cinstruction == "} else if (order_centered == 6 || boundaryNaway(3)){",
    delg[c,a,b]       == del6[c,g[a,b]],
    deldelg[a,b,c,d]  == deldel6[a,b,g[c,d]],
    delK[a,b,c]       == del6[a, K[b,c]],
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
    delgSST[a,b,c]      == Jac[d,a] delg[d,b,c],
    deldelgSST[a,b,c,d] == DJac[e,a,b] delg[e,c,d] + Jac[e,a] Jac[f,b] deldelg[e,f,c,d],
    delKSST[a,b,c]      == Jac[d,a] delK[d,b,c],
    
    (* now give the tmp derivative the same name again *)
    delg[a,b,c]         == delgSST[a,b,c],
    deldelg[a,b,c,d]    == deldelgSST[a,b,c,d],
    delK[a,b,c]         == delKSST[a,b,c],
    
  Cinstruction == "}",
  
  
  (* make transition to physical metric *)
  f == 1,
  Cif == usepsi, 
    f == psi^4,
    delf[a] == 4 f dpop[a],
    deldelf[a,b] == 4 f (ddpop[a,b] + 3 dpop[a] dpop[b]), 
    deldelg[a,b,c,d] == f deldelg[a,b,c,d] + deldelf[a,b] g[c,d] +
                        delf[a] delg[b,c,d] + delf[b] delg[a,c,d],  
    delg[c,a,b] == f delg[c,a,b] + g[a,b] delf[c],
  Cif == end,

  (* inverse physical metric *)
  detginvf == 1/(f matrixdet[g]),
  ginv[a,b] == detginvf matrixinvdet[g,a,b],

  (* connection of physical metric *) 
  gammado[c,a,b] == 1/2 (delg[a,b,c] + delg[b,a,c] - delg[c,a,b]),
  gamma[c,a,b] == ginv[c,d] gammado[d,a,b], 

  (* curvature of physical metric *)
  R[a,b] == ginv[c,d] ( 1/2 (
    - deldelg[c,d,a,b] - deldelg[a,b,c,d] + 
      deldelg[a,c,b,d] + deldelg[b,c,a,d]) +
      gamma[e,a,c] gammado[e,b,d] - 
      gamma[e,a,b] gammado[e,c,d]),
  R == ginv[a,b] R[a,b],

  (* extrinsic curvature terms *)
  Kud[a,b] == ginv[a,c] K[c,b],
  K == Kud[a,a],
  KudKud == Kud[a,b] Kud[b,a],
  codelK[a,b,c] == delK[a,b,c] - gamma[d,a,b] K[d,c] - gamma[d,a,c] K[b,d],
  cdKudd[a,b,c] == ginv[a,d] codelK[d,b,c],

  (* Hamiltonian constraint *)
  Cif == matter,
        ham == R + K K - KudKud - 16. PI rho,
  Cif == else,
        ham == R + K K - KudKud,
  Cif == end,

  (* momentum constraint *)
  momu[a] == ginv[a,b] cdKudd[c,b,c] - ginv[b,c] cdKudd[a,b,c],
  Cif == matter,
        mom[a] == f g[a,b] ( momu[b]  - 8. PI S[b] ),
  Cif == else,
        mom[a] == f g[a,b] momu[b],
  Cif == end,


  trK == K,

  Cif == normConstr,

    (* normalized Hamiltonian constraint *)
    (* normham == ham/(fabs[R] + K K + fabs[Kud[a,b] Kud[b,a]]), *)
    Cif == TermByTerm,
      (* compute fabs of some terms in R separately *)
      RA[a,b] == ginv[c,d] * (-deldelg[c,d,a,b]),
      RB[a,b] == ginv[c,d] * ( deldelg[a,c,b,d]),
      RC[a,b] == ginv[c,d] * ( gamma[e,a,c] gammado[e,b,d]),
      RD[a,b] == ginv[c,d] * (-gamma[e,a,b] gammado[e,c,d]),
      RA == ginv[a,b] RA[a,b],
      RB == ginv[a,b] RB[a,b],
      RC == ginv[a,b] RC[a,b],
      RD == ginv[a,b] RD[a,b],
      (* hamnum == RA + RB + RC + RD + K K - KudKud, *)
      denom  == fabs[RA]+fabs[RB]+fabs[RC]+fabs[RD] + K K + fabs[KudKud],
    Cif == else,
      denom  == fabs[R] + K K + fabs[KudKud],
    Cif == end,

    normham == Cal[denom <= 0.0, 0.0, ham/denom],

    (* normalized momentum constraint *)
    (* normmom[a] == mom[a]/( f*( fabs[ cdKudd[c,a,c] ] 
                                 +fabs[ ginv[b,c] codelK[a,b,c] ] ) ), *)
    Cif == TermByTerm,
      (* compute fabs of some terms in codelK separately *)
      codelKA[a,b,c] ==  delK[a,b,c],
      codelKB[a,b,c] == -gamma[d,a,b] K[d,c],
      codelKC[a,b,c] == -gamma[d,a,c] K[b,d],
      cdKdA[b] == ginv[a,c] codelKA[a,b,c],
      cdKdB[b] == ginv[a,c] codelKB[a,b,c],
      cdKdC[b] == ginv[a,c] codelKC[a,b,c],
      codelTrKA[a] == ginv[b,c] codelKA[a,b,c],
      codelTrKB[a] == ginv[b,c] codelKB[a,b,c],
      codelTrKC[a] == ginv[b,c] codelKC[a,b,c],
      (* note: mom[la] =  CD[ K[la,ue], le ] - CD[ TrK, la ] *)
      (* momnum[a] == cdKdA[a] + cdKdB[a] + cdKdC[a] -
                      codelTrKA[a] - codelTrKB[a] - codelTrKC[a], *)
      denom[a] ==  fabs[cdKdA[a]] + fabs[cdKdB[a]] + fabs[cdKdC[a]] +
                   fabs[codelTrKA[a]] + fabs[codelTrKB[a]] +
                   fabs[codelTrKC[a]],
    Cif == else,
      denom[a]  == fabs[cdKudd[c,a,c]] + fabs[ginv[b,c] codelK[a,b,c]],
    Cif == end,

    normmom[a] == Cal[denom[a] <= 0.0, 0.0, mom[a]/denom[a]],
 
  Cif == end
}


(* symmetries *)
g[a_,b_]       := g[b,a] /; !OrderedQ[{a,b}]
K[a_,b_]       := K[b,a] /; !OrderedQ[{a,b}]
R[a_,b_]       := R[b,a] /; !OrderedQ[{a,b}]
ginv[a_,b_]    := ginv[b,a] /; !OrderedQ[{a,b}]
ddpop[a_,b_]   := ddpop[b,a] /; !OrderedQ[{a,b}]
deldelf[a_,b_] := deldelf[b,a] /; !OrderedQ[{a,b}]

deldel[a_,b_,c_] := deldel[b,a,c] /; !OrderedQ[{a,b}]
deldel4[a_,b_,c_] := deldel4[b,a,c] /; !OrderedQ[{a,b}]

delg[c_,a_,b_] := delg[c,b,a] /; !OrderedQ[{a,b}]
delK[c_,a_,b_] := delK[c,b,a] /; !OrderedQ[{a,b}]
deldelg[a_,b_,c_,d_] := deldelg[b,a,c,d] /; !OrderedQ[{a,b}]
deldelg[a_,b_,c_,d_] := deldelg[a,b,d,c] /; !OrderedQ[{c,d}]

codelK[c_,a_,b_] := codelK[c,b,a] /; !OrderedQ[{a,b}]
cdKudd[c_,a_,b_] := cdKudd[c,b,a] /; !OrderedQ[{a,b}]

RA[a_,b_]         := RA[b,a] /; !OrderedQ[{a,b}]
RB[a_,b_]         := RB[b,a] /; !OrderedQ[{a,b}]
RC[a_,b_]         := RC[b,a] /; !OrderedQ[{a,b}]
RD[a_,b_]         := RD[b,a] /; !OrderedQ[{a,b}]
codelKA[c_,a_,b_] := codelKA[c,b,a] /; !OrderedQ[{a,b}]
codelKB[c_,a_,b_] := codelKB[c,b,a] /; !OrderedQ[{a,b}]
codelKC[c_,a_,b_] := codelKC[c,b,a] /; !OrderedQ[{a,b}]
cdKuddA[c_,a_,b_] := cdKuddA[c,b,a] /; !OrderedQ[{a,b}]
cdKuddB[c_,a_,b_] := cdKuddB[c,b,a] /; !OrderedQ[{a,b}]
cdKuddC[c_,a_,b_] := cdKuddC[c,b,a] /; !OrderedQ[{a,b}]

delgSST[c_,a_,b_]       := delgSST[c,b,a] /; !OrderedQ[{a,b}]
deldelgSST[a_,b_,c_,d_] := deldelgSST[b,a,c,d] /; !OrderedQ[{a,b}]
deldelgSST[a_,b_,c_,d_] := deldelgSST[a,b,d,c] /; !OrderedQ[{c,d}]
delKbSST[c_,a_,b_]      := delKbSST[c,b,a] /; !OrderedQ[{a,b}]

DJac[a_,b_,c_]       := DJac[a,c,b] /; !OrderedQ[{b,c}]


(************************************************************************)
(* information for C output *)

(* information about the function 
   the C file will be functionname.c 
*)
includeANDdefine[] := Module[{},
  pr["#include \"bam.h\"\n"];
  pr["#include \"adm.h\"\n"];
  pr["\n"];
  pr["#define Power(x,y) (pow((double) (x), (double) (y)))\n"];
  pr["#define Sqrt(x)    sqrt(x)\n"];
  pr["#define Log(x)     log((double) (x))\n"];
  pr["#define pow2(x)    ((x)*(x))\n"];
  pr["#define pow4(x)    ((x)*(x)*(x)*(x))\n"];
  pr["#define pow2inv(x) (1.0/((x)*(x)))\n\n"];
  pr["#define Tanh(x)    tanh(x)\n"];
  pr["#define Sech(x)    (1/cosh(x))\n"];
  pr["#define Cal(x,y,z) ((x)?(y):(z))\n\n"];
  pr["\n\n\n"];

];

functionname      = "adm_constraints_N";

functionarguments = "tVarList *ucon, tVarList *uadm";

functionbegin     = "bampi_openmp_start\n";

functionloop      = "forinnerpoints_ijk_openmp(level) {";

endfunctionloop   = "} endfor_ijk_openmp; /* loop i, j, k */";

endfunction       = "bampi_openmp_stop\n";




(* custom variable declaration
   we have to translate between the tensor names and the C variables  
   we may need special variables, for example for parameters
*)
variabledeclarations[] := Module[{},

  pr["tL *level = ucon->level;\n"];
  pr["\n"];
  prdecvl[{ham,mom[a], normham,normmom[a]}, "ucon"];
  prdecvl[{g[a,b], K[a,b], psi, dpop[a], ddpop[a,b], rho, S[a]}, "uadm"];
  pr["\n"];

  pr["const int normConstr = Getv(\"adm_normalizedConstraints\", \"yes\");\n"];
  pr["const int TermByTerm = Getv(\"adm_normalizedConstraints\", \"TermByTerm\");\n"];
  pr["const int usepsi = 0;\n"];
  pr["const int order_centered = Geti(\"order_centered\");\n"];
  
  pr["const int matter = ( (Getv(\"physics\",\"matter\")) && (!level->shells) );\n"];
  
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

