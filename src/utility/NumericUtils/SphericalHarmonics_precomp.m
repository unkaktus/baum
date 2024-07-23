(* mth 02/2011 *)

(* Mathematica script, to write Ynlm


*)  

functionname = "SphericalHarmonics_precomp";
lmax = 6;









(***************************************************************************)
(* compute Re Ylm and Im Ylm *)

Ylm[l_, m_] := Module[{}, 
    tmp = Simplify[SphericalHarmonicY[l,m,theta,phi]];
    valR = FullSimplify[ComplexExpand[Re[ tmp ]]];
    valI = FullSimplify[ComplexExpand[Im[ tmp ]]];
    
    pr["    *rY "];
    pr[" = "];
    pr[ CForm[valR] ];
    pr[";\n"];
    pr["    *iY "];
    pr[" = "];
    pr[ CForm[valI] ];
    pr[";\n"];
    pr["    return;\n"];
]





(***************************************************************************)
(* compute -2 spin weighted Re Ylm and Im Ylm *)
(*
    "Gauge - invariant Non - spherical Metric Perturbations of Schwarzschild Black - Hole Spacetimes"
    Authors : Alessandro Nagar, Luciano Rezzolla
    Class.Quant.Grav.22 (2005) R167; Erratum - ibid.23 (2006) 4297
    http : // arxiv.org/abs/gr - qc/0502064
    
    see mathematica script in src/utility/NumericUtils/math 
*)


(* Spherical Harmonics and drvts *)
Y[l_, m_, \[Theta]_, \[Phi]_] := 
  SphericalHarmonicY[l, m, \[Theta], \[Phi]] ;
dYd\[Theta][l_, m_, \[Theta]_, \[Phi]_] := 
  D[SphericalHarmonicY[l, m, \[Theta], \[Phi]] , \[Theta]];
d2Yd\[Theta]2[l_, m_, \[Theta]_, \[Phi]_] := 
  D[SphericalHarmonicY[l, m, \[Theta], \[Phi]] , {\[Theta], 2}];
dYd\[Phi][l_, m_, \[Theta]_, \[Phi]_] := 
  D[SphericalHarmonicY[l, m, \[Theta], \[Phi]] , \[Phi]];
d2Yd\[Phi]2[l_, m_, \[Theta]_, \[Phi]_] := 
  D[SphericalHarmonicY[l, m, \[Theta], \[Phi]] , {\[Phi], 2 }];
d2Yd\[Theta]\[Phi][l_, m_, \[Theta]_, \[Phi]_] := 
  D[D[SphericalHarmonicY[l, m, \[Theta], \[Phi]] , \[Theta]] , \[Phi]];


(* Angular functions W and X [Eq.(39-40)] *)
W[l_, m_, \[Theta]_, \[Phi]_] := 
  d2Yd\[Theta]2[l, m, \[Theta], \[Phi]] - 
   Cot[\[Theta]] dYd\[Theta][l, m, \[Theta], \[Phi]] - 
   d2Yd\[Phi]2[l, m, \[Theta], \[Phi]]/(Sin[\[Theta]]^2);
X[l_, m_, \[Theta]_, \[Phi]_] := 
  2 ( d2Yd\[Theta]\[Phi][l, m, \[Theta], \[Phi]] - 
     Cot[\[Theta]] dYd\[Phi][l, m, \[Theta], \[Phi]] );


(* Spin - weighted s = -2 spherical harmonics [Eq.(66)] *)
Ys[l_, m_, \[Theta]_, \[Phi]_] := 
  Sqrt[(l - 2)!/(l + 2)!] ( 
    W[l, m, \[Theta], \[Phi]]  - 
     I X[l, m, \[Theta], \[Phi]] /Sin[\[Theta]]);
YsNN[l_, m_, \[Theta]_, \[Phi]_] :=  
  W[l, m, \[Theta], \[Phi]]  - 
   I X[l, m, \[Theta], \[Phi]] /Sin[\[Theta]];


(* print function *)
m2Ylm[l_, m_] := Module[{}, 
    tmp = Simplify[Ys[l,m,theta,phi]];
    valR = FullSimplify[ComplexExpand[Re[ tmp ]]];
    valI = FullSimplify[ComplexExpand[Im[ tmp ]]];
    
    pr["    *rY "];
    pr[" = "];
    pr[ CForm[valR] ];
    pr[";\n"];
    pr["    *iY "];
    pr[" = "];
    pr[ CForm[valI] ];
    pr[";\n"];
    pr["    return;\n"];
]











(***************************************************************************)
(* loop over all needed l's and m's *)

collectall[lmax_,func_] :=  Module[{l,m}, 
    For[l=2, l<=lmax, l++,
        pr[ StringForm["if (l==``) {\n\n", l] ];
        For[m=-l, m<=l, m++,
            If[m==-l,pr["  "],];
            pr[ StringForm["if (m==``) {\n", m] ];
                func[l,m];
            If[m==l,
                pr["  }\n\n"], 
                pr["  } else "]
            ];
        ];
        If[l==lmax,
            pr["} else errorexit(\"increase lmax in mathematica script\");\n\n"], 
            pr["} else "]
        ];
    ];
]












(***************************************************************************)
(* write C code *)

Off[Part::partd]
Off[Part::pspec]

file = functionname <> ".c";
Print["Writing to ", file, "\n"];
DeleteFile[file];

filepointer = OpenAppend[file];
pr[x_] := Module[{s},
  WriteString[filepointer, x];
]


pr["/* "<>file<>" */\n"];
pr["/* mth, "<>
       ToString[Date[][[3]]]<>"."<>
       ToString[Date[][[2]]]<>"."<>
       ToString[Date[][[1]]]<>" */\n"];
pr["/* Produced with Mathematica */\n"];
pr["\n"];
pr["#include \"bam.h\"\n"];
pr["\n"];
pr["#define Power(x,y) (pow((double) (x), (double) (y)))\n"];
pr["#define Log(x)     log((double) (x)))\n"];
pr["#define pow2(x)    ((x)*(x))\n"];
pr["#define pow2inv(x) (1.0/((x)*(x)))\n"];
pr["#define Cal(x,y,z) ((x)?(y):(z))\n"];
pr["#define Pi         PI\n"];
pr["#define Cos(x)     cos(x)\n"];
pr["#define Sin(x)     sin(x)\n"];
pr["#define Sqrt(x)    sqrt(x)\n"];
pr["#define Csc(x)     (1./sin(x))\n"];
pr["\n\n\n"];




pr["\n\n\n\n"];
pr["void SphericalHarmonicYprecomp(double *rY, double *iY, int l, int m, double phi, double theta)\n{\n"];
pr["\n"];
collectall[lmax, Ylm];
pr["\n"];
pr["}\n"];


pr["\n\n\n\n"];
pr["void spinweightedSphericalHarmonicprecomp(double *rY, double *iY, int l, int m, double phi, double theta)\n{\n"];
pr["\n"];
collectall[lmax, m2Ylm];
pr["\n"];
pr["}\n"];



pr["/* "<>file<>" */\n"];


Close[filepointer];

prtimer[]


