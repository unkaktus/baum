(* writeC.m *)
(* Bernd Bruegmann, 2/96, 10/02 *)

(* Mathematica script, to be read by TensorEquationsToC.m:

   Given equations, write C code:
   - write to file functionname.c
   - write header of predefined C statements
   - start function definition with predefined argument list
   - write definitions of C variables
   - write equations as C assignments
   - end loop and close function
*)  


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
pr["/* Copyright (C) 1998 Bernd Bruegmann, "<>
       ToString[Date[][[3]]]<>"."<>
       ToString[Date[][[2]]]<>"."<>
       ToString[Date[][[1]]]<>" */\n"];
pr["/* Produced with Mathematica */\n"];
pr["\n"];


includeANDdefine[];

pr["void " <> functionname <> "(" <> functionarguments <> ")\n{\n"];

pr["\n"];



(***************************************************************************)
(* write definition of C variables *)

prdecvar[vars_, name_] := Module[{cvars, s, nf}, 
  cvars = vars /. expfreeindices /. gluevar;
  For[nf = 0, nf < Length[cvars], nf++,
    s = StringForm["double *`` = level->v[``+``];\n", 
                    cvars[[nf+1]], name, nf];
    pr[s];
  ];
]

prdecvarname[vars_, name_] := Module[{cvars, s, nf}, 
  cvars = vars /. expfreeindices /. gluevar;
  s = StringForm["int index_`` = Ind(\"``\");\n", name, name];
  pr[s];
  For[nf = 0, nf < Length[cvars], nf++,
      (* level should not be hardwired *)
      s = StringForm["double *`` = level->v[index_`` + ``];\n", 
                    cvars[[nf+1]], name, nf];
    pr[s];
  ];
]

prdecvl[vars_, name_] := Module[{cvars, s, nf}, 
  cvars = vars /. expfreeindices /. gluevar;
  For[nf = 0, nf < Length[cvars], nf++,
    s = StringForm["double *`` = vldataptr(``, ``);\n", 
                    cvars[[nf+1]], name, nf];
    pr[s];
  ];
]

prdecvlindices[vars_, name_] := Module[{cvars, s, nf},
  cvars = vars /. expfreeindices /. gluevar;
  For[nf = 0, nf < Length[cvars], nf++,
    s = StringForm["int index_`` = (``)->index[``];\n",
                    cvars[[nf+1]], name, nf];
    pr[s];
  ];
]

prdecvlNR[vars_, name_, nr_] := Module[{cvars, s, nf}, 
  cvars = vars /. expfreeindices /. gluevar;
  For[nf = 0, nf < Length[cvars], nf++,
    s = StringForm["double *`` = vldataptr(``, ``);\n", 
                    cvars[[nf+1]], name, nf+nr];
    pr[s];
  ];
]

prdecvlf[vars_, name_] := Module[{cvars, s, nf}, 
  s = StringForm["tVarList *vl_`` = ``(level);\n", name, name];
  pr[s];
  cvars = vars /. expfreeindices /. gluevar;
  For[nf = 0, nf < Length[cvars], nf++,
    s = StringForm["double *`` = vldataptr(vl_``, ``);\n", 
                    cvars[[nf+1]], name, nf];
    pr[s];
  ];
  s = StringForm["vlfree(vl_``);\n", name];
  pr[s];
]

writeCdef[vars_, auxvars_] := Module[{cvars, s, nf, nfields, nl},

  variabledeclarations[];
  pr["\n"];
  
  If[ValueQ[functionbegin]==True,pr["\n\n" <> functionbegin <> "\n\n"]];

  If[TensorEquationsDim === 4, pr["double dt = level->dt;\n"]];
  pr["double dx = level->dx;\n"];
  pr["double dy = level->dy;\n"];
  pr["double dz = level->dz;\n"];
  If[TensorEquationsDim === 4,
    pr["double oo2dt = 1/(2*dt);\n"];
    pr["double oodt2 = 1/(dt*dt);\n"];
    pr["double oo4dtdx = 1/(4*dt*dx);\n"];
    pr["double oo4dtdy = 1/(4*dt*dy);\n"];
    pr["double oo4dtdz = 1/(4*dt*dz);\n"]
  ];
  pr["double oo2dx = 1/(2*dx), oo2dy = 1/(2*dy), oo2dz = 1/(2*dz);\n"];
  pr["double oodx2 = 1/(dx*dx), oody2 = 1/(dy*dy), oodz2 = 1/(dz*dz);\n"];
  pr["double oo4dxdy = 1/(4*dx*dy);\n"];
  pr["double oo4dxdz = 1/(4*dx*dz);\n"];
  pr["double oo4dydz = 1/(4*dy*dz);\n"];
  pr["\n"];
  
  nfields = Length[auxvars];
  For[nf = 1, nf <= nfields, nf++,
    s = StringForm["double `` = 0.;\n", auxvars[[nf]]];
    pr[s];
  ];
  pr["\n\n\n"]
]

writeCdef[vars, auxvars]


(***************************************************************************)
(* begin loop over all space points *)

writeCloop[evolve_] := Module[ {nf, nfields},
  pr[functionloop <> "\n\n\n"];
];

writeCloop[True];



(***************************************************************************)
(* write equations *)

tocover = Join[constvars, {exp, log, pow, pow2, pow2inv, sqrt, Cal, sin, cos,
                           fabs, fmax, fmin}]
coverpsi   = Map[ (#[i__] :> #[{i}])&, tocover]
uncoverpsi = Map[ (#[[{i__}]] :> #[i])&, tocover]

writeCassign[e_] := Module[{x, nfields, nf},

  If[e[[1]] =!= Cif && e[[1]] =!= Cinstruction && e[[1]] =!= Cloop,
    x = e /. coverpsi;
    x = ToString[InputForm[x]];
    Print["Before: ",x];
    x = StringReplace[x, "]" -> "]]"];
    x = StringReplace[x, "[" -> "[["];
    x = ToExpression[x];
    Print["After: ",x];
    x = x /. uncoverpsi /. colpow;

    x = x /. Sqrt[[x_]]  -> Sqrt[x];
    x = x /. Power[x_,2] -> pow2[x];
    x = x /. Power[x_,3] -> pow3[x];
    x = x /. Power[x_,4] -> pow4[x];
    x = x /. Tan[[x_]]   -> Tan[x];        
    x = x /. Cos[[x_]]   -> Cos[x];        
    x = x /. Sin[[x_]]   -> Sin[x];        
    x = x /. Sec[[x_]]   -> 1/Sin[x];
    
    (* N[] introduces -1. since 5.0, which we do not want *)
    (* stopped working in 5.1: -1. x_ -> -x *)
    x = x /. -1 -> minus;   
    x = N[x, 20];
    x = x /. minus -> -1;


    PutAppend[CForm[x[[1]]], file];
    pr["=\n"];        (* to avoid assignment in CForm *)
    PutAppend[CForm[x[[2]]], file];
    pr[";\n\n"];
  ];

  If[e[[1]] === Cif,
    If[e[[2]] === end,
      pr["}\n"]; 
      pr[StringForm["/* if (``) */\n\n\n", lastif]]
    ];
    If[e[[2]] === else,
      pr["\n"]; 
      pr[StringForm["} else { /* if (!``) */\n\n", lastif]]
    ];
    If[e[[2]] =!= end && e[[2]] =!= else,
      pr["\n\n"]; 
      pr["/* conditional */\n"];
      pr[StringForm["if (``) {\n\n", e[[2]]]];
      lastif = e[[2]]
    ] 
  ];

  (* Cinstruction: a piece of C-code *)
  If[e[[1]] === Cinstruction,
      pr["\n"]; 
      pr[StringForm["`` \n\n", e[[2]]]];
  ];

]

Map[writeCassign, dif];


(***************************************************************************)
(* end loop over space points *)


pr["\n\n" <> endfunctionloop <> "\n\n"];
If[ValueQ[endfunction]==True,pr["\n\n" <> endfunction <> "\n\n"]];
pr["}  /* function */\n\n"];



(* count operations *)

char = ReadList[file, Character];
nmul = Count[char, "*"];
ndiv = Count[char, "/"];
nsum = Count[char, "+"] + Count[char, "-"];
count = StringForm[
 "/* nvars = ``, nauxs = ``, n* = ``,  n/ = ``,  n+ = ``, n = ``, O = `` */\n",
  Length[vars], Length[auxvars], nmul, ndiv, nsum, nmul+ndiv+nsum, 
  If [optimizeflag, 1, 0, 0]
];
Print[""];
Print[file];
Print[count];
pr["/* "<>file<>" */\n"];
pr[count];

Close[filepointer];

timer["after writeC.m"]
prtimer[]

