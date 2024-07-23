(* expandindices.m *)
(* Bernd Bruegmann, 2/96, 10/02 Nina Jansen 01.03*)

(* Mathematica script, to be read by TensorEquationsToC.m:

   Given equations in abstract tensor notation, expand indices:
   - expand contractions to sums with explicit indices
   - expand free indices by writing several equations with explicit indices
   - glue explicit indices to variable name, g[1,2] becomes g12
 This script expands 4 dimnesional indices, NJ
*)  




(*************************************************************************)
(* expand abstract indices *)

expcontraction = {
  s_[a___,b_Symbol,c___,d_Symbol,e___] \
    t_[f___,b_Symbol,g___] u_[h___,d_Symbol,i___] ->
    Sum[s[a,b,c,d,e] t[f,b,g] u[h,d,i], {b,1,3}, {d,1,3}],

  s_[a___,b_Symbol,c___] t_[d___,b_Symbol,e___] ->
    Sum[s[a,b,c] t[d,b,e], {b,1,3}],

  s_[a___,b_Symbol,c___,b_Symbol,d___] ->
    Sum[s[a,b,c,b,d], {b,1,3}]
}


expcontraction4 = {
  s_[a___,b_Symbol,c___,d_Symbol,e___] \
    t_[f___,b_Symbol,g___] u_[h___,d_Symbol,i___] ->
    Sum[s[a,b,c,d,e] t[f,b,g] u[h,d,i], {b,0,3}, {d,0,3}],

  s_[a___,b_Symbol,c___] t_[d___,b_Symbol,e___] ->
    Sum[s[a,b,c] t[d,b,e], {b,0,3}],

  s_[a___,b_Symbol,c___,b_Symbol,d___] ->
    Sum[s[a,b,c,b,d], {b,0,3}]
}


(* careful: Union sorts, so apply to elements for which this is ok or wanted *)
expfreesort[x_,a_Symbol] := Union[Flatten[Table[x, {a,1,3}]]]
expfreesort[x_,{a_Symbol}] := expfreesort[x,a]
expfreesort[x_,{a_Symbol,b__}] := 
  Union[Flatten[Table[expfreesort[x,{b}], {a,1,3}]]]

expfreesort4[x_,a_Symbol] := Union[Flatten[Table[x, {a,0,3}]]]
expfreesort4[x_,{a_Symbol}] := expfreesort4[x,a]
expfreesort4[x_,{a_Symbol,b__}] := 
  Union[Flatten[Table[expfreesort4[x,{b}], {a,0,3}]]]


expexpsign = -x_ :> x
expfreeexp[x_Symbol] := x
expfreeexp[x_[a__]] := 
  DeleteCases[Union[Flatten[expfreesort[x[a], {a}]/.expexpsign]], 0]

expexpsign4 = -x_ :> x
expfreeexp4[x_Symbol] := x
expfreeexp4[x_[a__]] := 
  DeleteCases[Union[Flatten[expfreesort4[x[a], {a}]/.expexpsign4]], 0]


(* careful: when called on equations, LHS and RHS must have same syms *)
expequsign = -x_ == y_ :> x == -y
expfreeequ[x_Symbol==y_] := x==y
expfreeequ[delt[x_Symbol]==y_] := delt[x] == y 
expfreeequ[x_[a__]==y_] := 
  DeleteCases[Union[Flatten[expfreesort[x[a]==y,{a}]/.expequsign]], True]
expfreeequ[delt[x_[a__]]==y_] := 
  DeleteCases[Union[Flatten[expfreesort[delt[x[a]]==y,{a}]/.expequsign]], True]

expfreeequ4[x_Symbol==y_] := x==y
expfreeequ4[delt[x_Symbol]==y_] := delt[x] == y 
expfreeequ4[x_[a__]==y_] := 
  DeleteCases[Union[Flatten[expfreesort4[x[a]==y,{a}]/.expequsign]], True]
expfreeequ4[delt[x_[a__]]==y_] := 
  DeleteCases[Union[Flatten[expfreesort4[delt[x[a]]==y,{a}]/.expequsign]], True]


(* expand free indices:
   for equations, expansion is based on lhs
   for expressions, only top level list is handled
*)
expfreeindices = {
  x_ :> Flatten[Map[expfreeequ,x]] /; !FreeQ[x,y_==z_],
  x_ :> Flatten[Map[expfreeexp,x]] /; FreeQ[x,y_==z_]
}

expfreeindices4 = {
  x_ :> Flatten[Map[expfreeequ4,x]] /; !FreeQ[x,y_==z_],
  x_ :> Flatten[Map[expfreeexp4,x]] /; FreeQ[x,y_==z_]
}


(* glue indices to tensor name, g[1,2] becomes g12 *)
glue[a_,b_] := (* SequenceForm[a,b]; *) 
 ToExpression[ToString[a] <> ToString[b]]; 
glueall[x_] := x //. y_Symbol[a_,b___] :> glue[y,a][b] /. y_[] :> y
gluevariables[variables_] := Module[{varexpanded, vartoglue, varglued},
  varexpanded = Flatten[Map[expfreeexp, variables]];
  vartoglue = Cases[varexpanded, x_[a__]];
  varglued = Map[glueall, vartoglue]; 
  Table[vartoglue[[i]] -> varglued[[i]], {i, 1, Length[varglued]}]
]
gluevariables4[variables_] := Module[{varexpanded, vartoglue, varglued},
  varexpanded = Flatten[Map[expfreeexp4, variables]];
  vartoglue = Cases[varexpanded, x_[a__]];
  varglued = Map[glueall, vartoglue]; 
  Table[vartoglue[[i]] -> varglued[[i]], {i, 1, Length[varglued]}]
]


(* utilities *)
colpow = {
  Power[x_,2] -> pow2[x], 
  Power[x_,-2] -> pow2inv[x]
}

prlist[a_,b_] := Block[{},
  Print[a, "\n"];
  Map[(Print[#,"\n"])&, b];
  Print["\n"];
];




(***************************************************************************)
(* transform equations into form that is to be differenced *)


(* substitutions *)
all4 = tocompute4 /. expcodel
prlist["all4 = tocompute4 /. expcodel", all4]

(* do as often as there are contractions *)
all4 = all4 /. expcontraction4
all4 = ExpandAll[ExpandAll[all4] //. expcontraction4]
all = ExpandAll[ExpandAll[all4] //. expcontraction4]
If[False, prlist["all4 = ExpandAll[all4] /. expcontraction4 ...", all4]]

(* expand free indices *)
all4 = all4 /. expfreeindices4
If[False, prlist["all4 = all4 /. expfreeindices4", all4]]

all3 = tocompute3 /. expcodel
prlist["all3 = tocompute3 /. expcodel", all3]

(* do as often as there are contractions *)
all3 = all3 /. expcontraction
all3 = ExpandAll[ExpandAll[all3] //. expcontraction]
all = ExpandAll[ExpandAll[all3] //. expcontraction]
If[False, prlist["all3 = ExpandAll[all3] /. expcontraction3 ...", all3]]

(* expand free indices *)
all3 = all3 /. expfreeindices
If[False, prlist["all3 = all3 /. expfreeindices", all3]]

(* various sets of variables *)

(* determine auxiliary variables that require a local definition:
   - these are all4 variables that appear on the lhs of an assignment
     but are not in the list of variables
   - treat the indices as abstract by using patterns
   - note that the level specification {2} also ensures that scalars do 
     not become patterns!
*)
(* old: auxvariables = Complement[Map[#[[1]]&, tocompute4], variables]; *)

lhsvariables4 = Map[#[[1]]&, tocompute4];
lhsvariables4 = DeleteCases[lhsvariables4, Cif | Cloop | Cinstruction];
patternofvariables4 = Map[(Pattern[#, Blank[]]) &, variables, {2}];
auxvariables4 = DeleteCases[lhsvariables4, 
                           patternofvariables4 /. List -> Alternatives];
auxvariables4 = Union[auxvariables4];   (* sort and remove duplicates *)

lhsvariables3 = Map[#[[1]]&, tocompute3];
lhsvariables3 = DeleteCases[lhsvariables3, Cif | Cloop | Cinstruction];
patternofvariables3 = Map[(Pattern[#, Blank[]]) &, variables, {2}];
auxvariables3 = DeleteCases[lhsvariables3, 
                           patternofvariables3 /. List -> Alternatives];
auxvariables3 = Union[auxvariables3];   (* sort and remove duplicates *)

constvariables = {r, lr[a]}  (* obsolete *)

gluevar = Join[gluevariables[variables], 
	       gluevariables4[auxvariables4], 
	       gluevariables[auxvariables3], 
	       gluevariables[constvariables]];

vars = variables /. expfreeindices /. gluevar;
auxvars4 = auxvariables4 /. expfreeindices4 /. gluevar;
auxvars3 = auxvariables3 /. expfreeindices /. gluevar;
constvars = constvariables /. expfreeindices /. gluevar;
Print["vars"];
Print[vars];
Print["auxvars3"];
Print[auxvars3];
Print["auxvars4"];
Print[auxvars4];
Print["constvars"];
Print[constvars];
auxvars = Join[auxvars3,auxvars4];

(* introduce names for tensor components *)

all =  Join[all3,all4]
all = all /. gluevar
If[False, prlist["all /. gluevar", all]]

(* remove duplicates assuming that there are no bugs in symmetries *)
remsymdups[all_] := Module[{n,r},
  r = {}; 
  For[n = 1, n < Length[all], n++, 
    If[all[[n,1]] =!= all[[n+1,1]]   || 
       all[[n,1]] === Cif            || 
       all[[n,1]] === Cinstruction, 
       r = Append[r, all[[n]]]]];
  Append[r, all[[n]]]
]
If[True, all = remsymdups[all]]

timer["after expandindices.m"]
