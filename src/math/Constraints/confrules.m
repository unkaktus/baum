(* This MathTensor script demonstrates how to derive the rules for the
conformal transformations from scratch *)

<<MathTensor.m
Dimension=3
Off[MetricgFlag];
On[SyntaxCheck]
On[Syntax::com];

Off[General::spell, General::spell1, PermWeight::sym, PermWeight::def];

psi (*conformal factor*);

DefineTensor[gphys,{ {2,1},1}]    
DefineTensor[Gphys,{{1,3,2},1}];
DefineTensor[Gtemp,{{1,3,2},1}];

(* 3D Riemann *)

permweights = Symmetries[RiemannR[la,lb,lc,ld]];

DefineTensor[Rietemp,permweights]
DefineTensor[Riephys,permweights] 

(* 3D Ricci *)

DefineTensor[Ricphys,{ {2,1}, 1 }]
DefineTensor[Rictemp,{ {2,1}, 1 }]
scaRphys
scaRtemp

Off[MetricgFlag]; (* never absorb the metric unless I use Absorbg *)
(* Trace of metric is 3 *)
DefUnique[Metricg[la_,ua_], Dimension, PairQ[la,ua] ]
DefUnique[gphys[ua_,lb_], Dimension, PairQ[ua,lb] ]
 (* Rules for the conformal transformation of the metric tensor*)

RuleUnique[MetricConfRule1,gphys[a_,b_], psi^4 Metricg[a,b] , LowerIndexQ[a]&& LowerIndexQ[b]];

RuleUnique[MetricConfRule2,gphys[a_,b_], psi^(-4) Metricg[a,b] ,
UpperIndexQ[a]&& UpperIndexQ[b]];

MetricConfRule := {MetricConfRule1,MetricConfRule2};

 (* Rule for affine connection *) 

Gtemp[ua, lb,lc] = Dum[1/2*gphys[ua,ud]*(OD[gphys[lb,ld],lc] +
OD[gphys[lc,ld],lb] - OD[gphys[lb,lc],ld])]; 
Gtemp[ua,lb,lc] = ApplyRules[Gtemp[ua,lb,lc],MetricConfRule];
Gtemp[ua,lb,lc] = Gtemp[ua,lb,lc] - AffineToMetric[AffineG[ua,lb,lc]] +
AffineG[ua,lb,lc]

(* insert output of above into rule: *)

RuleUnique[AffineConfRule,Gphys[a_,b_,c_], AffineG[a, b, c] -
 (2*Metricg[b, c]*Metricg[u1, a]*OD[psi, l1])/psi + (2*Metricg[l1,
 c]*Metricg[u1, a]*OD[psi, b])/psi + (2*Metricg[l1, b]*Metricg[u1,
 a]*OD[psi, c])/psi, UpperIndexQ[a]&&LowerIndexQ[b]&&
 LowerIndexQ[c]];

(*check to see if it works, remove ; to see output*)
ApplyRules[Gphys[ua,lb,lc],AffineConfRule];

(* Now: the Riemann tensor *)

Rietemp[ua,lb,lc,ld] = Dum[Gphys[u1, lb, ld]*Gphys[ua, l1, lc] -
  Gphys[u1, lb, lc]* Gphys[ua, l1, ld] - OD[Gphys[ua, lb, lc], ld] +
  OD[Gphys[ua, lb, ld], lc]];

Rietemp[ua,lb,lc,ld] =
ApplyRules[Rietemp[ua,lb,lc,ld],AffineConfRule];
Rietemp[ua,lb,lc,ld] = Rietemp[ua,lb,lc,ld] -
RiemannToAffine[RiemannR[ua,lb,lc,ld]] + RiemannR[ua,lb,lc,ld];

RuleUnique[ODtoCD1,OD[psi,a_],CD[psi,a],LowerIndexQ[a]];

RuleUnique[ODtoCD2,OD[psi,b_,c_],OD[CD[psi,b],c],LowerIndexQ[b]&&LowerIndexQ[c]];

RuleUnique[ODtoCD3,OD[CD[psi,b_],c_],CD[CD[psi,b],c] +AffineG[ua,b,c] CD[psi,la],LowerIndexQ[b]&&LowerIndexQ[c]];

RuleUnique[ODtoCD4,OD[Metricg[a_, b_], c_],CD[Metricg[a,b],c] + AffineG[ue,a,c]Metricg[le,b] + AffineG[ue,b,c] Metricg[a,le],LowerIndexQ[a]&&LowerIndexQ[b]&&LowerIndexQ[c]];

RuleUnique[ODtoCD5,OD[Metricg[a_, b_], c_],CD[Metricg[a,b],c] -
AffineG[a,le,c] Metricg[ue,b] - AffineG[b,le,c] Metricg[ue,a],UpperIndexQ[a]&&UpperIndexQ[b]&&LowerIndexQ[c]];

Rietemp[ua,lb,lc,ld] = ApplyRules[Rietemp[ua,lb,lc,ld],ODtoCD1];
Rietemp[ua,lb,lc,ld] = ApplyRules[Rietemp[ua,lb,lc,ld],ODtoCD2];
Rietemp[ua,lb,lc,ld] = ApplyRules[Rietemp[ua,lb,lc,ld],ODtoCD3];
Rietemp[ua,lb,lc,ld] = ApplyRules[Rietemp[ua,lb,lc,ld],ODtoCD4];

Rietemp[ua,lb,lc,ld] =Absorbg[ Rietemp[ua,lb,lc,ld]];
Rietemp[ua,lb,lc,ld]=Collect[Absorbg[ApplyRules[Rietemp[ua,lb,lc,ld],ODtoCD5]],psi]; 
(* insert output of above into rule: *)

(*Here's the rule we have to adjust for the fact that mathtensor
insists on reordering the indecees before applying the rule*) 

RuleUnique[RiemannConfRule,Riephys[b_,a_,c_,d_], -((2*CD[psi, d,a]*Metricg[b, c] - 2*CD[psi, c, a]*Metricg[b, d] - 2*CD[psi,b, d]*Metricg[c, a] + 2*CD[psi, b, c]*Metricg[d, a])/ psi +(-6*CD[psi, d]*CD[psi, a]*Metricg[b, c] + 6*CD[psi, c]*CD[psi,a]*Metricg[b, d] + 6*CD[psi, b]*CD[psi, d]* Metricg[c, a] -4*CD[psi, le]*CD[psi, ue]*Metricg[b, d]* Metricg[c, a] -6*CD[psi, b]*CD[psi, c]*Metricg[d, a] + 4*CD[psi, le]*CD[psi,ue]*Metricg[b, c]*Metricg[d, a])/psi^2 + RiemannR[a, b, c,d]), UpperIndexQ[a]&&LowerIndexQ[b]&&LowerIndexQ[c]&&LowerIndexQ[d]];

(*check to see if it works, remove ; to see output*)

ApplyRules[Riephys[ua,lb,lc,ld],RiemannConfRule];

Rictemp[lb,ld] =
Dum[Absorbg[Metricg[la,uc]*ApplyRules[Riephys[ua,lb,lc,ld],RiemannConfRule]]];

(* insert output of above into rule: *)

RuleUnique[RicciConfRule,Ricphys[a_,b_],(6*CD[psi, a]*CD[psi, b])/psi^2
 - (2*CD[psi, a, b])/psi - (2*CD[psi, l1]*CD[psi, u1]*Metricg[a,
 b])/psi^2 - (2*CD[psi, l1, u1]*Metricg[a, b])/psi + RicciR[a, b],LowerIndexQ[a]&&LowerIndexQ[b]]

(*check to see if it works, remove ; to see output*)

ApplyRules[Ricphys[la,lb],RicciConfRule];

(* and finally the ricci scalar *)

scaRtemp =Collect[ApplyRules[gphys[ua,ub],MetricConfRule]*ApplyRules[Ricphys[la,lb],RicciConfRule],psi];

scaRtemp = Collect[Dum[scaRtemp],psi];

RuleUnique[RicciSConfRule,scaRphys,ScalarR/psi^4 - (8*CD[psi, l1, u1])/psi^5];



