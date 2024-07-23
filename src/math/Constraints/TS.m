(* This script uses confrules.m to derive the conformal version of the constraint equations *) 

<<confrules.m

DefineTensor[uphys,{ {2,1},1 }]
DefineTensor[uphystemp,{ {2,1},1 }]
DefineTensor[u,{ {2,1},1 }]
DefineTensor[dotgphys,{ {2,1},1 }]
DefineTensor[dotgphystemp,{ {2,1},1 }]
DefineTensor[dotg,{ {2,1},1 }]
DefineTensor[beta,{ {1} }];
DefineTensor[LBphys,{ {2,1},1 }]
DefineTensor[LB,{ {2,1},1 }]
DefineTensor[Kphys,{ {2,1},1 }]
DefineTensor[Aphys,{ {2,1},1 }]
DefineTensor[A,{ {2,1},1 }]
DefineTensor[Atemp,{ {2,1},1 }]

trK
(* Hamiltonian and Momentum Constraints *)
H;
DefineTensor[P,{ {1} }];
DefineTensor[Ptemp,{ {2,1},1 }];

(* First we proove the conformal transmation of the metric *)

detgphys = -(gphys[-1,-3]*gphys[-1,-3]*gphys[-2,-2]) + 2*gphys[-1,-2]*gphys[-1,-3]*gphys[-2,-3] -gphys[-1,-1]*gphys[-2,-3]*gphys[-2,-3] - gphys[-1,-2]gphys[-1,-2]*gphys[-3,-3] + gphys[-1,-1]*gphys[-2,-2]*gphys[-3,-3]

(* transforming each component gives: *)

detg = (psi^4)^3 (-(Metricg[-1,-3]*Metricg[-1,-3]*Metricg[-2,-2]) + 2*Metricg[-1,-2]*Metricg[-1,-3]*Metricg[-2,-3] -Metricg[-1,-1]*Metricg[-2,-3]*Metricg[-2,-3] - Metricg[-1,-2]Metricg[-1,-2]*Metricg[-3,-3] + Metricg[-1,-1]*Metricg[-2,-2]*Metricg[-3,-3]);

(*which gives the rule that*)

RuleUnique[DetgConfRule,Detgphys, psi^(12) Detg];
(* NOTE: I use OD[ T, lo] to denote time deriv of T *)

RuleUnique[ ODgloTodotg, OD[Metricg[la_,lb_],lo], dotg[la,lb] ];

uphystemp[la,lb] = Detgphys^(1/3) OD[Detgphys^(-1/3) gphys[la,lb],lo];

uphystemp[la,lb] = ApplyRules[ApplyRules[uphystemp[la,lb],MetricConfRule],DetgConfRule];
RuleUnique[detg0,OD[Detg,lo],0];
uphystemp[la,lb] = ApplyRules[uphystemp[la,lb],detg0];

RuleUnique[dotgu, OD[Metricg[a_,b_],lo],u[a,b],LowerIndexQ[a]&&LowerIndexQ[b]];

uphystemp[la,lb] = ApplyRules[uphystemp[la,lb],dotgu];

(* which gives the transformation rule for u *)
RuleUnique[uConfRule1,uphys[a_,b_],psi^4 u[a,b],LowerIndexQ[a]&&LowerIndexQ[b]];
RuleUnique[uConfRule2,uphys[a_,b_],psi^(-4) u[a,b],UpperIndexQ[a]&&UpperIndexQ[b]];

uConfRule := {uConfRule1,uConfRule2};

(* Proove that u_ij is the tracefree part of dt g_ij *)

RuleUnique[smallp,gphys[a_,b_] OD[gphys[c_,d_],lo],OD[Detgphys,lo]/Detgphys,PairQ[a,c]&&PairQ[b,d]] 

uphystemp[la,lb] = Detgphys^(1/3) OD[Detgphys^(-1/3) gphys[la,lb],lo];

dotgphystemp[la,lb] = uphystemp[la,lb] + 1/3 gphys[la,lb] gphys[uc,ud] OD[gphys[lc,ld],lo];

dotgphystemp[la,lb] = ApplyRules[dotgphystemp[la,lb],smallp];

(* which checks out okay *)

(* now the equation for uphys *)
(* split the extrinsic curvature tensor *)

RuleUnique[KphysToAphys1,Kphys[a_,b_],Aphys[a,b] + 1/3 gphys[a,b] trK, LowerIndexQ[a]&& LowerIndexQ[b]]

RuleUnique[KphysToAphys2,Kphys[a_,b_],gphys[a,uc] gphys[b,ud] (Aphys[lc,ld] + 1/3 gphys[lc,ld] trK), UpperIndexQ[a]&&UpperIndexQ[b]]
KphysToAphys := {KphysToAphys1,KphysToAphys2};

RuleUnique[trAphys0,Aphys[a_, b_],0,PairQ[a,b]];

dotgphystemp[ua_,ub_] = -2 alpha Kphys[ua,ub] + CD[beta[ua],ub] + CD[beta[ub],ua];

(* subtract the trace part *)

uphystemp[ua,ub] = dotgphystemp[ua,ub] - 1/3 gphys[ua,ub] gphys[lc,ld] dotgphystemp[uc,ud] ;

uphystemp[ua,ub] = Absorb[ApplyRules[uphystemp[ua,ub],KphysToAphys],gphys]

uphystemp[ua,ub] = ApplyRules[uphystemp[ua,ub],trAphys0];

RuleUnique[defLB,CD[beta[b_],a_] + CD[beta[a_],b_] - 2/3 gphys[a_,b_] CD[beta[c_],d_],LBphys[a,b],UpperIndexQ[a]&&UpperIndexQ[b]&&PairQ[c,d]]

uphystemp[ua,ub] = ApplyRules[uphystemp[ua,ub],defLB];

(* which is euqation 46 in cooks paper *)

RuleUnique[LBConfRule,LBphys[a_,b_], psi^(-4) LB[a,b],UpperIndexQ[a]&&UpperIndexQ[b]];

Atemp[ua,ub] = psi^10 (LBphys[ua,ub] - uphys[ua,ub])/(2 alphaphys);
Atemp[ua,ub] = AtemppplyRules[AtemppplyRules[Atemp[ua,ub],LBConfRule],uConfRule];

RuleUnique[alphaConfRule,alphaphys,psi^6 alpha]

Atemp[ua,ub] = AppplyRules[Atemp[ua,ub],alphaConfRule];

RuleUnique[TSA, A[a_,b_], LB[a, b]/(2*alpha) - u[a, b]/(2*alpha),UpperIndexQ[a]&&UpperIndexQ[ub]];

(* which is expression 49 in Cooks paper *)

(* Define the hamiltonian constraint, note that we multiply by -1/8 * psi^5to get the same form as in cooks paper, NB: the extrinsic curvature is _not_ on conformal form *)

H = Expand[-1/8 * psi^5*(ApplyRules[scaRphys,RicciSConfRule] + (trK)^2 - Kphys[la,lb]*Kphys[ua,ub])]

(* split the extrinsic curvature tensor *)

RuleUnique[KphysToAphys1,Kphys[a_,b_],Aphys[a,b] + 1/3 gphys[a,b] trK, LowerIndexQ[a]&& LowerIndexQ[b]]

RuleUnique[KphysToAphys2,Kphys[a_,b_],gphys[a,uc] gphys[b,ud] (Aphys[lc,ld] + 1/3 gphys[lc,ld] trK), UpperIndexQ[a]&&UpperIndexQ[b]]
KphysToAphys := {KphysToAphys1,KphysToAphys2};

(* rule that Aphys is trace free. NB: this is a hack, it might not always work *)
RuleUnique[trAphys0,Aphys[a_, b_],0,PairQ[a,b]];

H = ApplyRules[H,KphysToAphys]

(* absorb physical metric*)

H = Absorb[H,gphys]

H = ApplyRules[H,trAphys0]

(* Define the momentum constraint, NB: the extrinsic curvature is _not_ on conformal form *)

Ptemp[ua_,ub_] = Kphys[ua,ub] - gphys[ua,ub] trK;
Ptemp[ua_,ub_] = ApplyRules[Ptemp[ua,ub],KphysToAphys];
Ptemp[ua_,ub_] = Absorb[Ptemp[ua,ub],gphys];

P[ua] = OD[Ptemp[ua,ub],lb] + Gphys[ua,ld,lb] Ptemp[ud,ub] + Gphys[ub,ld,lb] Ptemp[ua,ud];

P[ua] = Dum[Expand[P[ua]]];

(* this is about as good as it's going to get without doing the transformation of A_{ij} *)

(* Now, confromal transformation decomposition *)

RuleUnique[ExtrinsicConfRule1,Aphys[a_,b_], psi^(-2) A[a,b] , LowerIndexQ[a]&& LowerIndexQ[b]];

RuleUnique[ExtrinsicConfRule2,Aphys[a_,b_], psi^(-10) A[a,b] ,
UpperIndexQ[a]&& UpperIndexQ[b]];

ExtrinsicConfRule := {ExtrinsicConfRule1,ExtrinsicConfRule2};

(* hamiltonian constraint *)

H = ApplyRules[H,ExtrinsicConfRule];

(* momentum constraint *)

P[ua] = ApplyRules[P[ua],ExtrinsicConfRule];

P[ua] = ApplyRules[P[ua],MetricConfRule];
P[ua] = ApplyRules[P[ua],AffineConfRule];

RuleUnique[ODtoCD6,OD[A[a_,b_],c_],CD[A[a,b],c] - AffineG[a,ld,c] A[ud,b] - AffineG[b,ld,c] A[a,ud],UpperIndexQ[a]&&UpperIndexQ[b]&&LowerIndexQ[c]];
RuleUnique[ODtoCD7,OD[Metricg[a_,b_],c_],CD[Metricg[a,b],c] - AffineG[a,ld,c] Metricg[ud,b] - AffineG[b,ld,c] Metricg[a,ud],UpperIndexQ[a]&&UpperIndexQ[b]&&LowerIndexQ[c]];

RuleUnique[ODtoCD8,OD[Metricg[a_,b_],c_],CD[Metricg[a,b],c] + AffineG[ud,a,c] Metricg[ld,b] + AffineG[ud,b,c] Metricg[ld,a],LowerIndexQ[a]&&LowerIndexQ[b]&&LowerIndexQ[c]];

P[ua] = ApplyRules[P[ua],ODtoCD6];
P[ua] = ApplyRules[P[ua],ODtoCD7];
P[ua] = ApplyRules[P[ua],ODtoCD8];

P[ua] = Tsimplify[Absorbg[P[ua]]];
RuleUnique[trA0,A[a_, b_],0,PairQ[a,b]];

P[ua] = ApplyRules[P[ua],trA0];

(* multiply with psi^10 to get same form as in Cooks paper *)

P[ua] = Expand[psi^(10) P[ua]];

(* now we can input eqution 49 *)

P[ua] = ApplyRules[P[ua],TSA];

(* multiply with 2 alpha to get same for as Cook eq 50 *)

P[ua] = Expand[2 alpha P[ua]];

(* checking, that Cooks term is the same as mine *) 
alpha CD[(1/alpha) u[ua,ub],lb];

(* get equation for dotgphys *)

dotgphystemp[la,lb] = uphys[la,lb] + 1/3 gphys[la,lb] gphys[uc,ud] OD[gphys[lc,ld],lo]

RuleUnique[InsertEV,OD[gphys[a_,b_],lo],- 2 alphaphys Kphys[a,b] + CD[beta[b],a] + CD[beta[a],b],LowerIndexQ[a]&&LowerIndexQ[b]];

dotgphystemp[la,lb] = Absorb[ApplyRules[dotgphystemp[la,lb],InsertEV],gphys];

RuleUnique[trKphys,Kphys[a_, b_],trK,PairQ[a,b]];

dotgphystemp[la,lb] = ApplyRules[dotgphystemp[la,lb],trKphys]

RuleUnique[updown,CD[beta[a_],b_],CD[beta[b],a],PairQ[a,b]];

RuleUnique[CD2ODmixed,CD[beta[a_],b_],OD[beta[a],b] + Gphys[a,b,lc] beta[uc],UpperIndexQ[a]&&LowerIndexQ[b]];

dotgphystemp[la,lb] = ApplyRules[ApplyRules[dotgphystemp[la,lb],updown],CD2ODmixed]

dotgphystemp[la,lb] = ApplyRules[dotgphystemp[la,lb],AffineConfRule]
dotgphystemp[la,lb] = ApplyRules[dotgphystemp[la,lb],MetricConfRule];
dotgphystemp[la,lb] = ApplyRules[dotgphystemp[la,lb],uConfRule];
dotgphystemp[la,lb] = ApplyRules[dotgphystemp[la,lb],alphaConfRule];
dotgphystemp[la,lb] = ApplyRules[dotgphystemp[la,lb],ODtoCD1];

RuleUnique[ODtoCD9,OD[beta[a_], b_], CD[beta[a],b] - AffineG[a,b,lc] beta[uc],UpperIndexQ[a]&&LowerIndexQ[b]];

dotgphystemp[la,lb] = ApplyRules[dotgphystemp[la,lb],ODtoCD9];

dotgphystemp[la,lb] = Absorbg[dotgphystemp[la,lb]];

(* which is Cooks equation 52 *)

(* All in all we have *)


TSeq = {gphys[la,lb] == ApplyRules[gphys[la,lb],MetricConfRule], Kphys[ua,ub] == Absorbg[ApplyRules[ApplyRules[ApplyRules[Kphys[ua,ub],KphysToAphys],MetricConfRule],ExtrinsicConfRule]],A[ua,ub] == ApplyRules[A[ua,ub],TSA], ham == H, mom == Expand[P[ua]],dotgphys[la,lb] == dotgphystemp[la,lb]}



