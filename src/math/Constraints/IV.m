(* This script uses confrules.m to derive the conformal version of the constraint equations *) 

<<confrules.m

DefineTensor[Kphys,{ {2,1},1 }]
DefineTensor[Aphys,{ {2,1},1 }]
DefineTensor[A,{ {2,1},1 }]
trK

(* Hamiltonian and Momentum Constraints *)
H;
DefineTensor[P,{ {1} }];
DefineTensor[P2,{ {1} }];

(* TT decomposiontion variables *)
DefineTensor[V,{ {1} }];
DefineTensor[M,{ {2,1},1 }]
DefineTensor[LVphys,{ {2,1},1 }]
DefineTensor[LV,{ {2,1},1 }]

DefineTensor[Ptemp,{ {2,1},1 }];

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

(* Now, the conformal, traverse traceless decomposition *)

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
(* traverse traceless decomposition *)

RuleUnique[TT,A[a_,b_], Metricg[a,uc] CD[V[b],lc] + Metricg[b,uc] CD[V[a],lc] - 2/3 Metricg[a,b] CD[V[uc],lc]+ M[a,b],UpperIndexQ[a]&&UpperIndexQ[b]];

P[ua] = ApplyRules[P[ua],TT];

P[ua] = P[ua] /. CD[CD[V[u1],l2],l1] ->  - RiemannR[u1,l3,l2,l1] V[u3] + CD[CD[V[u1],l1],l2]

(* And that's it. We specify trK, M^{ij} and g_{ij}, and solve for V^{i} and psi *)

(* Now the physical transverse traceless decomposition *)

(* careful, CD is with respect to physical metric *)

P2[ua] = CD[Ptemp[ua,ub],lb]

RuleUnique[gphysmetric,CD[gphys[a_,b_],c_],0];
P2[ua] = Absorb[ApplyRules[P2[ua], gphysmetric],gphys];
RuleUnique[TTphys,Aphys[a_,b_], LVphys[a,b] + Mphys[a,b],UpperIndexQ[a]&&UpperIndexQ[b]];
P2[ua] = ApplyRules[P2[ua], TTphys];

RuleUnique[Qtransv,CD[Q[a_,b_],c_],0];
P2[ua] = ApplyRules[P2[ua],Qtransv]
P2[ua] = CDtoOD[CDtoOD[P2[ua]]] /. {AffineG[a_,b_,c_] -> Gphys[a,b,c], Metricg[a_,b_] -> gphys[a,b]}

RuleUnique[LVConfRule,LVphys[a_,b_], psi^(-4) LV[a,b],UpperIndexQ[a]&&UpperIndexQ[b]];
RuleUnique[MConfRule1,Mphys[a_,b_], psi^(-2) M[a,b] , LowerIndexQ[a]&& LowerIndexQ[b]];

RuleUnique[MConfRule2,Mphys[a_,b_], psi^(-10) M[a,b] ,UpperIndexQ[a]&& UpperIndexQ[b]];

MConfRule := {MConfRule1,MConfRule2};

P2[ua] = ApplyRules[P2[ua],LVConfRule];
P2[ua] = ApplyRules[P2[ua],MConfRule];
P2[ua] = ApplyRules[P2[ua],AffineConfRule];
P2[ua] = ApplyRules[P2[ua],MetricConfRule]

RuleUnique[ODtoCD9,OD[LV[a_,b_],c_],CD[LV[a,b],c] - AffineG[a,ld,c] LV[ud,b] - AffineG[b,ld,c] LV[ud,a], UpperIndexQ[a]&&UpperIndexQ[b]&&LowerIndexQ[c]]
RuleUnique[ODtoCD10,OD[trK,a_],CD[trK,a],LowerIndexQ[a]];
RuleUnique[ODtoCD11,OD[M[a_,b_],c_],CD[M[a,b],c] - AffineG[a,ld,c] M[ud,b] - AffineG[b,ld,c] M[ud,a], UpperIndexQ[a]&&UpperIndexQ[b]&&LowerIndexQ[c]]

P2[ua] = ApplyRules[P2[ua],ODtoCD9];
P2[ua] = ApplyRules[P2[ua],ODtoCD10];
P2[ua] = ApplyRules[P2[ua],ODtoCD11];
P2[ua] = ApplyRules[P2[ua],ODtoCD1];
P2[ua] = Absorbg[P2[ua]];

RuleUnique[LVtracefree,LV[a_,b_],0,PairQ[a,b]]
RuleUnique[Mtracefree,M[a_,b_],0,PairQ[a,b]]

P2[ua] = ApplyRules[P2[ua],LVtracefree];
P2[ua] = ApplyRules[P2[ua],Mtracefree];

(* multiply by psi^4 to get same form as Cook *)

P2[ua] = Expand[psi^4 P2[ua]];

RuleUnique[LVinsertdef,LV[a_,b_],Metricg[a,uc] CD[V[b],lc] + Metricg[b,uc] CD[V[a],lc] - 2/3 Metricg[a,b] CD[V[uc],lc], UpperIndexQ[a]&&UpperIndexQ[b]];

P2[ua] = ApplyRules[P2[ua], LVinsertdef]

P2[ua] = Expand[P2[ua] /. CD[CD[V[u1],l2],l1] ->  - RiemannR[u1,l3,l2,l1] V[u3] + CD[CD[V[u1],l1],l2]]

P2[ua] = Collect[P2[ua],CD[psi,l1]]

(*and there we have it *)

confTTeq = {gphys[la,lb] == ApplyRules[gphys[la,lb],MetricConfRule], Kphys[ua,ub] == Absorbg[ApplyRules[ApplyRules[ApplyRules[Kphys[ua,ub],KphysToAphys],MetricConfRule],ExtrinsicConfRule]],A[ua,ub] == ApplyRules[A[ua,ub],TT], ham == H, mom == Expand[P[ua]]}

physTTeq = {gphys[la,lb] == ApplyRules[gphys[la,lb],MetricConfRule], Kphys[ua,ub] == ApplyRules[Absorb[ApplyRules[Kphys[ua,ub],KphysToAphys],gphys],MetricConfRule],Aphys[ua,ub] == Collect[ApplyRules[ApplyRules[ApplyRules[ApplyRules[Aphys[ua,ub],TTphys],LVConfRule],MConfRule],LVinsertdef],psi],ham==H,mom==Expand[P2[ua]]}
