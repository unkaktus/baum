(* for physical TT transformation, prove equation 35 in Cooks paper *)

<<confrules.m

DefineTensor[LW,{ {2,1},1 }]

(*confirm rule for transforming LW *)

LW[ua,ub] = CD[W[ub],ua] + CD[W[ua],ub] - 2/3 gphys[ua,ub] CD[W[uc],lc]

LW[ua,ub] = CDtoOD[LW[ua,ub]]/. {AffineG[a_,b_,c_] -> Gphys[a,b,c],Metricg[a_,b_] -> gphys[a,b]}

LW[ua,ub] = ApplyRules[LW[ua,ub],AffineConfRule];
LW[ua,ub] = ApplyRules[LW[ua,ub],MetricConfRule];

RuleUnique[afCD1,AffineG[a_,b_,c_] Metricg[e_,f_] W[g_], CD[W[a],c] Metricg[e,f] - OD[W[a],c] Metricg[e,f],PairQ[b,g]];
RuleUnique[afCD2,AffineG[a_,b_,c_] Metricg[e_,f_] W[g_], CD[W[a],b] Metricg[e,f] - OD[W[a],b] Metricg[e,f],PairQ[c,g]];

LW[ua,ub] = ApplyRules[LW[ua,ub],afCD1];
LW[ua,ub] = ApplyRules[LW[ua,ub],afCD2];

LW[ua,ub] = Absorbg[LW[ua,ub]]

(* This is the correct result *)
