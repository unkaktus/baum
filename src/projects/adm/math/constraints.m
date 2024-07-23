<<init.m
Dimension=3
Off[MetricgFlag];
On[SyntaxCheck]

Off[General::spell, General::spell1, PermWeight::sym, PermWeight::def];

<<rules.m

DefineTensor[g,{ {2,1},1 }]
DefineTensor[gu,{ {2,1},1 }]
DefineTensor[dg,{ {2,1,3},1 }]
DefineTensor[dgu,{ {2,1,3},1 }]
DefineTensor[ddg,{ {2,1,4,3},1 }]

RicS = ScalarRtoAffine[ScalarR]
RicS = AffineToMetric[RicS]
RicS = RicS /. Metricg[u1, u2]*OD[Metricg[l3, l4], l2]*OD[Metricg[u3, u4], l1] -> Metricg[u2, u1]*OD[Metricg[l3, l4], l1]*OD[Metricg[u3, u4], l2]

RicS = ApplyRules[RicS,adMetricgRicS]

RicS = CanNonInvert[RicS /. Metricg[u1, u2]*Metricg[u3, u4]*Metricg[u5, u6]*OD[Metricg[l1, l4], l6]*OD[Metricg[l2, l3], l5] -> Metricg[u5, u6]*OD[Metricg[u2, u3], l6]*OD[Metricg[l2, l3], l5]];
RicS = CanNonInvert[RicS /. Metricg[u1, u2]*Metricg[u3, u4]*Metricg[u5, u6]*OD[Metricg[l1, l6], l4]*OD[Metricg[l2, l3], l5] -> Metricg[u3, u4]*OD[Metricg[u2, u5], l4]*OD[Metricg[l2, l3], l5]];
RicS = CanNonInvert[RicS /. Metricg[u1, u2]*Metricg[u3, u4]*Metricg[u5, u6]*OD[Metricg[l1, l4], l6]*OD[Metricg[l2, l5], l3] -> Metricg[u5, u6]*OD[Metricg[u2, u3], l6]*OD[Metricg[l2, l5], l3]];

RicS = CanNonInvert[RicS /. Metricg[u1, u2]*Metricg[u3, u4]*Metricg[u5, u6]*OD[Metricg[l1, l6], l4]*OD[Metricg[l2, l5], l3] -> Metricg[u3, u4]*OD[Metricg[u2, u5], l4]*OD[Metricg[l2, l5], l3]];

RicS = CanNonInvert[RicS /. Metricg[u1, u2]*Metricg[u3, u4]*OD[Metricg[l3, l4], l1, l2] -> Metricg[u1, u2]*Metricg[u3, u4]*OD[Metricg[l1, l2], l3, l4]]

RicS = CanNonInvert[RicS /. Metricg[u1, u2]*Metricg[u3, u4]*OD[Metricg[l1, l4], l2, l3] -> Metricg[u1, u2]*Metricg[u3, u4]*OD[Metricg[l2, l3], l1, l4]]

RuleUnique[ MetricgTOg, Metricg[a_,b_], g[a,b],LowerIndexQ[a]&&LowerIndexQ[b]]
RuleUnique[ MetricgTOgu, Metricg[a_,b_], gu[a,b],UpperIndexQ[a]&&UpperIndexQ[b]]
(* Rules to replace all derivs *)
RuleUnique[ ODgTOdg, OD[Metricg[a_,b_],c_], dg[a,b,c],LowerIndexQ[a]&&LowerIndexQ[b]&&LowerIndexQ[c]]
RuleUnique[ ODguTOdgu, OD[Metricg[a_,b_], c_], dgu[a,b,c],UpperIndexQ[a]&&UpperIndexQ[b]&&LowerIndexQ[c]]

RuleUnique[ ODODgTOddg, OD[Metricg[a_,b_],c_,d_], ddg[a,b,c,d],LowerIndexQ[a]&&LowerIndexQ[b]&&LowerIndexQ[c]&&LowerIndexQ[d]]

namerules := {ODODgTOddg,ODgTOdg,ODguTOdgu,MetricgTOg,MetricgTOgu}


DummyToBAMDummy = { u1->p, u2->q, u3->r, u4->s, u5->t, u6->u, 
			   u7->v, u8->w, u9->x, u10->y, u11->z,
		    l1->p, l2->q, l3->r, l4->s, l5->t, l6->u,
			   l7->v, l8->w, l9->x, l10->y, l11->z   }
IndexToBAMIndex = { ui->i, uj->j, uk->k,
		    li->i, lj->j, lk->k }

RicS = ApplyRules[ApplyRules[RicS,namerules],DummyToBAMDummy]

Req = {{RicSdef == RicS}}

Req >> Requation.m
