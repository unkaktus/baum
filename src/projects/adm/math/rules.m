
RuleUnique[aMetricg1, Metricg[ua_,uc_] Metricg[lc_,lb_] , Metricg[ua,lb], PairQ[uc,lc] ] 
RuleUnique[aMetricg2, Metricg[ua_,uc_] Metricg[lb_,lc_] , Metricg[ua,lb], PairQ[uc,lc] ] 
RuleUnique[aMetricg3, Metricg[uc_,ua_] Metricg[lc_,lb_] , Metricg[ua,lb], PairQ[uc,lc] ]
RuleUnique[aMetricg4, Metricg[uc_,ua_] Metricg[lb_,lc_] , Metricg[ua,lb], PairQ[uc,lc] ]

aMetricg := {aMetricg1,aMetricg2,aMetricg3,aMetricg4}

RuleUnique[adMetricg1, Metricg[ua_,ub_]*Metricg[uc_,ud_]*Metricg[ue_,uf_]*OD[Metricg[la_,le_],lf_]*OD[Metricg[lb_,ld_],lc_],-Metricg[ue,uf]*OD[Metricg[ua,uc],lc]*OD[Metricg[la,le],lf],PairQ[ua,la]&&PairQ[ub,lb]&&PairQ[uc,lc]&&PairQ[ud,ld]&&PairQ[ue,le]&&PairQ[uf,lf]]

RuleUnique[adMetricg2, Metricg[ua_,ub_]*Metricg[uc_,ud_]*Metricg[ue_,uf_]*OD[Metricg[la_,lf_],le_]*OD[Metricg[lb_,ld_],lc_],-Metricg[ue,uf]*OD[Metricg[ua,uc],lc]*OD[Metricg[la,lf],le],PairQ[ua,la]&&PairQ[ub,lb]&&PairQ[uc,lc]&&PairQ[ud,ld]&&PairQ[ue,le]&&PairQ[uf,lf]]

RuleUnique[adMetricg3, Metricg[ua_,ub_]*Metricg[uc_,ud_]*Metricg[ue_,uf_]*OD[Metricg[la_,le_],lf_]*OD[Metricg[lb_,lc_],ld_],-Metricg[ue,uf]*OD[Metricg[la,le],lf]*OD[Metricg[ua,ud],ld],PairQ[ua,la]&&PairQ[ub,lb]&&PairQ[uc,lc]&&PairQ[ud,ld]&&PairQ[ue,le]&&PairQ[uf,lf]]

RuleUnique[adMetricg4, Metricg[ua_,ub_]*Metricg[uc_,ud_]*Metricg[ue_,uf_]*OD[Metricg[la_,lf_],le_]*OD[Metricg[lb_,lc_],ld_],-Metricg[ue,uf]*OD[Metricg[la,lf],le]*OD[Metricg[ua,ud],ld],PairQ[ua,la]&&PairQ[ub,lb]&&PairQ[uc,lc]&&PairQ[ud,ld]&&PairQ[ue,le]&&PairQ[uf,lf]]

RuleUnique[adMetricg5, Metricg[ua_,ub_]*Metricg[uc_,ud_]*Metricg[ue_,uf_]*OD[Metricg[lb_,lc_],ld_]*OD[Metricg[le_,lf_],la_],-Metricg[ue,uf]*OD[Metricg[ua,ud],ld]*OD[Metricg[le,lf],la],PairQ[ua,la]&&PairQ[ub,lb]&&PairQ[uc,lc]&&PairQ[ud,ld]&&PairQ[ue,le]&&PairQ[uf,lf]]

RuleUnique[adMetricg6, Metricg[ua_,ub_]*Metricg[uc_,ud_]*Metricg[ue_,uf_]*OD[Metricg[lb_,ld_],lc_]*OD[Metricg[le_,lf_],la_],-Metricg[ue,uf]*OD[Metricg[ua,uc],lc]*OD[Metricg[le,lf],la],PairQ[ua,la]&&PairQ[ub,lb]&&PairQ[uc,lc]&&PairQ[ud,ld]&&PairQ[ue,le]&&PairQ[uf,lf]]

RuleUnique[adMetricg7, Metricg[ua_,ub_]*Metricg[uc_,ud_]*Metricg[ue_,uf_]*OD[Metricg[la_,ld_],lf_]*OD[Metricg[lc_,le_],lb_],-Metricg[ue,uf]*OD[Metricg[ub,uc],lf]*OD[Metricg[lc,le],lb],PairQ[ua,la]&&PairQ[ub,lb]&&PairQ[uc,lc]&&PairQ[ud,ld]&&PairQ[ue,le]&&PairQ[uf,lf]]

RuleUnique[adMetricg8, Metricg[ua_,ub_]*Metricg[uc_,ud_]*Metricg[ue_,uf_]*OD[Metricg[la_,lf_],ld_]*OD[Metricg[lc_,le_],lb_],-Metricg[uc,ud]*OD[Metricg[ub,ue],ld]*OD[Metricg[lc,le],lb],PairQ[ua,la]&&PairQ[ub,lb]&&PairQ[uc,lc]&&PairQ[ud,ld]&&PairQ[ue,le]&&PairQ[uf,lf]]

RuleUnique[adMetricg9, Metricg[ua_,ub_]*Metricg[uc_,ud_]*Metricg[ue_,uf_]*OD[Metricg[lb_,lc_],le_]*OD[Metricg[ld_,lf_],la_],-Metricg[ue,uf]*OD[Metricg[ua,ud],le]*OD[Metricg[ld,lf],la],PairQ[ua,la]&&PairQ[ub,lb]&&PairQ[uc,lc]&&PairQ[ud,ld]&&PairQ[ue,le]&&PairQ[uf,lf]]

RuleUnique[adMetricg10, Metricg[ua_,ub_]*Metricg[uc_,ud_]*Metricg[ue_,uf_]*OD[Metricg[lb_,le_],lc_]*OD[Metricg[ld_,lf_],la_],-Metricg[uc,ud]*OD[Metricg[ua,uf],lc]*OD[Metricg[ld,lf],la],PairQ[ua,la]&&PairQ[ub,lb]&&PairQ[uc,lc]&&PairQ[ud,ld]&&PairQ[ue,le]&&PairQ[uf,lf]]

RuleUnique[adMetricg11, Metricg[ua_,ub_]*Metricg[uc_,ud_]*Metricg[ue_,uf_]*OD[Metricg[lc_,le_],lb_]*OD[Metricg[ld_,lf_],la_],-Metricg[ua,ub]*OD[Metricg[ud,uf],lb]*OD[Metricg[ld,lf],la],PairQ[ua,la]&&PairQ[ub,lb]&&PairQ[uc,lc]&&PairQ[ud,ld]&&PairQ[ue,le]&&PairQ[uf,lf]]

RuleUnique[adMetricg12, Metricg[ua_,ub_]*Metricg[uc_,ud_]*Metricg[ue_,uf_]*OD[Metricg[la_,le_],lf_]*OD[Metricg[lc_,ld_],lb_],-Metricg[uc,ud]*OD[Metricg[ub,uf],lf]*OD[Metricg[lc,ld],lb],PairQ[ua,la]&&PairQ[ub,lb]&&PairQ[uc,lc]&&PairQ[ud,ld]&&PairQ[ue,le]&&PairQ[uf,lf]]

RuleUnique[adMetricg13, Metricg[ua_,ub_]*Metricg[uc_,ud_]*Metricg[ue_,uf_]*OD[Metricg[la_,lf_],le_]*OD[Metricg[lc_,ld_],lb_],-Metricg[uc,ud]*OD[Metricg[ub,ue],le]*OD[Metricg[lc,ld],lb],PairQ[ua,la]&&PairQ[ub,lb]&&PairQ[uc,lc]&&PairQ[ud,ld]&&PairQ[ue,le]&&PairQ[uf,lf]]

adMetricgRicS := {adMetricg1,adMetricg2,adMetricg3,adMetricg4,adMetricg5,adMetricg6,adMetricg7,adMetricg8,adMetricg9,adMetricg10,adMetricg11,adMetricg12,adMetricg13}

