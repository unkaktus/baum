(* andw 04/11
   mth  05/11 
   Compute Jacobians for transformation from sphericalshells to cartesian coordinates 
   do it for all 6 patches, opposit patches have the same transformation rule
   the transformation is from src/main/amr/shells.c
*)

JacobianMatrix[f_List?VectorQ, x_List] := Outer[D, f, x]; 
DJacobianMatrix[f_List?VectorQ, x_List] := Outer[D, Outer[D, f, x], x];
DDJacobianMatrix[f_List?VectorQ, x_List] := Outer[D, Outer[D, Outer[D, f, x], x], x];





(* coordinate transfo for fisheye coordinates *)
(* Pollney/Reisswig 
fisheyef3[r_] := a*(r-R)+b*Sqrt[1+(r-R)^2/e];
fisheyef2[r_] := fisheyef3[r] /.b->(s-a)*Sqrt[e];
fisheyef1[r_] := fisheyef2[r] /.a->(e+R(R+Sqrt[e+R^2]S))/(e+R(R+Sqrt[e+R^2]));
fisheyef[r_]  := fisheyef1[r] /.{s->shellsS,R->shellsR,e->shellsE};
*)

(* mth *)
fisheyef3[r_] := b*r + a*e*Log[Cosh[(r-R)/e]];
fisheyef2[r_] := fisheyef3[r] /.b->(s-a);
fisheyef1[r_] := fisheyef2[r] /.a->1/2*(1+Exp[-2*R/e])* (s-1);
fisheyef[r_]  := fisheyef1[r] /.{s->shellsS,R->shellsR,e->shellsE};

(* now the derivatives needed inside the Jacobians *)
fisheyeR[r_] := fisheyef[r] - fisheyef[0]; 
DRDr  = Simplify[ D[fisheyeR[rr],rr] ];
DDRDr = Simplify[ D[fisheyeR[rr],{rr,2}] ];









(* for RHS  transform from box-spherical coordinates to cartesian coordinates *)

(* -+x *)
rp[xp_, yp_, zp_] = F[Sqrt[xp^2 + yp^2 + zp^2]];
pp[xp_, yp_, zp_] = ArcTan[yp/xp];
tp[xp_, yp_, zp_] = ArcTan[zp/xp];

JacobB01   = JacobianMatrix[{rp[xp, yp, zp], pp[xp, yp, zp], tp[xp, yp, zp]}, {xp, yp, zp}];
JacobB01   = Simplify[JacobB01/. xp^2 + yp^2 + zp^2 -> rp^2, rp>0];
JacobB01   = Simplify[JacobB01/. Derivative[1][F][rp] -> 1/dRdr];
JacobB01   = Simplify[JacobB01];

DJacobB01  = DJacobianMatrix[{rp[xp, yp, zp], pp[xp, yp, zp], tp[xp, yp, zp]}, {xp, yp, zp}];
DJacobB01  = Simplify[DJacobB01/. xp^2 + yp^2 + zp^2 -> rp^2, rp>0];
DJacobB01  = Simplify[DJacobB01/. Derivative[1][F][rp] -> 1/dRdr];
DJacobB01  = Simplify[DJacobB01/. Derivative[2][F][rp] -> -ddRdr/dRdr^3];
DJacobB01  = Simplify[DJacobB01];


(* -+y *)
rp[xp_, yp_, zp_] = F[Sqrt[xp^2 + yp^2 + zp^2]];
pp[xp_, yp_, zp_] = ArcTan[xp/yp];
tp[xp_, yp_, zp_] = ArcTan[zp/yp];

JacobB23   = JacobianMatrix[{rp[xp, yp, zp], pp[xp, yp, zp], tp[xp, yp, zp]}, {xp, yp, zp}];
JacobB23   = Simplify[JacobB23/. xp^2 + yp^2 + zp^2 -> rp^2, rp>0];
JacobB23   = Simplify[JacobB23/. Derivative[1][F][rp] -> 1/dRdr];
JacobB23   = Simplify[JacobB23];

DJacobB23  = DJacobianMatrix[{rp[xp, yp, zp], pp[xp, yp, zp], tp[xp, yp, zp]}, {xp, yp, zp}];
DJacobB23  = Simplify[DJacobB23/. xp^2 + yp^2 + zp^2 -> rp^2, rp>0];
DJacobB23  = Simplify[DJacobB23/. Derivative[1][F][rp] -> 1/dRdr];
DJacobB23  = Simplify[DJacobB23/. Derivative[2][F][rp] -> -ddRdr/dRdr^3];
DJacobB23  = Simplify[DJacobB23];


(* -+z *)
rp[xp_, yp_, zp_] = F[Sqrt[xp^2 + yp^2 + zp^2]];
pp[xp_, yp_, zp_] = ArcTan[xp/zp];
tp[xp_, yp_, zp_] = ArcTan[yp/zp];

JacobB45   = JacobianMatrix[{rp[xp, yp, zp], pp[xp, yp, zp], tp[xp, yp, zp]}, {xp, yp, zp}];
JacobB45   = Simplify[JacobB45/. xp^2 + yp^2 + zp^2 -> rp^2, rp>0];
JacobB45   = Simplify[JacobB45/. Derivative[1][F][rp] -> 1/dRdr];
JacobB45   = Simplify[JacobB45];

DJacobB45  = DJacobianMatrix[{rp[xp, yp, zp], pp[xp, yp, zp], tp[xp, yp, zp]}, {xp, yp, zp}];
DJacobB45  = Simplify[DJacobB45/. xp^2 + yp^2 + zp^2 -> rp^2, rp>0];
DJacobB45  = Simplify[DJacobB45/. Derivative[1][F][rp] -> 1/dRdr];
DJacobB45  = Simplify[DJacobB45/. Derivative[2][F][rp] -> -ddRdr/dRdr^3];
DJacobB45  = Simplify[DJacobB45];



















(* for boundary  transform from cartesian coordinates to spherical coordinates *)
rp[xp_, yp_, zp_] = Sqrt[xp^2 + yp^2 + zp^2];
pp[xp_, yp_, zp_] = ArcTan[yp/xp];
tp[xp_, yp_, zp_] = ArcCos[zp/ Sqrt[xp^2 + yp^2 + zp^2]];

Jacob   = JacobianMatrix[{rp[xp, yp, zp], pp[xp, yp, zp], tp[xp, yp, zp]}, {xp, yp, zp}];
Jacob   = FullSimplify[Jacob, xp^2 + yp^2 + zp^2 > 0];

DJacob  = DJacobianMatrix[{rp[xp, yp, zp], pp[xp, yp, zp], tp[xp, yp, zp]}, {xp, yp, zp}];
DDJacob = DDJacobianMatrix[{rp[xp, yp, zp], pp[xp, yp, zp], tp[xp, yp, zp]}, {xp, yp, zp}];

xp[rp_, pp_, tp_] = rp Sin[tp]Cos[pp];
yp[rp_, pp_, tp_] = rp Sin[tp]Sin[pp];
zp[rp_, pp_, tp_] = rp Cos[tp];

rule      = {rp -> rp[xp, yp, zp], pp -> pp[xp, yp, zp], tp -> tp[xp, yp, zp]};
JacobI   = JacobianMatrix[{xp[rp, pp, tp], yp[rp, pp, tp], zp[rp, pp, tp]}, {rp, pp, tp}];
JacobI   = FullSimplify[JacobI /. rule, {xp^2 + yp^2 + zp^2 > 0, xp > 0}];
detJI    = Simplify[Det[JacobI]];

DJacobI   = DJacobianMatrix[{xp[rp, pp, tp], yp[rp, pp, tp], zp[rp, pp, tp]}, {rp, pp, tp}];
DJacobI   = FullSimplify[DJacobI /. rule, {xp^2 + yp^2 + zp^2 > 0, xp > 0}];







(*
rule      = {rp -> rp[xp, yp, zp], pp -> pp[xp, yp, zp] - \[Pi], tp -> tp[xp, yp, zp]};
JacobImp  = JacobianMatrix[{xp[rp, tp, pp], yp[rp, tp, pp], zp[rp, tp, pp]}, {rp, tp, pp}];
JacobImp  = FullSimplify[JacobImp /. rule, {xp^2 + yp^2 + zp^2 > 0, xp < 0}];
detJImp   = Simplify[Det[JacobImp]];


rule      = {rp -> rp[xp, yp, zp], pp -> pp[xp, yp, zp] + \[Pi], tp -> tp[xp, yp, zp]};
JacobImm  = JacobianMatrix[{xp[rp, tp, pp], yp[rp, tp, pp], zp[rp, tp, pp]}, {rp, tp, pp}];
JacobImm  = FullSimplify[JacobImm /. rule, {xp^2 + yp^2 + zp^2 > 0, xp < 0}];
detJImm   = Simplify[Det[JacobImm]];
*)






(*
DJacobIp  = DJacobianMatrix[{xp[rp, tp, pp], yp[rp, tp, pp], zp[rp, tp, pp]}, {rp, tp, pp}];
DJacobIm  = DJacobianMatrix[{xp[rp, tp, pp], yp[rp, tp, pp], zp[rp, tp, pp]}, {rp, tp, pp}];
DDJacobIp = DDJacobianMatrix[{xp[rp, tp, pp], yp[rp, tp, pp], zp[rp, tp, pp]}, {rp, tp, pp}];
DDJacobIm = DDJacobianMatrix[{xp[rp, tp, pp], yp[rp, tp, pp], zp[rp, tp, pp]}, {rp, tp, pp}];
*)









timer["after jacobians.m"]
prtimer[]






