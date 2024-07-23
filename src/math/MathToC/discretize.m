(* ::Package:: *)

(* discretize.m *)
(* Bernd Bruegmann, 2/96, 10/02 *)

(* Mathematica script, to be read by TensorEquationsToC.m:

   Given equations with explicit indices, discretize for numerical grid:
   - put variables onto grid, g11 becomes g11[ijk] 
   - replace derivatives by finite differences
*)  


(***************************************************************************)
(* finite differencing *)

(* centered stencils, see FiniteDifferences.nb *)
d1stencil2nd = {{-1, 0, 1}, {-1, 0, 1}, 1/(2h)};
d2stencil2nd = {{-1, 0, 1}, { 1,-2, 1}, 1/h^2};  

d1stencil4th = {{-2,-1, 0, 1, 2}, { 1,-8,  0, 8,-1}, 1/(12h)};
d2stencil4th = {{-2,-1, 0, 1, 2}, {-1,16,-30,16,-1}, 1/(12h^2)};

d1stencil6th = {{-3,-2,-1, 0, 1, 2, 3},
                {-1,  9,-45,   0, 45, -9, 1}, 1/(60h)};
d2stencil6th = {{-3,-2,-1, 0, 1, 2, 3}, 
                { 2,-27,270,-490,270,-27, 2}, 1/(180h^2)};

d1stencil8th = {{-4, -3, -2, -1, 0, 1, 2, 3, 4},
  {3, -32, 168, -672, 0, 672, -168, 32, -3}, 1/(840 h)};

d2stencil8th = {{-4, -3, -2, -1, 0, 1, 2, 3, 4},
  {-9, 128, -1008, 8064, -14350, 8064, -1008, 128, -9}, 1/(5040 h^2)};

d1stencil10th = {{-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5},
  {-2, 25, -150, 600, -2100, 0, 2100, -600, 150, -25, 2}, 1/(2520 h)};

d2stencil10th = {{-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5},
  {8, -125, 1000, -6000, 42000, -73766, 42000, -6000, 1000, -125, 8},
  1/(25200 h^2)};


d4stencil2nd = {   {-2,-1, 0, 1, 2},       { 1, -4,  6, -4, 1},     1/h^4};
d6stencil2nd = {{-3,-2,-1, 0, 1, 2, 3}, { 1,-6, 15,-20, 15, -6, 1}, 1/h^6};
d8stencil2nd = {{-4, -3, -2, -1, 0, 1, 2, 3, 4}, 
  {1, -8, 28, -56, 70, -56, 28, -8, 1}, 1/h^8};
d10stencil2nd = {{-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5},
  {1, -10, 45, -120, 210, -252, 210, -120, 45, -10, 1}, 1/h^10};
d12stencil2nd = {{-6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6}, 
  {1, -12, 66, -220, 495, -792, 924, -792, 495, -220, 66, -12, 1}, 1/h^12};




(* direction dependent finite difference *)
dxyz = {dx, dy, dz};
dijk = {di, dj, dk};
findif[d_, {st_,sw_,sc_}, f_[ijk_]] := 
  (sc /. h->dxyz[[d]]) *
  Sum[ sw[[i]] f[ijk + dijk[[d]] st[[i]]], {i, 1, Length[st]}]
findif2[d_, {st_,sw_,sc_}, d2_, {st2_,sw2_,sc2_}, f_[ijk_]] := 
  (sc /. h->dxyz[[d]]) (sc2 /. h->dxyz[[d2]]) *
  Sum[ sw[[i]] sw2[[j]] f[ijk + dijk[[d]] st[[i]] + dijk[[d2]] st2[[j]]], 
       {i, 1, Length[st]}, {j, 1, Length[st2]}]





d1stencilNthMOS = {{"sMOS0","sMOS1","sMOS2","sMOS3","sMOS4","sMOS5","sMOS6"}, 
                   {"vMOS0","vMOS1","vMOS2","vMOS3","vMOS4","vMOS5","vMOS6"}, 1/(h)};
d2stencilNthMOS = {{"sMOS0","sMOS1","sMOS2","sMOS3","sMOS4","sMOS5","sMOS6"}, 
                   {"VMOS0","VMOS1","VMOS2","VMOS3","VMOS4","VMOS5","VMOS6"}, 1/(h)};
direction = {"i","j","k"};

findifMOS[d_, o_, {st_,sw_,sc_}, f_[ijk_]] := 
  (sc /. h->dxyz[[d]]) *
  Sum[ ToExpression[sw[[i]]<>direction[[d]]] f[ijk + dijk[[d]] ToExpression[st[[i]]<>direction[[d]]]], {i, 1, o}]

findif2MOS[d_, d2_, o_, {st_,sw_,sc_}, f_[ijk_]] := 
  (sc /. h->dxyz[[d]]) (sc /. h->dxyz[[d2]]) *
  Sum[ ToExpression[sw[[i]]<>direction[[d]]] ToExpression[sw[[j]]<>direction[[d2]]] f[ijk + dijk[[d]] ToExpression[st[[i]]<>direction[[d]]] + dijk[[d2]] ToExpression[st[[j]]<>direction[[d2]]]], 
    {i, 1, o}, {j, 1, o}]

MOSdiff = {
    delMOS[m_, x_[ijk]]  :> findifMOS[m, 3, d1stencilNthMOS, x[ijk]],
    delMOS4[m_, x_[ijk]]  :> findifMOS[m, 5, d1stencilNthMOS, x[ijk]],
    deldelMOS[m_, m_, x_[ijk]]  :> findifMOS[m, 3, d2stencilNthMOS, x[ijk]],
    deldelMOS4[m_, m_, x_[ijk]]  :> findifMOS[m, 5, d2stencilNthMOS, x[ijk]],
    deldelMOS[m_, n_, x_[ijk]]  :> findif2MOS[m, n, 3, d1stencilNthMOS, x[ijk]],
    deldelMOS4[m_, n_, x_[ijk]]  :> findif2MOS[m, n, 5, d1stencilNthMOS, x[ijk]]
}


(* standard centered differences *)
centereddiff = {
  del[m_,0] :> 0,
  del4[m_,0] :> 0,
  del6[m_,0] :> 0,
  del8[m_,0] :> 0,
  del10[m_,0] :> 0,
  deldel[m_,n_,0] :> 0,
  deldel4[m_,n_,0] :> 0,
  deldel6[m_,n_,0] :> 0,
  deldel8[m_,n_,0] :> 0,
  deldel10[m_,n_,0] :> 0,
  del[m_, x_[ijk]]  :> findif[m, d1stencil2nd, x[ijk]],
  del4[m_, x_[ijk]] :> findif[m, d1stencil4th, x[ijk]],
  del6[m_, x_[ijk]] :> findif[m, d1stencil6th, x[ijk]],
  del8[m_, x_[ijk]] :> findif[m, d1stencil8th, x[ijk]],
  del10[m_, x_[ijk]] :> findif[m, d1stencil10th, x[ijk]],
  deldel[m_, m_, x_[ijk]] :> findif[m, d2stencil2nd, x[ijk]],
  deldel4[m_, m_, x_[ijk]]:> findif[m, d2stencil4th, x[ijk]],
  deldel6[m_, m_, x_[ijk]]:> findif[m, d2stencil6th, x[ijk]],
  deldel8[m_, m_, x_[ijk]]:> findif[m, d2stencil8th, x[ijk]],
  deldel10[m_, m_, x_[ijk]]:> findif[m, d2stencil10th, x[ijk]],
  deldel[m_, n_, x_[ijk]] :> findif2[m, d1stencil2nd, n, d1stencil2nd, x[ijk]],
  deldel4[m_, n_, x_[ijk]]:> findif2[m, d1stencil4th, n, d1stencil4th, x[ijk]],
  deldel6[m_, n_, x_[ijk]]:> findif2[m, d1stencil6th, n, d1stencil6th, x[ijk]],
  deldel8[m_, n_, x_[ijk]]:> findif2[m, d1stencil8th, n, d1stencil8th, x[ijk]],
  deldel10[m_, n_, x_[ijk]]:> findif2[m, d1stencil10th, n, 
                                      d1stencil10th, x[ijk]]
}


(* dissipation term constructed from high-order derivatives (centered) *)
dissipationdiff = {
  dissipation4[0] :> 0,
  dissipation6[0] :> 0,
  dissipation8[0] :> 0,
  dissipation10[0] :> 0,
  dissipation12[0] :> 0,
  dissipation4[x_[ijk]] :> 
    - Sum[findif[m, d4stencil2nd, x[ijk]] dxyz[[m]]^3, {m,1,3}],
  dissipation6[x_[ijk]] :>
    + Sum[findif[m, d6stencil2nd, x[ijk]] dxyz[[m]]^5, {m,1,3}],
  dissipation8[x_[ijk]] :>
    - Sum[findif[m, d8stencil2nd, x[ijk]] dxyz[[m]]^7, {m,1,3}],
  dissipation10[x_[ijk]] :>
    + Sum[findif[m, d10stencil2nd, x[ijk]] dxyz[[m]]^9, {m,1,3}],
  dissipation12[x_[ijk]] :>
    - Sum[findif[m, d12stencil2nd, x[ijk]] dxyz[[m]]^11, {m,1,3}],
    
  dissipation4v[m_, x_[ijk]] :>  
    - findif[m, d4stencil2nd, x[ijk]] dxyz[[m]]^3,
  dissipation6v[m_, x_[ijk]] :>  
    - findif[m, d6stencil2nd, x[ijk]] dxyz[[m]]^5,
  dissipation8v[m_, x_[ijk]] :>  
    - findif[m, d8stencil2nd, x[ijk]] dxyz[[m]]^7,
  dissipation10v[m_, x_[ijk]] :>  
    - findif[m, d10stencil2nd, x[ijk]] dxyz[[m]]^9,
  dissipation12v[m_, x_[ijk]] :>  
    - findif[m, d12stencil2nd, x[ijk]] dxyz[[m]]^11
} 


(* differences using weights determined by for example advection 
   seems to be obsolete
*)
weighteddiff = {
  wdel[1, u_[ijk]] :> 
  w1M*u[ccM] + w1m*u[ccm] + w1c*u[ccc] + w1p*u[ccp] + w1P*u[ccP], 
  wdel[2, u_[ijk]] :> 
  w2M*u[ccM] + w2m*u[ccm] + w2c*u[ccc] + w2p*u[ccp] + w2P*u[ccP], 
  wdel[3, u_[ijk]] :> 
  w3M*u[ccM] + w3m*u[ccm] + w3c*u[ccc] + w3p*u[ccp] + w3P*u[ccP] 
}


(* advection with precomputed weights and indices *)

advectiondiff = {
  adv[0] :> 0,

  adv[u_[ijk]] :>
          advu0 u[advi0] + advv0 u[advj0] + advw0 u[advk0] +
          advu1 u[advi1] + advv1 u[advj1] + advw1 u[advk1] +
          advu2 u[advi2] + advv2 u[advj2] + advw2 u[advk2],

  adv2[u_[ijk]] :>
          advu0 u[advi0] + advv0 u[advj0] + advw0 u[advk0] +
          advu1 u[advi1] + advv1 u[advj1] + advw1 u[advk1] +
          advu2 u[advi2] + advv2 u[advj2] + advw2 u[advk2],

  adv4[u_[ijk]] :>
          advu0 u[advi0] + advv0 u[advj0] + advw0 u[advk0] +
          advu1 u[advi1] + advv1 u[advj1] + advw1 u[advk1] +
          advu2 u[advi2] + advv2 u[advj2] + advw2 u[advk2] +
          advu3 u[advi3] + advv3 u[advj3] + advw3 u[advk3] +
          advu4 u[advi4] + advv4 u[advj4] + advw4 u[advk4],

  adv6[u_[ijk]] :>
          advu0 u[advi0] + advv0 u[advj0] + advw0 u[advk0] +
          advu1 u[advi1] + advv1 u[advj1] + advw1 u[advk1] +
          advu2 u[advi2] + advv2 u[advj2] + advw2 u[advk2] +
          advu3 u[advi3] + advv3 u[advj3] + advw3 u[advk3] +
          advu4 u[advi4] + advv4 u[advj4] + advw4 u[advk4] +
          advu5 u[advi5] + advv5 u[advj5] + advw5 u[advk5] +
          advu6 u[advi6] + advv6 u[advj6] + advw6 u[advk6],
    
  adv8[u_[ijk]] :>
          advu0 u[advi0] + advv0 u[advj0] + advw0 u[advk0] +
          advu1 u[advi1] + advv1 u[advj1] + advw1 u[advk1] +
          advu2 u[advi2] + advv2 u[advj2] + advw2 u[advk2] +
          advu3 u[advi3] + advv3 u[advj3] + advw3 u[advk3] +
          advu4 u[advi4] + advv4 u[advj4] + advw4 u[advk4] +
          advu5 u[advi5] + advv5 u[advj5] + advw5 u[advk5] +
          advu6 u[advi6] + advv6 u[advj6] + advw6 u[advk6] +
          advu7 u[advi7] + advv7 u[advj7] + advw7 u[advk7] +
          advu8 u[advi8] + advv8 u[advj8] + advw8 u[advk8],

  adv10[u_[ijk]] :>
          advu0 u[advi0] + advv0 u[advj0] + advw0 u[advk0] +
          advu1 u[advi1] + advv1 u[advj1] + advw1 u[advk1] +
          advu2 u[advi2] + advv2 u[advj2] + advw2 u[advk2] +
          advu3 u[advi3] + advv3 u[advj3] + advw3 u[advk3] +
          advu4 u[advi4] + advv4 u[advj4] + advw4 u[advk4] +
          advu5 u[advi5] + advv5 u[advj5] + advw5 u[advk5] +
          advu6 u[advi6] + advv6 u[advj6] + advw6 u[advk6] +
          advu7 u[advi7] + advv7 u[advj7] + advw7 u[advk7] +
          advu8 u[advi8] + advv8 u[advj8] + advw8 u[advk8] +
          advu9 u[advi9] + advv9 u[advj9] + advw9 u[advk9] +
          advu10 u[advi10] + advv10 u[advj10] + advw10 u[advk10]
}


(* works for arbitrary vector with more multiplications *)
(*
  adv[v_, u_[ijk]] :>
    glue[v,1][ccc] (wcaa u[ccc] + wmaa u[maa] + wMaa u[Maa]) +
    glue[v,2][ccc] (waca u[ccc] + wama u[ama] + waMa u[aMa]) +
    glue[v,3][ccc] (waac u[ccc] + waam u[aam] + waaM u[aaM])
*)


(* temporary: first compute linear index stencils, then translate to 
   local node format: xyz = ijk *)
linindtonode = {
  +di +dj +ijk -> ppc,
  -di +dj +ijk -> mpc,
  +di -dj +ijk -> pmc,
  -di -dj +ijk -> mmc,

  +di +dk +ijk -> pcp,
  -di +dk +ijk -> mcp,
  +di -dk +ijk -> pcm,
  -di -dk +ijk -> mcm,

  +dj +dk +ijk -> cpp,
  -dj +dk +ijk -> cmp,
  +dj -dk +ijk -> cpm,
  -dj -dk +ijk -> cmm,

  +di +ijk -> pcc,
  -di +ijk -> mcc,
  +dj +ijk -> cpc,
  -dj +ijk -> cmc,
  +dk +ijk -> ccp,
  -dk +ijk -> ccm,

  ijk -> ccc
}



(* collect and rename some common factors *)
colfactors = {
  x_ / (2dx) -> x oo2dx,
  x_ / (2dy) -> x oo2dy,
  x_ / (2dz) -> x oo2dz,
  Power[dx, -2] -> oodx2,
  Power[dy, -2] -> oody2,
  Power[dz, -2] -> oodz2,
  1 / (4 dx dy) -> oo4dxdy,
  1 / (4 dx dz) -> oo4dxdz,
  1 / (4 dy dz) -> oo4dydz,
  -1 / (4 dx dy) -> -oo4dxdy,
  -1 / (4 dx dz) -> -oo4dxdz,
  -1 / (4 dy dz) -> -oo4dydz,
  Sqrt[2] -> sqrt2,
  Power[2, Rational[-1, 2]] -> oosqrt2,

  (* unfinished, these are left over from time derivatives *)
  dx -> 1/oo2dx/2,  
  dy -> 1/oo2dy/2,  
  dz -> 1/oo2dz/2
}
(*  maybe later we also want:  x_ / 12 -> x oo12,  x_ / 144 -> x oo144  *)



(* MathTensor uses OD *)
ODtodel = {
  OD[u_, a_] :> del[a, x],
  OD[u_, a_, b_] :> deldel[a, b, u]
}



(* time derivatives are a special case *)
timederiv = {
  del[0, x_] :> oo2dt (glue["p", x] - glue["m", x]),
  deldel[0, 0, x_] :> oodt2 (glue["p", x] + glue["m", x] - 2 x),
  deldel[0, n_, x_] :> oo2dt (del[n, glue["p", x]] - del[n, glue["m", x]])
}
timederivdot = {
  del[0, x_] :> glue["dot", x],
  deldel[0, 0, x_] :> 0,
  deldel[0, n_, x_] :> del[n, glue["dot", x]]
}




(* put variables on grid *)
gridvars = Join[vars];
putvarsgrid = Map[(# -> (#)[ijk])&, gridvars]
dif = all /. putvarsgrid;
If[False, prlist["dif = all /. putvarsgrid", dif]]



(* difference equations *)
If[TensorEquationsDim === 4, 
  If[DoTimeDerivs =!= False, dif = dif //. timederiv, 
                             dif = dif //. timederivdot]
]

dif = dif //. MOSdiff 
dif = dif //. centereddiff 
dif = dif //. advectiondiff
dif = dif //. dissipationdiff

dif = DeleteCases[dif, True]



dif = dif //. colfactors
If[False, prlist["dif = dif/.putvarsgrid/.diffrule", dif]]

timer["after discretize.m"]


