Off[General::spell]
Off[General::spell1]

<<Graphics`Graphics3D`
<<Graphics`Legend`;
<<Graphics`Animation`
<<Graphics`MultipleListPlot`

<<NumericalMath`NIntegrateInterpolatingFunct`

<<Statistics`DataManipulation`
<<Statistics`NonlinearFit`
<<Statistics`LinearRegression`

<< MapLookup.m;  (* Part of Kranc, requires Error.m from Kranc *)  

<< SphereData.m;


Redshift[m_,rsource_,robserve_]:=Sqrt[(1-2 m/rsource)/(1-2m/robserve)];

LumDistance[m_,r_]:=r  (1+ m/(2r))^2;

Psi4Corr[r_,m_]:= (LumDistance[m,r]/r) *(1-2 m/LumDistance[m,r]);

TortoiseR[r_,m_]:=r + 2 m Log[r/(2m) -1];

TimeDelay[r1_,r2_,m_]:=TortoiseR[r2,m]-TortoiseR[r1,m];

AvgOfSliceNumber[data_, slicenum_] := Module[{slice},

      slice = SliceValues@GridSlice[data, slicenum];

      Total[ (Flatten@slice) ]/ (Length@Flatten@slice)
      ];


Print["Finished Init"];

<< sphere.in.m;

Do[
 CompoundExpression[
  name2D = names2D[[j]];

  Timing[data2D = ReadSphereData[dirnames[[1]] <> "/" <> name2D];];

  Print["Finished Reading Data"];


  data = CombineColumns[Times /. data2D, 
      Table[AvgOfSliceNumber[data2D, i], {i, 1, Length[Times /. data2D]}]   ];

  idata = CombineColumns[Times /. data2D, SurfaceIntegrate[data2D] ];   

  Export[name2D <> ".dat", data, "Table"];
  Export[name2D <> ".I.dat", idata, "Table"];
], {j, 1, Length@names2D} ];
