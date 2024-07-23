Off[General::spell]
Off[General::spell1]

<< "BarCharts`"; 
<< "PlotLegends`";

<< psi2h.in.m;



CombineColumns[list1_,list2_]:=Module[{len,lendiff},

lendiff=Length@list1-Length@list2;
If[lendiff < 0, Print["list 1 shorter by ", -lendiff, "elements!"]];
If[lendiff > 0, Print["list 1 longer by ", lendiff, "elements!"]];

len=Min[Length@list1,Length@list2];

Table[{list1[[i]],list2[[i]]},{i,1,len}]
];




RK2Int[rhs_,dt_]:=Module[{i, sol},

sol = 0 * rhs; (* create array of correct dimensions *) 

i = 1;
While[i< Length@rhs,

sol[[i+1]]=sol[[i]]+(dt/2)(rhs[[i]]+rhs[[i+1]]);

i=i+1;
];

sol
];



ReadColData[file_,cols_,lines_]:=  (* version that skips first lines_ lines *)
Module[ {inStream,table,form,ii,rec, throwaway},

form = Table[Number,{ii,1,cols}];

inStream=OpenRead[file];

(* we do not need the first lines *)
Do[ throwaway=Read[inStream,String],{lines}];


(* loop over time steps and read data *)
rec={};
table={};

While[!TrueQ[rec==  EndOfFile],

rec=Read[inStream, form
];

If[!TrueQ[rec==  EndOfFile],AppendTo[table, rec]];
];

(* print amount of data read & close stream *)
Print["read ",N[ByteCount@table / 1024], " KByte of data;"];
Close[inStream];

table
];




hFrom0DPsi4[repsi4_,impsi4_]:=Module[{tmp,dt,offset,ilist,im,re,psi4,times},

avg[list_]:=Total@list/Length@list;

times=Map[First,repsi4];

re=Map[Last,repsi4];
im=Map[Last,impsi4];

psi4=re + I im;

dt=times[[2]]-times[[1]];

tmp=RK2Int[psi4,dt];

offset=avg@tmp;
tmp = tmp - offset;
tmp=RK2Int[tmp,dt];

offset=avg@tmp;
tmp=tmp - offset;

CombineColumns[times,tmp]
];





lmFixed[] := Module[{},

 fnames = {rfname, ifname};

 hString   = "h." <> fname;

 psi = Map[Union, Table[ReadColData[fnames[[j]], 2, 1], {j,1,Length@fnames}]];

 times = Map[First, psi[[1]]];
 rpsi  = Map[Last,  psi[[1]]];
 ipsi  = Map[Last,  psi[[2]]];
 cpsi  = rpsi + I ipsi;

 re = CombineColumns[times, Re@cpsi];
 im = CombineColumns[times, Im@cpsi];


 h = Map[Last, hFrom0DPsi4[psi[[1]], psi[[2]]]];

 h3C = Table[{times[[j]], Re@h[[j]], -Im@h[[j]]}, {j,1,Length@times}];

 Print["Writing to File ", hString];
 Export[hString, h3C, "Table"];
];



(* MAIN LOOP  *)
Do[

fname = filenames[[j]];
Print["Reading from files: [ri]", fname];

rfname = "r" <> fname;
ifname = "i" <> fname;

lmFixed[], {j, 1, Length@filenames}];
