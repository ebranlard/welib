(* ::Package:: *)

MyDrop[L_,I_]:=Module[{LL,II},
II=Sort[I,Greater];
LL=L;
Do[
LL=Drop[LL,{II[[i]]}];
,{i,1,Length[II]}
];
LL
]

WhichInf[x_,x0_]:=Select[Transpose[{Range[Length[x]],x}],#[[2]]<x0&][[All,1]];
WhichSup[x_,x0_]:=Select[Transpose[{Range[Length[x]],x}],#[[2]]>x0&][[All,1]];
WhichEq[x_,x0_]:=Select[Transpose[{Range[Length[x]],x}],#[[2]]==x0&][[All,1]];
WhichDiff[x_,x0_]:=Select[Transpose[{Range[Length[x]],x}],#[[2]]!=x0&][[All,1]];


