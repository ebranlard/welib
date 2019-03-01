(* ::Package:: *)

MyLegend[Texts_,LineStyles_,LineWidths_,Colors_,Markers_,MarkersColors_,Scale_]:=Module[{W,H,n,nC,nL,nLW,nMC,nM,h,fs,ff,thR,ImageScale,LinesT,ls},
(* number of elements*)
n=Length[Texts];
nL=Length[LineStyles];
nLW=Length[LineWidths];
nC=Length[Colors];
nM=Length[Markers];
nMC=Length[MarkersColors];

(* Sizes *)
W=Max[Map[StringLength,Texts]]*0.45+1.5+0.5;
H=n+1/(2);
h=(n-Range[0,n]-1/(4));
ImageScale=40;

(* Font options *)
fs=1/W;
(* ff="Monospace"; *)
(* line thickness ratio *)
thR=H/n/W/20;

(* Lines *)
LinesT={Text[""]};
Do[
If[LineWidths[[Mod[i-1,nLW]+1]]>0,
	(*If[LineStyles[[Mod[i-1,nL]+1]]==1,
		ls=Dashing[1],
		ls=LineStyles[[Mod[i-1,nL]+1]]
	];*)
	ls=Dashing[LineStyles[[Mod[i-1,nL]+1]]];
	LinesT=Join[LinesT,	{
CapForm["Round"],
	ls,
	Colors[[Mod[i-1,nC]+1]],
	Thickness[LineWidths[[Mod[i-1,nLW]+1]]*thR],
	Line[{{0.2,h[[i]]},{1.2,h[[i]]}}]
}],Text[""]];
,{i,1,n}];

Framed[
Graphics[
{
 {Thick,White,Rectangle[{0,0},{W,H}]},

(* Texts *)
Table[
Text[
Style[Texts[[i]],FontSize->Scaled[fs],FontProperties->{"FontMonospaced" -> True}]
,{1.5,h[[i]]},{-1,0} ],{i,1,n}]
,
LinesT
,
(* Markers *)
Table[
Text[
Style[Markers[[Mod[i-1,nM]+1]],MarkersColors[[Mod[i-1,nMC]+1]],FontSize->Scaled[fs/2]],{0.7,h[[i]]},{0,0}]
,{i,1,n}]
}
,
ImageSize->{W*ImageScale,H*ImageScale}*Scale,
Frame->False,PlotRange->{{0,W},{0,H}},AspectRatio->Automatic]
]
]

blue=RGBColor[63/255,63/255,153/255];
red=RGBColor[153/255,61/255,113/255];
green=RGBColor[61/255,153/255,86/255];
yellow=RGBColor[152/255,139/255,61/255];





MyListPlot[L_,o___]:=ListPlot[L
,Frame->True 
,GridLines->Automatic
,GridLinesStyle->Directive[Gray,Dotted]
,o
]


MyListLogPlot[L_,o___]:=ListLogPlot[L
,Joined->True
,Frame->True 
,GridLines->Automatic
,GridLinesStyle->Directive[Gray,Dotted]
,o
]


MyListLogLogPlot[L_,o___]:=ListLogLogPlot[L
,Joined->True
,Frame->True 
,GridLines->Automatic
,GridLinesStyle->Directive[Gray,Dotted]
,o
]


MyExport[fname_]:=Export[FigsPath<>fname<>".pdf",ToExpression[fname]];
