(* ::Package:: *)

(* ::Input:: *)
(**)


SpectrumZeroFilled[x_]:=Module[{},
Zeros=0*Range[0,255];
Zeros[[x[[1]]+1]]=x[[2]];
Zeros
]

Resample[x_,fOld_,fNew_,T_]:=Module[{},
Interpolation[ 
Transpose[{  (Range[Length[x]] -1) *(T-1/fOld)/(Length[x]-1),x} ]
] [ Range[0,T-1/fNew,1/fNew] ]
];

ResampleAndAverage[x_,fOld_,fNew_,T_]:=Module[{},
Mean/@Partition[
Interpolation[
Transpose[{  (Range[Length[x]] -1) *(T-1/fOld)/(Length[x]-1),x} ]][Range[0,T-1/fOld,1/(fNew*Floor[fOld/fNew])]],Floor[fOld/fNew]
]
]

sowingBar[{{x0_,x1_},{y0_,y1_}},__]:=(Sow[{(x1+x0)/2,(y1-y0)}];Rectangle[{x0,y0},{x1,y1}])
(*[{x0,y0},{x1,y1}])*)


(* ::Input:: *)
(**)


(* If Wlidar_,Wsonic are returned by getResults* the nBoxes simulations' statistical moment are compared  *)
PlotComparison[Wsonic_,Wlidar_,yext___]:=Module[
{model,modelAx,xProjValues,xProjBins,yProjValues,yProjBins,xm,xM,ym,yM,mx,Mx,my,My,dx,dy,nBins,sigmaLidar,sigmaSonic,meanLidar,meanSonic},
(* Dealing with plot range *)
my=Min[Wlidar];
My=Max[Wlidar];
mx=Min[Wsonic];
Mx=Max[Wsonic];
my=my-0.1(My-my);
My=My+0.1(My-my);
mx=mx-0.1(Mx-mx);
Mx=Mx+0.1(Mx-mx);
xm=mx;
xM=Mx;

If[Quiet[Check[Length[yext]>0,False]],ym=yext[[1]];yM=yext[[2]],ym=my;yM=My];
xm=ym;
xM=yM;


model=LinearModelFit[Transpose[{Wsonic,Wlidar}], x,x];
modelAx=Fit[Transpose[{Wsonic,Wlidar}], {x},x];

(* Density Projections *)
nBins=16;
dx=(xM-xm)/nBins;
dy=(yM-ym)/nBins;
(*
xProjValues=BinCounts[Wsonic,{xm,xM,(xM-xm)/nBins}];
xProjValues=xProjValues/Total[xProjValues]*(yM-ym)/1.5+ym+0.005(yM-ym); 
xProjBins=Table[x,{x,xm+dx,xM,dx}];
yProjValues=BinCounts[Wlidar,{ym,yM,(yM-ym)/nBins}];
yProjValues=yProjValues/Total[yProjValues]*(xM-xm)/1.5+xm+0.005(xM-xm); 
yProjBins=Table[y,{y,ym+dy,yM,dy}];
*)
(* Standard deviations and mean *)
sigmaLidar=StandardDeviation[Wlidar];
sigmaSonic=StandardDeviation[Wsonic];

meanLidar=Mean[Wlidar];
meanSonic=Mean[Wsonic];


Show[
(* Plotting data : *) 
ListPlot[Transpose[{Wsonic,Wlidar}],Joined->False,PlotStyle->Gray,PlotMarkers->{Automatic,6}], 
Plot[x,{x,xm,xM},PlotStyle->Black],
Plot[model["BestFit"], {x,Min[Wsonic], Max[Wsonic]}],

(* Plotting projections : *)
(*ListPlot[Transpose[{xProjBins,xProjValues}],Joined->True,PlotStyle->Green],
ListPlot[Transpose[{yProjValues,yProjBins}],Joined->True,PlotStyle->Red],*)

(*Now the text : *)
(*Graphics[Text[TextCell[Row[
{ExpressionCell[\[Sigma]==SetPrecision[sigmaLidar,3],"Output"]}],"Text"]
,{xm+0.5dx,(yM+ym)/2},{Left,Bottom}],PlotRange->{{xm,xM},{ym,yM}}],
Graphics[Text[TextCell[Row[
{ExpressionCell[\[Sigma]==SetPrecision[sigmaSonic,3],"Output"]}],"Text"]
,{(xM+xm)/2+0.5dx,ym+0.5dy},{Left,Bottom}],PlotRange->{{xm,xM},{ym,yM}}],*)

Graphics[Text[Style[TextCell[Column[
{ExpressionCell[Row[{"f=",SetPrecision[model["BestFit"],3]}],"Output"],
ExpressionCell[Row[{"R"^2,"=",SetPrecision[model["RSquared"],3]}],"Output"],
ExpressionCell[""],
ExpressionCell[Row[{"f=",SetPrecision[modelAx,3]}],"Output"]
(*ExpressionCell[Mean (y)/Mean(x)==SetPrecision[meanLidar/meanSonic,2],"Output"]*)
} ],FontColor->ColorData[3,"ColorList"]][[1]],"Text"]
,{xm+0.5dx,yM-0.5dy},{Left,Top}],PlotRange->{{xm,xM},{ym,yM}}]

,Frame->True
,Axes->False
,GridLines->Automatic
,GridLinesStyle->Directive[Gray,Dotted]
,PlotRange->{{ym,yM},{ym,yM}}
,AspectRatio->1
 ]
];


(* ::Input:: *)
(**)


(* ::Input:: *)
(**)
