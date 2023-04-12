(* ::Package:: *)

fcen[ss_]:= N[Total[ss Range[Length[ss]]]/Total[ss]]
fmax[ss_]:=First[ Ordering[ss,-1]]
fmed[ss_]:= Module[{cumspec, acc=Accumulate[ss],half,p},
\[NonBreakingSpace] half = 0.5 acc[[-1]];
\[NonBreakingSpace] cumspec=Transpose[{Range[Length[ss]],0.5(acc+Prepend[Drop[acc,-1],0])}];
\[NonBreakingSpace] p =Flatten[ Position[(#[[2]]>half)& /@ cumspec,True,1,1]][[1]];
\[NonBreakingSpace] cumspec[[p-1,1]]+(cumspec[[p,1]]-cumspec[[p-1,1]])(half-cumspec[[p-1,2]])/(cumspec[[p,2]]-cumspec[[p-1,2]])
\[NonBreakingSpace] ]



NormalizeSpectrum[u_,p_]:=
     p/Integrate[
		Interpolation[Transpose[{u,p}],InterpolationOrder->3][x]
		,{x,Min[u],Max[u]}
	];


getStatsFromPDF[u_,p_]:=Module[{xmean},
	xmean=Integrate[
		Interpolation[Transpose[{u,u*p}],InterpolationOrder->1][x]
		,{x,Min[u],Max[u]}
	];
	{xmean
	,
		Sqrt[Integrate[
			Interpolation[Transpose[{u,(u-xmean)^2*p}],InterpolationOrder->1][x]
			,{x,Min[u],Max[u]}
		]]
	}
]



PDFInterp[p_]:=Table[{x,Interpolation[p,InterpolationOrder->1][x]},{x,0,19,0.05}];
CumDistr[p_]:=Table[{p[[i,1]],Trapz[p[[Range[i],1]],p[[Range[i],2]]]},{i,1,Length[p]}];
