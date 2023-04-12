(* ::Package:: *)

Trapz[x_,y_]:=Total[N[Differences[x]*(Drop[y,1]+Drop[y,-1]) /2  ]]
LeastSquare[x_,y_]:=Total[(x-y)^2];
