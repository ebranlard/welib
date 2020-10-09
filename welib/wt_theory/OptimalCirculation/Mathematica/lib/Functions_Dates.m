(* ::Package:: *)

(* Creating Dates sequences *)
DateSeqSplit[t0_,tend_,length_]:=Module[{dt},
dt=DateDifference[t0,tend]/(length-1);
Table[DatePlus[t0,i*dt],{i,0,length-1}]
]
DateSeqBy[t0_,tend_,by_]:=Module[{length},
length=First[Floor[DateDifference[t0,tend,"Second"]/by]];
Table[DatePlus[t0,{i*by,"Second"}],{i,0,length}]
]
DateSeqFrom[t0_,by_,length_]:=Module[{},
Table[DatePlus[t0,{i*by,"Second"}],{i,0,length-1}]
]
(* Examples
d1=DateList["May 5 2009-11:10am"]//DateString 
d1=DateList["May 6 2009-11:10am"]//DateString 
DateSeqSplit[d1,d2,3]
DateSeqBy[d1,d2,600]
DateSeqFrom[d1,600,8]  *)


(* Dealing with timestamps *)
Timestamp2Date[d_]:=If[Length[d]==0,DeTimestamp[d],
Table[DateString[{StringInsert[ToString[d[[i]]],"-",{5,7,9,11}],{"Year","-","Month","-","Day","-","Hour24","-","Minute"}}],{i,1,Length[d]}]]

Timestamp[t_]:=FromDigits[
DateString[t,{"Year", "Month","Day","Hour","Minute" }]
]

DeTimestamp[ts_]:=DateString[{Mod[Floor[ts/10^8],10^4],Mod[Floor[ts/10^6],10^2],Mod[Floor[ts/10^4],10^2],Mod[Floor[ts/10^2],10^2],Mod[ts,10^2]}]

