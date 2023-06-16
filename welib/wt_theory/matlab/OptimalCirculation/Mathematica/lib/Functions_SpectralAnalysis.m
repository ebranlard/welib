(* ::Package:: *)

(* Double-sided spectra : *)
(* Frequencies are cyclic, not angular!  *)


Spectrum[ll_, n_, SamplFreq_, SpecOnly_:False] :=
 Module[ {lreal = Partition[ll, n], ft},
  ft = (Plus @@ (Map[ (Abs[Fourier[#, FourierParameters->{0,1}]]^2)&, lreal]))/(Length[lreal]*SamplFreq);
  If[ SpecOnly, Take[ft, n/2+1],
    Transpose[ {Range[0,n/2]*N[SamplFreq/n], Take[ft, n/2+1]} ]]  ]

CrossSpectrum[ll1_, ll2_, n_, SamplFreq_, SpecOnly_:False] :=
 Module[ {l1 = Partition[ll1, n], l2 = Partition[ll2, n], xft},
  xft = (Plus @@ (Map[ Fourier[#, FourierParameters->{0,1}]&, l1]*Conjugate[Map[ Fourier[#, FourierParameters->{0,1}]&, l2]]));
  xft = xft/(Length[l1]*SamplFreq);
  If[ SpecOnly, Take[xft, n/2+1],
    Transpose[ {Range[0,n/2]*N[SamplFreq/n], Take[xft, n/2+1]} ]]  ]

Coherence[ll1_, ll2_, n_, SamplFreq_] :=
  Transpose[ {Range[0,n/2]*N[SamplFreq/n],
  Abs[ CrossSpectrum[ ll1, ll2, n, SamplFreq, True] ]^2/(
      Spectrum[ ll1, n, SamplFreq, True]*
      Spectrum[ ll2, n, SamplFreq, True]) } ]

(* In this logarithmic smoothing function the frequencies must be
   equally spaced   *)
LogSmooth[spec_List, npd_Integer] :=
   Module[ {deltaf = spec[[2,1]] - spec[[1,1]],
            nlo, nhi = 0, nmax, lloc,
            spc = If[ spec[[1,1]] == 0.0, Drop[spec,1], spec],
            result},
       nmax = Floor[ npd Log[10, spc[[-1,1]]/spc[[1,1]]] + 1];
       result = Table[ nlo = nhi+1;
                  nhi = If[i == nmax, Length[spc], Floor[N[ 10^(i/npd)]]];
                  If[nlo > nhi, {},
                    lloc = Take[ spc, {nlo, nhi}];
                    {(Plus @@ lloc)/Length[lloc]}], {i, nmax}];
       Flatten[result, 1]
   ]


KSTest[F1_,F2_]:=Module[{MaxTarget},
MaxTarget=Max[Abs[F1[[All,2]] -F2[[All,2]] ] ] ;
Select[Transpose[ { Range[Length[F1[[All,1]]]],F1[[All,1]] , Abs[F1[[All,2]] -F2[[All,2]]  ] }] ,#[[3]]== MaxTarget&]][[1]]
