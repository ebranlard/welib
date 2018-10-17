function testAeroCoeff
global PATH
global VERSIONNUM
global VERSION
require('WTlib',VERSION,1);

if(VERSIONNUM>1)
    WT=fInitWT('SB1','xblade',PATH.DATA_WT);
else
    WT=fInitWT('SB1','xblade',PATH.DATA_IN);
end
WT=fSetRotorGrid(20,WT);



% not rough
bReInterp=0; bRough=0; bThicknessInterp=0;
ClCdCm = fAeroCoeff(10,WT.Profiles,WT.Rotor.ProfileSet(:,4),WT.Rotor.thickness_rel(4),60,bReInterp,bThicknessInterp,bRough);
assertElementsAlmostEqual(ClCdCm ,[ 1.0922    0.0391   -0.2295] , 'relative',10^-3);


% not rough
bReInterp=0; bRough=0;
ClCdCm = fAeroCoeff(10,WT.Profiles,WT.Rotor.ProfileSet(:,4),WT.Rotor.thickness_rel(4),60,bReInterp,1,bRough);
assertElementsAlmostEqual(ClCdCm ,[ 1.2123 0.0496 -0.2323] , 'relative',10^-3);
% 
% rough
bReInterp=0; bRough=1;
ClCdCm = fAeroCoeff(10,WT.Profiles,WT.Rotor.ProfileSet(:,4),WT.Rotor.thickness_rel(4),60,bReInterp,1,bRough);
assertElementsAlmostEqual(ClCdCm ,[ 0.6481    0.1832   -0.2518 ] , 'relative',10^-3);

% Re interp
bReInterp=1; bRough=0;
ClCdCm = fAeroCoeff(10,WT.Profiles,WT.Rotor.ProfileSet(:,4),WT.Rotor.thickness_rel(4),3,bReInterp,1,bRough);
assertElementsAlmostEqual(ClCdCm ,[1.1701    0.0537   -0.2298] , 'relative',10^-3);


% Re interp and several alpha
bReInterp=1; bRough=0;
ClCdCm = fAeroCoeff(10:12,WT.Profiles,WT.Rotor.ProfileSet(:,4:6),WT.Rotor.thickness_rel(4:6),2:4,bReInterp,1,bRough);
assertElementsAlmostEqual(ClCdCm ,[ 1.1521    0.0554   -0.2274 ; 1.1962    0.0432   -0.2188;  1.3708    0.0404   -0.2084] , 'relative',10^-2);

