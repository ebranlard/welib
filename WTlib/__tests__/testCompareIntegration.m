function testCompareIntegration
% purely an exercise of testing
global PATH
global VERSIONNUM
global VERSION

verref=1;
VERSIONREF='v01';
vermax=[];

%% some checks
if isempty(vermax) vermax=max([VERSIONNUM verref])+1; end
if(VERSIONNUM == verref) fprintf('S'); return; end
if(max(VERSIONNUM)>vermax) fprintf('S'); return; end
if(min(VERSIONNUM)<verref) fprintf('S'); return; end
if(length(VERSIONNUM)==1) % we will compare with verref
    vers=[VERSIONNUM verref];
    versname={VERSION,VERSIONREF};
else
    vers=VERSIONNUM;
    versname=VERSION;
end

%% now the test

% init that are not version specific
vx=linspace(0,1,100);
vx=vx(2:end-1);
vy=vx.*(vx-0.5);


for iv=1:length(vers)
    nver=vers(iv);
    require('WTlib',versname{iv},1);
    %version specific stuff if any
%     switch (nver)
%         case :
%         case :
%         otherwise:
%     end
    Q(iv)=getTorqueFromBlade(vx,vy,1);
    T(iv)=getThrustFromBlade(vx,vy,1);
end

assertElementsAlmostEqual(Q(1) ,Q(2) , 'relative',10^-8);
assertElementsAlmostEqual(T(1) ,T(2) , 'relative',10^-8);



end

