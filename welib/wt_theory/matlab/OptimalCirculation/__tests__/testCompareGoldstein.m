function testCompareOneRing
% purely an exercise of testing
global PATH
global VERSIONNUM
global VERSION
global DEBUG
global PLOT
global STOP
if(DEBUG)
    fprintf('--- testCompareGoldstein --- Version: %s %d \n',VERSION,VERSIONNUM );
end

verref=0; VERSIONREF='v00';
vermax=[];
%% some checks
if isempty(vermax) vermax=max([VERSIONNUM verref])+1; end
if(VERSIONNUM == verref) fprintf('S'); return; end
if(max(VERSIONNUM)>vermax) fprintf('S'); return; end
if(min(VERSIONNUM)<verref) fprintf('S'); return; end
if(length(VERSIONNUM)==1) % we will compare with verref
    vers=[VERSIONNUM verref]; versname={VERSION,VERSIONREF};
else
    vers=VERSIONNUM; versname=VERSION;
end

%% now the test

% init that are not version specific
Vl_bar_inv=[2 15];
vl_bar=1./Vl_bar_inv/(2*pi);
VB=[3];
VN=[100 100 100 100 100];
s=-1; % s=-1: WT - s=+1: Propellers
for iB=1:length(VB)
    B=VB(iB);
    Vr=linspace(0,1,VN(iB)); 
    for il=1:length(Vl_bar_inv)
        l_bar=vl_bar(il);        
        % loop on versions
        for iv=1:length(vers)
            nver=vers(iv);
            require('OPTIMCIRC',versname{iv},~DEBUG);
            %version specific stuff if any
            G(iv,:)=fGoldsteinFactor( l_bar,B,Vr );
        end
        assertElementsAlmostEqual(G(1,:) ,G(2,:) , 'relative',10^-14);
        if(PLOT)
            figure(121212),clf,hold all, plot(Vr,G(1,:),'b-',Vr,G(2,:),'k+')
        end
    end
end



%%
if(DEBUG)
    fprintf('--- End of test ---\n');
end

if(STOP)
    keyboard
end
