function testGoldstein
global PATH
global VERSIONNUM
global VERSION
global DEBUG
global PLOT
global STOP
if(DEBUG)
    fprintf('--- testGoldstein --- Version: %s %d \n',VERSION,VERSIONNUM );
end


require('OPTIMCIRC',VERSION,~DEBUG);

% should be a line vector
vx=[0 0.1 0.5 0.9 1];
G=fGoldsteinFactor(1/12,3,vx);
assertEqual(G(:)',G);

% zeros at extremities
vx=[0 0.1 0.5 0.9 1];
G=fGoldsteinFactor(1/12,3,vx);
assertEqual([G(1) G(end)],[0 0]);


%% test vs Tibery
vlbar=[1 1/2 1/4 1/8];
vB=[4 3 2];
N=100;  % Number of helical vortex per blade
% vx=(1/N:1/N:1);   % [.] 
vx=linspace(0,1,N);
R=1;
for il=1:length(vlbar)
    l_bar=vlbar(il);
    h=2*pi*R*l_bar;
    for iB=1:length(vB)
        B=vB(iB);
%         fprintf('B %d - 1/l %d\n',B,1/l_bar);
        GammaGoldstein(iB,:) = fGoldsteinFactor( l_bar,B,vx);
        [rt kt] = fTiberyWrench( l_bar , B );
        G=interp1(vx,GammaGoldstein(iB,:),rt); % interpolate on points where tibery is known
        assertElementsAlmostEqual(G,kt,'relative',0.07); % 7% error that's not the best...

    end  
end




%%
if(DEBUG)
    fprintf('--- End of test ---\n');
end

if(STOP)
    keyboard
end
