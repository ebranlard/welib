function testIntegration
global PATH
global VERSIONNUM
global VERSION
require('WTlib',VERSION,1);


vx=linspace(0,1,100);
vx=vx(2:end-1);

vy=vx.*(vx-0.5);

[Q]=getTorqueFromBlade(vx,vy,1);
[T]=getThrustFromBlade(vx,vy,1);
assertElementsAlmostEqual(Q ,0.080825338323027 , 'relative',10^-8);
assertElementsAlmostEqual(T ,0.080850078171780, 'relative',10^-8);


% the interface will and should change
