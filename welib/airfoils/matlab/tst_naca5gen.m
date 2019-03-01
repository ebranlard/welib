clear all
clc

iaf.designation='23012';
iaf.n=56;
iaf.HalfCosineSpacing=1;
iaf.wantFile=1;
iaf.datFilePath='./'; % Current folder
iaf.is_finiteTE=0;

af = naca5gen(iaf);


% set(af.hAF,'Marker','none')
hold on

plot(af.x,af.z,'bo-')

plot(af.xU,af.zU,'bo-')
plot(af.xL,af.zL,'ro-')

plot(af.xC,af.zC,'r--')

title(af.name)
axis equal
