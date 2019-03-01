iaf.designation='4415';
% designation='0008';
iaf.n=30;
iaf.HalfCosineSpacing=1;
iaf.wantFile=1;
iaf.datFilePath='./'; % Current folder
iaf.is_finiteTE=0;

af = naca4gen(iaf);

% plot(af.x,af.z,'bo-')

plot(af.xU,af.zU,'bo-')
hold on
plot(af.xL,af.zL,'ro-')

axis equal
