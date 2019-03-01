function plotAirfoilForces(CP,N,T,dForce,scale,scale2,flag)
n= length(CP(1,:));
figure(1213)
hold on, box on
axis equal

Vdelta_thickness=(norm(CP(:,1)-CP(:,2)))*2;
plot(CP(1,[1:end 1]),CP(2,[1:end 1]),'k+-','LineWidth',2,'MarkerSize',5)
for(i=1:n)
    vectarrowb(CP(:,i)+N(:,i)*Vdelta_thickness ,CP(:,i)+N(:,i)*Vdelta_thickness+dForce(:,i)*scale2,scale,'b',1,0.005,20);
end
end