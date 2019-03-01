function plotAirfoilForces(CP,N,T,dForce,Force,scale,scale2,scale3,flag)
n= length(CP(1,:));
figure
hold on, box on
axis equal
plot(CP(1,[1:end 1]),CP(2,[1:end 1]),'k-','LineWidth',1.8)
for(i=1:n)
      if(flag)
          if(dot(dForce(:,i),N(:,i))<0)
              vectarrowb(CP(:,i)-dForce(:,i)*scale2*scale ,CP(:,i),1,'r',1,0.02,20);
          else
              vectarrowb(CP(:,i) ,CP(:,i)+dForce(:,i)*scale2,scale,'b',1,0.02,20);
          end
      else 
          if(dot(dForce(:,i),N(:,i))<0)
              vectarrowb(CP(:,i) ,CP(:,i)+dForce(:,i)*scale2,scale,'r',1,0.02,20);
          else
              vectarrowb(CP(:,i) ,CP(:,i)+dForce(:,i)*scale2,scale,'b',1,0.02,20);
          end
      end
end
   vectarrowb([0,0], Force'*scale3,scale,'k',2,0.05,20);
%    vectarrowb([0,0], Lift'*scale3,scale,'b',2,0.05,20);
%    vectarrowb([0,0], Drag'*scale3,scale,'r',2,0.05,20);
xlabel('$x/c$ [.]')
ylabel('$y/c$ [.]')
end