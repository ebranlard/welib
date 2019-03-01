function fPlotAirfoilScene(CP,P,N,T,scale)
n= size(P,1);
figure(12),clf
hold on
axis equal
plot(P(:,1),P(:,2),'r.--','LineWidth',1,'MarkerSize',13)
plot(CP(:,1),CP(:,2),'k+-','LineWidth',2,'MarkerSize',13)

% quiver(CP(:,1),CP(:,2),N(:,1),N(:,2))
% quiver(CP(:,1),CP(:,2),T(:,1),T(:,2))
quiver([CP(:,1);CP(:,1)],[CP(:,2);CP(:,2)],[N(:,1);T(:,1)],[N(:,2);T(:,2)])
% for(i=1:(n-1))
%       vectarrow(CP(i,:) ,CP(i,:)+N(i,:),scale,'b',1);
%       vectarrow(CP(i,:) ,CP(i,:)+T(i,:),scale,'b',1);
%    end
% end
