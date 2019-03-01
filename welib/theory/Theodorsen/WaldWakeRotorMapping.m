%%
R=1;
Rw=1.5;
rh=0.1*R
x=linspace(rh/R,1,5);
xh=rh/R;

xw_linear=(x-xh)/(1-xh);
xw_quad=sqrt((x.^2-xh.^2)/(1-xh.^2));


figure,hold all
plot(x,xw_linear)
plot(x,xw_quad)

%%
% n=10;
% 
% zbar=linspace(0,n);
% a=0.3;
% Rwbar=sqrt((1-a)./(1-a.*(1+zbar./sqrt(1+zbar.^2)))); %from vortex rings / vortex cylinder analogy
% 
% 
% for i=1:n
%     Mr_linear(i,:)=(x-xh)/(1-xh)*Rwbar(i);
%     Mr_quad(i,:)=sqrt((x.^2-xh.^2)/(1-xh.^2))*Rwbar(i);
%     Mz(i,:)=x*0+zbar(i);
% end
% 
% figure,
% plot(zbar,Rwbar);
% 
% figure,hold all
% plot(Mz',Mr_linear');


%%
Mz=[x*0; x*0+3];
Mr=[x; xw_linear*1.3];
Mr_quad=[x; xw_quad*1.3];

figure
hold all
plot(Mz,Mr,'k')

plot(Mz,Mr_quad,'b')


%%
