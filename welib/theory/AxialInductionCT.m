InitClear;
setFigureFont('15');
setFigurePath({'./'});
setMatFigure(0);

%%
a=0:0.01:1;

a13=a(a>=1/3);
ac=0.2;
aac=a(a>=ac);

CT=4*a.*(1-a);
Cp=4*a.*(1-a).^2;

CTg=CT;
CTgs=CT;
CTg(a>=1/3)=4*a13.*(1-1/4*(5-3*a13).*a13);
% ac=0.2;
% aac=a(a>=ac);
% CTgs(a>=ac)=4*(ac^2 +(1-2*ac).*aac);
ac=0.3;
aac=a(a>=ac);
CTgs2=CT;
CTgs2(a>=ac)=4*(ac^2 +(1-2*ac).*aac);

%% wind energy explained (same as Hibbs and Radkey)
CTwee=CT;
CTwee(a>0.4)=((a(a>0.4)-0.143).^2-0.0203+0.6427*0.889)/0.6427;
% CtHibbs
% CtHibbs=0.889 - 0.444*a + 1.55*a.^2; 

%% Madsen
vCTMadsen = linspace(0,4,100);
k3        = 0.089207;
k2        = 0.054496;
k1        = 0.251163;
% k0        = -0.001701;
k0        = 0;
aMadsen   = k0+k1*vCTMadsen+k2*vCTMadsen.^2 + k3*vCTMadsen.^3;
b=vCTMadsen>2.5;
aMadsen(b)=2.19627925 * (vCTMadsen(b)-2.5)+2.3606623;

% (k0+k1*2.5+k2*2.5^2 + k3*2.5^3)
% k1+2*k2*2.5 + 3*k3*2.5^2

figure
hold on
grid on
box on
plot(a      ,CTg      ,'--' ,'LineWidth',2,'Color',fColrs(1))
plot(a      ,CTwee    ,':'  ,'LineWidth',2,'Color',fColrs(1))
plot(a      ,CTgs2    ,'-.'  ,'LineWidth',2,'Color',fColrs(3)) % Spera
plot(aMadsen,vCTMadsen,'-'  ,'LineWidth',2,'Color',fColrs(2))
% plot(a      ,CTgs     ,'-'  ,'LineWidth',2,'Color',fColrs(1))
% plot(a      ,CtHibbs  ,'-.'  ,'LineWidth',2,'Color',fColrs(4))
plot(a      ,CT       ,'k-' ,'LineWidth',2,'Color','k')
plot(a      ,Cp       ,'k--','LineWidth',2,'Color','k')
xlabel('a [-]')
ylabel('C_t, C_p [-]')
% legend('C_t Glauert','C_t Glauert empirical','C_t Spera (a_c=0.3)','C_t Madsen','C_t STT','C_p STT','Location','NorthWest')
legend('C_t Glauert','C_t Glauert empirical','C_t Spera (a_c=0.3)','C_t Madsen','Location','NorthWest');
% ,'C_t STT','C_p STT','Location','NorthWest')
text(0.73,0.84,'C_t STT','FontSize',15)
text(0.73,0.27,'C_p STT','FontSize',15)

title('CTCpInduction')
% xlim([0 1])
% ylim([0 2.5])
%%

% load 'P2MW.mat'
% 
% 
% plot(P2MW(:,1),P2MW(:,2),'+')

%%
% V=P2MW(:,1);
% P=P2MW(:,2)*1000;
% Cp=P./(0.5*1.125*V.^3*5027 )
% lambda=16.7*2*pi/60*39./V;
% figure
% plot(lambda,Cp,'+')
% figure
% plot(V,Cp,'+')
% 
% %%
% 
% VO=10;
% gamma=2/3*V0;
% R=40;
% r=0;
% ui=gamma*R/(2*(R-r)).*(1+x./sqrt((R-r)^2+x.^2));
% V=V0+ui
% 
% 
