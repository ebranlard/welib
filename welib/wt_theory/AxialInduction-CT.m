InitClear;

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
ac=0.2;
aac=a(a>=ac);
CTgs(a>=ac)=4*(ac^2 +(1-2*ac).*aac);
ac=0.293;
aac=a(a>=ac);
CTgs2=CT;
CTgs2(a>=ac)=4*(ac^2 +(1-2*ac).*aac);

%wind energy explained
CTwee=CT;
CTwee(a>0.4)=((a(a>0.4)-0.143).^2-0.0203+0.6427*0.889)/0.6427;

figure
hold on
grid on
box on
plot(a,CT,'k-','LineWidth',2)
plot(a,CTg,'b--','LineWidth',2)
plot(a,CTwee,'b:','LineWidth',2)
plot(a,CTgs,'b-','LineWidth',2)
plot(a,CTgs2,'b-','LineWidth',2)
plot(a,Cp,'k--','LineWidth',2)
xlabel('a [.]')
ylabel('C_T,C_p')
legend('C_T','C_T Glauert','C_T Glauert empirical','C_T Spera , a_c=0.2','C_T','C_p',0)
title('CTCpInduction')
%%

load 'P2MW.mat'


plot(P2MW(:,1),P2MW(:,2),'+')

%%
V=P2MW(:,1);
P=P2MW(:,2)*1000;
Cp=P./(0.5*1.125*V.^3*5027 )
lambda=16.7*2*pi/60*39./V;
figure
plot(lambda,Cp,'+')
figure
plot(V,Cp,'+')

%%

VO=10;
gamma=2/3*V0;
R=40;
r=0;
ui=gamma*R/(2*(R-r)).*(1+x./sqrt((R-r)^2+x.^2));
V=V0+ui


