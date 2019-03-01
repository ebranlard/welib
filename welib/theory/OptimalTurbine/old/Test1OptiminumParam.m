%% Question 1
InitClear
InitRotor
InitSystemProperties
R=Rotor.R;
r=Rotor.r;

lambda=8;
alpha_d=6;  % alpha(17)=6 deg
[a aprime phi c beta]=getOptimizedParameters(alpha_d,lambda,3);


%% table
[r/R a aprime c phi beta]



%% plotting
figure(14)
hold on
plot([0.1 1], [0.3333 0.3333],'r-.')
plot(r/R,a,'r');
plot(r/R,aprime);

xlabel('r/R')
grid()
box()

legend('Axial induction factor without wake rotation (a=1/3)','Axial induction factor (a)','Tangential induction factor (a`)')
ylabel('Induction factors')
%%
figure(15)
hold on
plot(r/R,phi);
grid()
box()
xlabel('r/R [.]')
ylabel('Flow angle phi [deg]')

figure(17)
hold on
grid()
box()
plot(r/R,c);
xlabel('r/R [.]')
ylabel('Chord [m]')
%ylim([0 0.14])


%aa=0:0.01:0.5;plot(aa,getRelationAX(aa,x(1)))
