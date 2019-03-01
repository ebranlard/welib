%% Question 1
loading;
B=3;
R=3;
lambda=8;
alpha_d=6;  % alpha(17)=6 deg
r= [0.3 0.5 1.0 1.5 2.0 2.5 2.6 2.7 2.8 2.9 2.95 ];
x=r/R*lambda;   % this is actually lambda_r
[a aprime phi c beta]=getOptimizedParameters(Profile,alpha_d,lambda,B,R,r,x);


%% table
[r'/R a' aprime' c' phi'*180/pi beta']



%% fdgdg

figure(11)
hold on
plot(alpha_origin,Cl_origin,'g+')
plot(alpha,Cl,'b')
xlabel('alpha [deg]')
ylabel('C_l [.]')
xlim([-90 90])
legend('Original data','Interpolated data')
grid
box

figure(12)
hold on;
plot(alpha_origin,Cd_origin,'g+')
plot(alpha,Cd,'b')
xlabel('alpha [deg]')
ylabel('C_d [.]')
xlim([-90 90])
legend('Original data','Interpolated data')
grid
box

figure(13)
plot(alpha,Cl./Cd)
xlabel('alpha [deg]')
ylabel('C_l / C_d [.]')
xlim([-10 20])
grid
box

%% Q1 step1 : calculating a

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
plot(r/R,phi*180/pi);
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
