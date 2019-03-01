function R=fLoadShen(file,Blade)
load('data/ShenFn35')
load('data/ShenFn60')
load('data/ShenFn82')
load('data/ShenFn92')
load('data/ShenFt35')
load('data/ShenFt60')
load('data/ShenFt82')
load('data/ShenFt92')
npsi=20;
R.psiT=linspace(0,360,npsi);
R.IT=1:length(R.psiT);
R.Fnc=zeros(npsi,5);
R.Ftc=zeros(npsi,5);


R.Fnc(:,1)=70;
A=ShenFn35; R.Fnc(:,2)=interp1(A(:,1),A(:,2),R.psiT,'linear','extrap');
A=ShenFn60; R.Fnc(:,3)=interp1(A(:,1),A(:,2),R.psiT,'linear','extrap');
A=ShenFn82; R.Fnc(:,4)=interp1(A(:,1),A(:,2),R.psiT,'linear','extrap');
A=ShenFn92; R.Fnc(:,5)=interp1(A(:,1),A(:,2),R.psiT,'linear','extrap');
A=ShenFt35; R.Ftc(:,2)=interp1(A(:,1),A(:,2),R.psiT,'linear','extrap');
A=ShenFt60; R.Ftc(:,3)=interp1(A(:,1),A(:,2),R.psiT,'linear','extrap');
A=ShenFt82; R.Ftc(:,4)=interp1(A(:,1),A(:,2),R.psiT,'linear','extrap');
A=ShenFt92; R.Ftc(:,5)=interp1(A(:,1),A(:,2),R.psiT,'linear','extrap');


r_bar_ref=[0.25 0.35 0.60 0.82 0.92];
R.vr=r_bar_ref*2.25;


% R.ipsi0   = whichvalue(R.psiT,0) ; 
% R.ipsi90  = whichvalue(R.psiT,90) ; 
% R.ipsi180 = whichvalue(R.psiT,180) ; 
% R.ipsi270 = whichvalue(R.psiT,270) ; 
% R.ipsiT=[R.ipsi0 R.ipsi90  R.ipsi180 R.ipsi270 ];
