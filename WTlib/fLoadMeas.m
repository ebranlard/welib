function R=fLoadMeas(file,Blade)
R.dat=load(file);
R.psiT=R.dat(:,1);
R.IT=1:length(R.psiT);
R.Fnc=R.dat(:,[4 6 8 10 12]);
R.Ftc=R.dat(:,[4 6 8 10 12]+1);
r_bar_ref=[0.25 0.35 0.60 0.82 0.92];
R.vr=r_bar_ref*2.25;


R.ipsi0   = whichvalue(R.psiT,0) ; 
R.ipsi90  = whichvalue(R.psiT,90) ; 
R.ipsi180 = whichvalue(R.psiT,180) ; 
R.ipsi270 = whichvalue(R.psiT,270) ; 
R.ipsiT=[R.ipsi0 R.ipsi90  R.ipsi180 R.ipsi270 ];
