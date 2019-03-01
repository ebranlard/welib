%%
% InitClear
% setFigurePath({'figs_aero/'})

n=100;
%----------------------------------------------------------------------------------------------------   
%% Testing Low level functions
%----------------------------------------------------------------------------------------------------   
if false
    %% Karman Treffz
    n=100;
    xc=0.2;
    yc=0.2;
    tau=20;
    [ P PS SS ] = fProfileKarmanTrefftz( xc,yc,tau,n );

    %% Van de Vooren
    n=100;
    t_rel=0.15;
    tau=17;
    flagMethod=2;
    [ P PS SS Cp u v xg yg chord X Y] = fProfileVanDeVooren( t_rel,tau,1,n,flagMethod);



    % [P PS SS TE chord]=fProfileStandardize(X(end:-1:1),Y(end:-1:1),1,0)


    %                 case 'naca4'
    nSecondaryDirection=90; param.number='0012';
    [ P CP N T PS SS D Dcp L ] = fProfileCoordinates('naca4',param, nSecondaryDirection+1);
    P(:,1)=P(:,1)-0.5; % point from [0 1] to [-0.5 0.5]
    P(:,1)=-P(:,1); % due to my convention of axes
    P=[P zeros(nSecondaryDirection+1,1)]; %TODO inconsistency between methods it seems
    size(P)
    %                 case 'naca5'
    % [ P CP N T PS SS D Dcp L ] = fProfileCoordinates('naca5','' ,nSecondaryDirection+1);
    % P(:,1)=P(:,1)-0.5; % point from [0 1] to [-0.5 0.5]
    % P(:,1)=-P(:,1); % due to my convention of axes
    % P=[P zeros(nSecondaryDirection+1,1)]; %TODO inconsistency between methods it seems

    %                 case 'karman-trefftz'
     param.xc=-0.2; param.yc=0; param.tau=20;
    [ P CP N T PS SS D Dcp L ] = fProfileCoordinates('karman-trefftz',param, nSecondaryDirection);
    P(:,1)=P(:,1)-0.5; % point from [0 1] to [-0.5 0.5]
    P(:,1)=-P(:,1); % due to my convention of axes
    P=[P zeros(nSecondaryDirection+1,1)]; %TODO inconsistency between methods it seems
    size(P)
    %                 case 'vandevooren'
    param=[]; param.t_rel=0.15; param.tau=17; 
    [ P CP N T PS SS D Dcp L ] = fProfileCoordinates('vandevooren',param, nSecondaryDirection+1);
    P(:,1)=P(:,1)-0.5; % point from [0 1] to [-0.5 0.5]
    P(:,1)=-P(:,1); % due to my convention of axes
    P=[P zeros(nSecondaryDirection+1,1)]; %TODO inconsistency between methods it seems
    size(P)


    [ P CP N T PS SS D Dcp L ] = fProfileCoordinates('cylinder',[],nSecondaryDirection+1 );
    P=[P zeros(nSecondaryDirection+1,1)]; %TODO inconsistency between methods it seems
    size(P)

end

%%

Files=dir('data/*.dat');

for i=2:length(Files)
    filename=Files(i).name
    [folder,base,ext]=fileparts(filename);
    sProfile=base; param=[];
    M=load(['data/' base ext]);
    X=M(:,1);
    Y=M(:,2);
    %[P PS SS TE chord IPin IPout] =fProfileStandardize(X,Y,[],false,200);
    [P PS SS TE chord IPin IPout] =fProfileStandardize(X,Y,[],false,400);

%     splines=cscvn(P');
%     plot(P(:,1),P(:,2),'ro')
%     plot(P(1,1),P(1,2),'ko')
%     plot(PS(:,1),PS(:,2),'ro')
%     plot(SS(:,1),SS(:,2),'bd')
%     fnplt(splines)
    size(P)
    figure, box on,grid on,axis square, hold on,
    plot(X,Y,'k.')
    fplotDegrad(PS(:,1),PS(:,2),2), fplotDegrad(SS(:,1),SS(:,2),1),xlim([0 1]),title([sProfile ' ' struct2str(param,' ')]);

end



%%
%----------------------------------------------------------------------------------------------------   
%% Testing Higher level function
%----------------------------------------------------------------------------------------------------   
sProfile='cylinder'; param=[]; 
[ P CP N T PS SS D Dcp L] = fProfileCoordinates( sProfile,param,n);
figure, box on,grid on,axis square, fplotDegrad(PS(:,1),PS(:,2),2), fplotDegrad(SS(:,1),SS(:,2),1),xlim([0 1]),title([sProfile ' ' struct2str(param,' ')]);

sProfile='naca4'; param=[]; param.number='4415';
[ P CP N T PS SS D Dcp L] = fProfileCoordinates( sProfile,param,n);
figure, box on,grid on,axis square, fplotDegrad(PS(:,1),PS(:,2),2), fplotDegrad(SS(:,1),SS(:,2),1),xlim([0 1]),title([sProfile ' ' struct2str(param,' ')]);

sProfile='naca5'; param=[]; param.number='23012';
[ P CP N T PS SS D Dcp L] = fProfileCoordinates( sProfile,param,n);
figure, box on,grid on,axis square, fplotDegrad(PS(:,1),PS(:,2),2), fplotDegrad(SS(:,1),SS(:,2),1),xlim([0 1]),title([sProfile ' ' struct2str(param,' ')]);

sProfile='joukowski'; param=[]; param.xc=-0.2; param.yc=0.2;
[ P CP N T PS SS D Dcp L] = fProfileCoordinates( sProfile,param,n);
figure, box on,grid on,axis square, fplotDegrad(PS(:,1),PS(:,2),2), fplotDegrad(SS(:,1),SS(:,2),1),xlim([0 1]),title([sProfile ' ' struct2str(param,' ')]);

sProfile='karman-trefftz'; param=[]; param.xc=-0.2; param.yc=0.2; param.tau=20; 
[ P CP N T PS SS D Dcp L] = fProfileCoordinates( sProfile,param,n);
figure, box on,grid on,axis square, fplotDegrad(PS(:,1),PS(:,2),2), fplotDegrad(SS(:,1),SS(:,2),1),xlim([0 1]),title([sProfile ' ' struct2str(param,' ')]);

sProfile='karman-trefftz'; param=[]; param.xc=0; param.yc=0.2; param.tau=20; 
[ P CP N T PS SS D Dcp L] = fProfileCoordinates( sProfile,param,n);
figure, box on,grid on,axis square, fplotDegrad(PS(:,1),PS(:,2),2), fplotDegrad(SS(:,1),SS(:,2),1),xlim([0 1]),title([sProfile ' ' struct2str(param,' ')]);

%% In Kerwin
sProfile='karman-trefftz'; param=[]; param.xc=-0.1; param.yc=0.0; param.tau=5; 
[ P CP N T PS SS D Dcp L] = fProfileCoordinates( sProfile,param,n);
figure, box on,grid on,axis square, fplotDegrad(PS(:,1),PS(:,2),2), fplotDegrad(SS(:,1),SS(:,2),1),xlim([0 1]),title([sProfile ' ' struct2str(param,' ')]);
Pref=load('data/geom-karman-trefftz_-0.1_0.0_5_.dat');
hold all
plot(Pref(:,1),Pref(:,2),'-')
%%
sProfile='karman-trefftz'; param=[]; param.xc=-0.2; param.yc=0; param.tau=20; 
[ P CP N T PS SS D Dcp L] = fProfileCoordinates( sProfile,param,n);
figure, box on,grid on,axis square, fplotDegrad(PS(:,1),PS(:,2),2), fplotDegrad(SS(:,1),SS(:,2),1),xlim([0 1]),title([sProfile ' ' struct2str(param,' ')]);
[~,s]=struct2str(param); textlegend(s,'North-East');
%%
sProfile='vandevooren'; param=[]; param.t_rel=0.15; param.tau=17; 
[ P CP N T PS SS D Dcp L] = fProfileCoordinates( sProfile,param,n);
figure, box on,grid on,axis square, fplotDegrad(PS(:,1),PS(:,2),2), fplotDegrad(SS(:,1),SS(:,2),1),xlim([0 1]),title([sProfile ' ' struct2str(param,' ')]);
[~,s]=struct2str(param); textlegend(s,'North-East'); 
