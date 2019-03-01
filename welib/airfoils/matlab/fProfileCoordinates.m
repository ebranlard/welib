function [ P CP N T PS SS D Dcp L ] = fProfileCoordinates( name, param,varargin )
% !!! See file MainProfiloeCoordinatesExamples for examples of how to call this function


% Arguments:
% 1: profile name
% 2: parameters for airfoil generation
% 3: number of points
% 4: boolean: make sure there is a point exactly at the Trailing edge
% 5: make sure there is a point exactly at the LE

% Ouputs: matrices of coordinates n x 2
% P arifoil points: the first and last points are the same (should be something like 1,0) corresponding to the TE. Points go from TE,PS,LE,SS,TE
% PS : points on the pressure side from LE to TE (both TE and LE points are included if they exist)
% SS : points on the suction side from TE to LE (both TE and LE points are included if they exsist)
% !!!! see fProfileStandardize for more on the conventions
% CP, control points between each P
% N Normal at CPs
% T tangent at CPs
% D: distance from TE at CPs startng from TE and along the SS
% L: total circumference

% [ P CP N T PS SS D Dcp L] = fProfileCoordinates( 'naca4','4415');

% [ P CP N T PS SS D Dcp L] = fProfileCoordinates( 'naca5','23012');

% param.xc=-0.2; param.yc=0.2; param.tau=20; [ P CP N T PS SS D Dcp L] = fProfileCoordinates( 'karman-trefftz',param,1');

n=100;
bPointAtTE=1;
bPointAtLE=1;
if nargin>2
    n=varargin{1};
end
if nargin==4
    bPlot=varargin{2};
else
    bPlot=0;
end

switch name
    case 'joukowski'
        tau=0; % ie lambda =2; cupsed trailing edge
        if param.yc<0 || abs(param.xc)>1 || abs(param.yc)>1 
            warning('weird param for Joukowski airfoil');
        end
        [ P PS SS ] = fProfileKarmanTrefftz( param.xc,param.yc,tau,n ) ;
    case 'karman-trefftz'
        if param.yc<0 || abs(param.xc)>1 || abs(param.yc)>1 
            warning('weird param for Karman-Trefftz airfoil');
        end
        [ P PS SS ] = fProfileKarmanTrefftz( param.xc,param.yc,param.tau,n );
    
    case 'vandevooren'
        if param.t_rel>1 
            warning('weird param for Van De Vooren airfoil');
        end
        flagMethod=2;
        chord=1;
        [ P PS SS ] = fProfileVanDeVooren( param.t_rel,param.tau,chord,n,flagMethod );

    case 'cylinder'
        theta=linspace(0,-2*pi,n);
%         thetaSS=[theta(theta>-pi) -pi];
%         thetaPS=[theta(theta<-pi) -pi];
        P=[0.5*cos(theta')+0.5 0.5*sin(theta)'];
%         PS=[0.5*cos(thetaPS')+0.5 0.5*sin(thetaPS)'];
%         SS=[0.5*cos(thetaSS')+0.5 0.5*sin(thetaSS)'];

        [P PS SS TE chord]=fProfileStandardize(P(:,1),P(:,2),1);

    case 'naca4'
        iaf.n=floor((n-1)/2);
%         if mod(n,2)==0
%             iaf.n=floor(n/2)-1;
%         end
        iaf.designation=param.number;
        % designation='0008';
         iaf.HalfCosineSpacing=0; iaf.wantFile=0; iaf.datFilePath='./'; % Current folder
        iaf.is_finiteTE=0;
        af = naca4gen(iaf);
%         P=[af.x af.z];
        [P PS SS TE chord]=fProfileStandardize(af.x,af.z,1);
%         PS=[af.xU af.zU];
%         SS=[af.xL af.zL];
%         P=P(end:-1:1,:);
%         PS=PS(end:-1:1,:);
%         SS=SS(end:-1:1,:);
    case 'naca5'
        iaf.n=floor(n/2);
        if mod(n,2)==0
            iaf.n=floor(n/2)-1;
        end
        iaf.designation=param.number;
        iaf.HalfCosineSpacing=0; iaf.wantFile=0; iaf.datFilePath='./'; % Current folder
        iaf.is_finiteTE=0;
        af = naca5gen(iaf);
%         P=[af.x af.z];
        [P PS SS TE chord]=fProfileStandardize(af.x,af.z,1);
%         PS=[af.xU af.zU];
%         SS=[af.xL af.zL];
%         P=P(end:-1:1,:);
%         PS=PS(end:-1:1,:);
%         SS=SS(end:-1:1,:);
    otherwise
        error('unknown profile');
end
% plot(PS(:,1),PS(:,2),'bo-')
% hold on
% plot(SS(:,1),SS(:,2),'ro-')
% 
% figure,
% hold all
% fplotDegrad(PS(:,1),PS(:,2),1)
% fplotDegrad(SS(:,1),SS(:,2),2)
% fplotDegrad(P(:,1),P(:,2)+1,2)

%% Now for the generation of points, normal
Rot=[0 1 ;-1 0];

n=size(P,1);
CP=zeros(n-1,2); % CP control points!!
T=zeros(n-1,2);  %
N=zeros(n-1,2);
Dcp=zeros(1,n-1);
D=zeros(1,n);

CP=(P(1:(end-1),:)+P(2:end,:))/2; % control points between points
L=0;
D(1)=0;
for i=1:(n-1)
    T(i,:)= P(i+1,:)-P(i,:);
    D(i+1)=D(i)+norm(T(i,:)); % distance parcourue since the LE
    Dcp(i)=D(i)+norm(T(i,:))/2;

    T(i,:)= T(i,:)/norm(T(i,:));
    N(i,:)=T(i,:)*Rot;
end
L=D(end);

if bPlot
    fPlotAirfoilScene(CP,P,N,T,.010)
end



end

