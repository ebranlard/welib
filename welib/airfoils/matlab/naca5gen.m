function af = naca5gen(iaf)
% 
% "naca5gen" Generates the NACA 5 digit airfoil coordinates with desired no.
% of panels (line elements) on it.
%      Author : Divahar Jayaraman (j.divahar@yahoo.com)
% 
% INPUTS-------------------------------------------------------------------
%       iaf.designation = NACA 5 digit designation (eg. '23012') - STRING !
%                 iaf.n = no of panels (line elements) PER SIDE (upper/lower)
% iaf.HalfCosineSpacing = 1 for "half cosine x-spacing" 
%                       = 0 to give "uniform x-spacing"
%          iaf.wantFile = 1 for creating airfoil data file (eg. 'naca2412.dat')
%                       = 0 to suppress writing into a file
%       iaf.datFilePath = Path where the data  file has to be created
%                         (eg. 'af_data_folder/naca5digitAF/') 
%                         use only forward slash '/' (Just for OS portability)
% 
% OUTPUTS------------------------------------------------------------------
% Data:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
%       af.x = x cordinate (nx1 array)
%       af.z = z cordinate (nx1 array)
%      af.xU = x cordinate of upper surface (nx1 array)
%      af.zU = z cordinate of upper surface (nx1 array)
%      af.xL = x cordinate of lower surface (nx1 array)
%      af.zL = z cordinate of lower surface (nx1 array)
%      af.xC = x cordinate of camber line (nx1 array)
%      af.zC = z cordinate of camber line (nx1 array)
%    af.name = Name of the airfoil
%  af.header = Airfoil name ; No of panels ; Type of spacing
%              (eg. 'NACA23012 : [50 panels,Uniform x-spacing]')
% 
% 
% File:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% First line : Header eg. 'NACA23012 : [50 panels,Half cosine x-spacing]'
% Subsequent lines : (2*iaf.n+1) rows of x and z values
% 
% Typical Inputs:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% iaf.designation='23012';
% iaf.n=56;
% iaf.HalfCosineSpacing=1;
% iaf.wantFile=1;
% iaf.datFilePath='./'; % Current folder
% iaf.is_finiteTE=0;


% % [[Calculating key parameters-----------------------------------------]]
cld=str2num(iaf.designation(1))*(3/2)/10;
p=0.5*str2num(iaf.designation(2:3))/100;
t=str2num(iaf.designation(4:5))/100;

a0= 0.2969;
a1=-0.1260;
a2=-0.3516;
a3= 0.2843;

if iaf.is_finiteTE ==1
    a4=-0.1015; % For finite thick TE
else
    a4=-0.1036;  % For zero thick TE
end

% % [[Giving x-spacing---------------------------------------------------]]
if iaf.HalfCosineSpacing==1
    beta=linspace(0,pi,iaf.n+1)';
    x=(0.5*(1-cos(beta))); % Half cosine based spacing
    iaf.header=['NACA' iaf.designation ' : [' num2str(2*iaf.n) 'panels,Half cosine x-spacing]'];
else
    x=linspace(0,1,iaf.n+1)';
    iaf.header=['NACA' iaf.designation ' : [' num2str(2*iaf.n) 'panels,Uniform x-spacing]'];
end

yt=(t/0.2)*(a0*sqrt(x)+a1*x+a2*x.^2+a3*x.^3+a4*x.^4);

P=[  0.05     0.1     0.15    0.2     0.25  ];
M=[  0.0580   0.1260  0.2025  0.2900  0.3910];
K=[361.4     51.64   15.957   6.643   3.230 ];

m=spline(P,M,p);
k1=spline(M,K,m);

xc1=x(find(x<=p));
xc2=x(find(x>p));
xc=[xc1 ; xc2];

if p==0
    xu=x;
    yu=yt;

    xl=x;
    yl=-yt;
    
    zc=zeros(size(xc));
else
    yc1=(1/6)*k1*( xc1.^3-3*m*xc1.^2+m^2*(3-m)*xc1 );
    yc2=(1/6)*k1*m^3*(1-xc2);
    zc=(cld/0.3)*[yc1 ; yc2];

    dyc1_dx=(1/6)*k1*( 3*xc1.^2-6*m*xc1+m^2*(3-m) );
    dyc2_dx=repmat((1/6)*k1*m^3,size(xc2));
    dyc_dx=[dyc1_dx ; dyc2_dx];
    theta=atan(dyc_dx);

    xu=x-yt.*sin(theta);
    yu=zc+yt.*cos(theta);

    xl=x+yt.*sin(theta);
    yl=zc-yt.*cos(theta);
end
af.name=['NACA ' iaf.designation];

af.x=[flipud(xu) ; xl(2:end)];
af.z=[flipud(yu) ; yl(2:end)];

indx1=1:min( find(af.x==min(af.x)) );  % Upper surface indices
indx2=min( find(af.x==min(af.x)) ):length(af.x); % Lower surface indices
af.xU=af.x(indx1); % Upper Surface x
af.zU=af.z(indx1); % Upper Surface z
af.xL=af.x(indx2); % Lower Surface x
af.zL=af.z(indx2); % Lower Surface z

af.xC=xc;
af.zC=zc;

lecirFactor=0.8;
af.rLE=0.5*(a0*t/0.2)^2;

le_offs=0.5/100;
dyc_dx_le=(1/6)*k1*( 3*le_offs.^2-6*m*le_offs+m^2*(3-m) );
theta_le=atan(dyc_dx_le);
af.xLEcenter=af.rLE*cos(theta_le);
af.yLEcenter=af.rLE*sin(theta_le);

% % [[Writing iaf data into file------------------------------------------]]
if iaf.wantFile==1
    F1=iaf.header;
    F2=num2str([af.x af.z]);
    F=strvcat(F1,F2);
    fileName=[iaf.datFilePath 'naca' iaf.designation '.dat'];
    dlmwrite(fileName,F,'delimiter','')
end