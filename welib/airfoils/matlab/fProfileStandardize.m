function [P PS SS TE chord IPin IPout] =fProfileStandardize(X,Y,chord_input,bFlatBack,varargin)
% Function that standardize the profile coordinates:
% - Make sure the coordinates goes in a clockwise direction as: TE PS LE SSTE
% - Make sure the Trailing Edge point (1,0) appears twice
% - Normalize to 1 if required
% - Returns Coordinates of Pressure Side, Suction Side and Trailing Edge (if flat back)
%          P:  airfoil points: the first and last points are (1,0) corresponding to the TE. Points go from TE,PS,LE,SS,TE
%          PS : points on the pressure side from LE to TE (both TE and LE points are included if they exist)
%          SS : points on the suction side from TE to LE (both TE and LE points are included if they exsist)
% - For now this script is not concerned about the number of points inputed
% 
% ^y
% |
% |        .>>>SS>>>>,
% | (0,0)=PLE          PTE=(1,0) (point in double. if flat back, several points)
% |        '<<<<PS<<<'
% | -------------------------->x
% 
%
% Emmanuel Branlard, November 2012

bDebug=0;

%
ninit=length(X);
IPin=1:ninit; % this vectors gives the index so that  [ X(IPin) Y(IPin) ] ~= P(IPout,:)
IPout=[];

% TODO flat airfoil would crash..



%% Remove duplicate points, only TE can be duplicate
% fix round off
Y(abs(Y)<10^-14)=0;
P=[X(:) Y(:)];
[~,I1]=unique(P,'rows');
Isuppr=setdiff(1:ninit,I1);
Iunique=setdiff(1:ninit,Isuppr); % his can probably be done in a better way, but I had to keep the orderring and unique sort
X=P(Iunique,1);
Y=P(Iunique,2);
IPin=IPin(Iunique);

n=length(X);
%% Rotate indexes to avoid starting or finishing by the TE or LE... (avoid dealing with some exceptions..)
imax=find(X==max(X));
imin=find(X==min(X));
istart=floor((imax(1)+imin(1))/2);
Inew=loop(istart+(1:n),n);
X2=X(Inew);
Y2=Y(Inew);
IPin=IPin(Inew);
if bDebug
    figure, fplotDegrad(X,Y,1)
    figure, fplotDegrad(X2,Y2,1)
end
X=X2;
Y=Y2;

%% Make profile rotate clockwise
imax=find(X==max(X)); % TE
imin=find(X==min(X)); % LE
imin=imin(1);
imax=imax(1);
if imin<imax
    I1=imin:imax;
else
    I1=imax:imin;
end
I2=setdiff(1:n,I1);
meanY1=mean(Y(I1));
meanY2=mean(Y(I2));
X3=X;
Y3=Y;
if (imin<imax && meanY1<meanY2) || (imin>imax && meanY1>meanY2)
    X3=X(end:-1:1);
    Y3=Y(end:-1:1);
    IPin=IPin(end:-1:1);
else
    % do nothing
end
if bDebug
    figure, fplotDegrad(X3,Y3,1)
end
X=X3;
Y=Y3;




%% First find the Leading edge Point. Whatever the scaling the leading edge point is assumed to be at Y==0.
iLEsup=[]; iLEinf=[]; iLE=[];
IXLE=find(X==min(X));
iplus =loop(IXLE(1)+1,n);
iminus=loop(IXLE(1)-1,n);
if length(IXLE)==1 && Y(IXLE)==0
    % fine, the LE point is part of the input data
    iLE=IXLE;
else
    if length(IXLE) ==1
        % let's find the index of the positive y and negative y points around the LE (+/- 1 around the IXLE)
        if Y(IXLE)>0;
            iLEsup=IXLE;
            if Y(iplus)<0 ,
                iLEinf=iplus;
            else, 
                iLEinf=iminus;
            end
        else
            iLEinf=IXLE;
            if Y(iplus)>0 ,
                iLEsup=iplus;
            else, 
                iLEsup=iminus;
            end
        end
    else
        % flat LE, we assumed that this is not correct, so we'll interpolate later.
        if Y(IXLE(1))>0;
            iLEsup=IXLE(1);
            iLEinf=IXLE(2);
        else
            iLEsup=IXLE(2);
            iLEinf=IXLE(1);
        end
    end
end

%% Similarly find the trailing edge Point. Whatever the scaling the leading edge point is assumed to be at Y==0.
iTEsup=[]; iTEinf=[]; iTE=[];
IXTE=find(X==max(X));
iplus =loop(IXTE(1)+1,n);
iminus=loop(IXTE(1)-1,n);
if sum(abs(Y(IXTE)))==0
    % fine, the TE point is part of the input data
    iTE=IXTE;
else
    if length(IXTE)==1
        % let's find the index of the positive y and negative y points around the TE (+/- 1 around the IXTE)
        if Y(IXTE)>0;
            iTEsup=IXTE;
            if Y(iplus)<0 ,
                iTEinf=iplus;
            else, 
                iTEinf=iminus;
            end
        else
            iTEinf=IXTE;
            if Y(iplus)>0 ,
                iTEsup=iplus;
            else, 
                iTEsup=iminus;
            end
        end
    else
        % TODO review that
        % flat TE -> This depends now, if flat back or not
        [ymin,iy]=min(abs(Y(IXTE)));
        IYTE=IXTE(iy);
        if ymin==0
            % perfect the TE is part of the input data
            iTE=IYTE;
        else
            % let's find the index of the positive y and negative y points around the TE (+/- 1 around the IYTE)
            iplus =loop(IYTE+1,n);
            iminus=loop(IYTE-1,n);
            if Y(IXTE)>0;
                iTEsup=IYTE;
                if Y(iplus)<0 ,
                    iTEinf=iplus;
                else, 
                    iTEinf=iminus;
                end
            else
                iTEinf=IYTE;
                if Y(iplus)>0 ,
                    iTEsup=iplus;
                else, 
                    iTEsup=iminus;
                end
            end
        end
    end
end


%% If needed, interpolation of the LE
if isempty(iLE)
    I=loop((iLEinf-2):(iLEsup+2),n);
    xLE=interp1(Y(I),X(I),0,'spline');
else
    xLE=X(iLE);
end
if isempty(iTE)
    I=loop((iTEinf-2):(iTEsup+2),n);
    xTE=interp1(Y(I),X(I),0,'spline');
else
    xTE=X(iTE);
end



if xLE>min(X)
    warning('It''s an airfoil whose min point in x is not at y=0')
    xLEmin=min(X);
else
    xLEmin=xLE;
end
if xTE<max(X)
    warning('It''s an airfoil whose max point in x is not at y=0')
    [xTEmax imax]=max(X);
    iTE=imax;
    xTE=xTEmax;

else
    xTEmax=xTE;
end


%% Normalizing coordinates with respect to chord and garanteeing that the profile has a real LE at x=0
if isempty(chord_input)
    chord_input=xTEmax-xLEmin; %should be >0
end
chord=chord_input;
if chord_input<0
    warning('Oh no, I should deal with airfoil in front direction')
    kbd
end
xx=(X-xLEmin)./chord_input;
yy=Y./chord_input;


if min(xx)<0
    warning('bad LE point')
end
if max(xx)>1
    warning('bad TE point')
end
if bDebug
    figure, fplotDegrad(xx,yy,1)
end


%% Small index manip
xxLE=0;
xxTE=(xTE-xLEmin)/chord_input;

xxTE=1;

if isempty(iLE)

else
    PLE=[xxLE 0];
    iLEinf=iLE;
    iLEsup=iLE;
end
if isempty(iTE)
    yyTE=0;
    PTEbefore=[xxTE yyTE];
else
    yyTE=Y(iTE)/chord_input; % should be zero most of the time, but not for bad profiles
    PTEbefore=[];
    iTEinf=iTE;
    iTEsup=iTE-1;
end
PTEafter=[xxTE yyTE];





%%  Now assembling the final airfoil coordinates

Iall=loop(iTEinf-1+(1:n),n);
P=[PTEbefore; X(Iall) Y(Iall); PTEafter];

if ~isempty(PTEbefore)
    IPoutOffset=2;
else
    IPoutOffset=1;
end
if ~isempty(PTEafter)
    IPoutOffset2=1;
else
    IPoutOffset2=0;
end
IPin=IPin(Iall);
IPout=IPoutOffset:(length(Iall)+IPoutOffset2);



if iLEinf> iTEinf
    IPS=loop((iTEinf):(iLEinf),n);
else
    IPS=loop((iTEinf):(iLEinf+n),n);
end
if iLEsup<iTEsup
    ISS=loop(iLEsup:iTEsup,n);
else
    ISS=loop(iLEsup:(iTEsup+n),n);
end
P=[PTEbefore; xx(Iall) yy(Iall); PTEafter];
PS=[PTEbefore; xx(IPS) yy(IPS) ];
SS=[xx(ISS) yy(ISS); PTEafter];
TE=[]; % TODO
if(nargin==5)
    % Interpolation
    % We use a monotonaly increasing function to do the interpolation. We use the curvilinear length along the airfoil
    %
    n=varargin{1};

    s=zeros(length(P(:,1)),1);
    s(1)=0;
    for i=2:length(P(:,1))
        s(i)=sqrt((P(i-1,1)-P(i,1))^2+(P(i-1,2)-P(i,2))^2);
    end
    s=cumsum(s);
    s0=linspace(min(s),max(s),n);
    x0=interp1(s,P(:,1),s0);
    y0=interp1(s,P(:,2),s0);
    P=[x0' y0'];
end

if bDebug
    figure, fplotDegrad(P(:,1),P(:,2),1)
    figure, fplotDegrad(PS(:,1),PS(:,2),2), fplotDegrad(SS(:,1),SS(:,2),1)
end

