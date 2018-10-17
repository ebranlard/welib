function varargout = wind_rose(D,F,varargin)
%WIND_ROSE   Wind rose of direction and intensity
% 
%   Syntax:
%      [HANDLES,DATA] = WIND_ROSE(D,I,VARARGIN)
%
%   Inputs:
%      D   Directions
%      I   Intensities
%      VARARGIN:
%       -dtype, type of input directions D, standard or meteo,
%            if meteo, the conversion dnew=mod(-90-D,360) is done;
%            if not meteo, standard is used (default)
%       -n, number of D subdivisons
%       -di, intensities subdivisons, default is automatic
%       -ci, percentage circles to draw, default is automatic
%       -labtitle, main title
%       -lablegend, legend title
%       -cmap, colormap [jet]
%       -colors, to use instead of colormap, for each di
%       -quad, Quadrant to show percentages [1]
%       -ri, empty internal radius, relative to size of higher
%            percentage [1/30]
%       -legtype, legend type: 1, continuous, 2, separated boxes [2]
%       -bcolor, full rectangle border color ['none']
%       -lcolor, line colors for axes and circles ['k']
%       -percbg, percentage labels bg ['w']
%       -ax, to place wind rose on pervious axes, the input for ax
%            must be [theax x y width], where theax is the previous
%            axes, x and y are the location and width is the wind
%            rose width relative to theax width (default=1/5)
%       -parent, by default a new axes is created unless parent is
%                given, ex, parent may be a subplot
%       -iflip, flip the intensities as they go outward radially, ie,
%                highest values are placed nearest the origin [{0} 1]
%       -inorm, normalize intensities, means all angles will have 100%
%       -incout, if 0, data outside di limits will not be used [0 {1}]
%
%   Output:
%      HANDLES   Handles of all lines, fills, texts
%      DATA   Wind rose occurences per direction and intensity
%
%   Examle:
if nargin==0
     d=0:10:350;
     D=[];
     V=[];
     for i=1:length(d)
       n=d(i)/10;
       D=[D ones(1,n)*d(i)];
       V=[V 1:n];
     end
%
     figure
     wind_rose(D,V)
%
     figure
     wind_rose(D,V,'iflip',1)
%
     figure
     wind_rose(D,V,'ci',[1 2 7],'dtype','meteo')
%
     % place it on a previous axes:
%      ax=axes;
%      plot(lon,lat)
%      wind_rose(D,V,'ax',[ax x y 1/3])
end
%
%   MMA 26-11-2007, mma@odyle.net
%
%   IEO, Instituto Español de Oceanografía
%   La Coruña, España

%   10-12-2007 - Added varargin ci and n (nAngles removed as input)
%   17-12-2007 - Added varargin ax, colors
%   22-02-2008 - Added varargin dtype
%   08-05-2008 - Bug fix (bar at dir=0 could be incorrect in some cases)
%   14-05-2008 - Added varargin iflip
%   16-06-2008 - Added varargin parent
%   10-06-2009 - Added varargin incout
%   27-04-2010 - Added output DATA
%   17-06-2010 - Bug fix (E(i,end)=length(find(b>=Ag(end-1))),
%                previously was ...b>Ag...). So the percentages where
%                wrong only when using intensities equal to the lower
%                value of the highest intensity subdivision, basically
%                an academic case.

handles=[];

% varargin options:
dtype='meteo';
nAngles=36;
ri=1/30;
quad=1;
legType=2;
percBg='w';
titStr='';
legStr='';
cmap=jet;
colors=[];
Ag=[]; % intensity subdivs.
ci=[]; % percentage circles
lineColors='k';
borderColor='none';
onAxes=false;
iflip=0;
inorm=0;
parent=0;
IncHiLow=1; % include values higher and lower that the limits of Ag.

vin=varargin;
for i=1:length(vin)
  if isequal(vin{i},'dtype')
    dtype=vin{i+1};
  elseif isequal(vin{i},'n')
    nAngles=vin{i+1};
  elseif isequal(vin{i},'ri')
    ri=vin{i+1};
  elseif isequal(vin{i},'quad')
    quad=vin{i+1};
  elseif isequal(vin{i},'legtype')
    legType=vin{i+1};
  elseif isequal(vin{i},'percbg')
    percBg=vin{i+1};
  elseif isequal(vin{i},'labtitle')
    titStr=vin{i+1};
  elseif isequal(vin{i},'lablegend')
    legStr=vin{i+1};
  elseif isequal(vin{i},'cmap')
    cmap=vin{i+1};
  elseif isequal(vin{i},'colors')
    colors=vin{i+1};
  elseif isequal(vin{i},'di')
    Ag=vin{i+1};
  elseif isequal(vin{i},'ci')
    ci=vin{i+1};
  elseif isequal(vin{i},'lcolor')
    lineColors=vin{i+1};
  elseif isequal(vin{i},'bcolor')
    borderColor=vin{i+1};
  elseif isequal(vin{i},'ax')
    ax=vin{i+1};
    try
      onAxes=ax(1);
      onAxesX=ax(2);
      onAxesY=ax(3);
      onAxesR=ax(4);
    catch
      disp(':: cannot place wind rose on axes, bad argument for ax')
      return
    end
  elseif isequal(vin{i},'iflip')
    iflip=vin{i+1};
  elseif isequal(vin{i},'inorm')
    inorm=vin{i+1};
  elseif isequal(vin{i},'parent')
    parent=vin{i+1};
  elseif isequal(vin{i},'incout')
    IncHiLow=vin{i+1};
  end
end

% other options:
% size of the full rectangle:
rs=1.2;
rl=1.7;

% directions conversion:
if isequal(dtype,'meteo')
  D=mod(90-D,360);
end


% angles subdivisons:
D=mod(D,360);
Ay=linspace(0,360,nAngles+1)-0.5*360/nAngles;

% calc instensity subdivisions:
if isempty(Ag)
  % gen Ag:
  f=figure('visible','off');
  plot(F); axis tight;
  yl=get(gca,'ytick');
  close(f)
  dyl=diff(yl); dyl=dyl(1);
  if min(F)<yl(1),   yl=[yl(1)-dyl yl];   end
  if max(F)>yl(end), yl=[yl yl(end)+dyl]; end
  Ag=yl;
end

for i=1:length(Ay)-1
  if i==1
     I=find( (D>=Ay(i) & D<Ay(i+1)) | D>=Ay(end));
  else
    I=find(D>=Ay(i) & D<Ay(i+1));
  end
  b=F(I);

  for j=1:length(Ag)-1
    if j==length(Ag)-1
      J=find(b>=Ag(j) & b<=Ag(j+1)); % include data with last Agg
    else
      J=find(b>=Ag(j) & b<Ag(j+1));
    end
    E(i,j)=length(J);
  end

  if IncHiLow
    E(i,1)=length(find(b<Ag(2)));
    E(i,end)=length(find(b>=Ag(end-1)));
  end
end
b=sum(E,2)/length(D)*100;

% normalize data:
if inorm
  n=sum(E,2);
  for i=1:length(n)
    E(i,:)=E(i,:)/n(i);
  end
  b=100*ones(size(b));
end

% check if has values higher or lower than the Ag limits
hasH=length(find(F>Ag(end)));
hasL=length(find(F<Ag(1)));

% calc number of percentage circles to draw:
if isempty(ci)
  if inorm
    ci=[25 50 75];
    g=120;
    ncircles=3;
  else
    dcircles=[1 2 5 10 15 20 25 30 50];
    ncircles=3;
    d=abs(1./(dcircles/max(b))-ncircles);
    i=find(d==min(d));
    d=dcircles(i(1));
    if d*ncircles<max(b)
      ncircles=ncircles+1;
    end
    ci=[1:ncircles]*d;
    g=ncircles*d;
  end
else
  ncircles=length(ci);
  g=max(max(ci),max(b));
end

% plot axes, percentage circles and percent. data:
if parent
  wrAx=parent;
  set(wrAx,'units','normalized');
else
  wrAx=axes('units','normalized');
end
ri=g*ri;
handles(end+1)=fill([-rs*g rl*g rl*g -rs*g],[-rs*g -rs*g rs*g rs*g],'w',...
                     'EdgeColor',borderColor);
if onAxes
  set(handles(end),'facecolor','none')
end
hold on
handles(end+1)=plot([-g-ri -ri nan ri g+ri nan 0 0 nan 0 0],...
                    [0 0 nan 0 0 nan -g-ri -ri nan ri g+ri],':','color',lineColors);
t0=[0:360]*pi/180;
labs=[];
Ang=[1/4 3/4 5/4 7/4]*pi;
Valign={'top' 'top' 'bottom' 'bottom'};
Halign={'right' 'left' 'left' 'right'};
for i=1:ncircles
  x=(ci(i)+ri)*cos(t0);
  y=(ci(i)+ri)*sin(t0);

  circles(i)=plot(x,y,':','color',lineColors);
  handles(end+1)=circles(i);

  labs(i)=text((ci(i)+ri)*cos(Ang(quad)),(ci(i)+ri)*sin(Ang(quad)),[num2str(ci(i)),'%'],...
      'VerticalAlignment',Valign{quad},'HorizontalAlignment',Halign{quad},...
      'BackgroundColor',percBg,'FontSize',8);
end
handles=[handles labs];

% calc colors:
if isempty(colors)
  cor={};
  for j=1:length(Ag)-1
    cor{j}=caxcolor(Ag(j),[Ag(1) Ag(end-1)],cmap);
  end
else
  cor=colors;
end

% fill data:
n=sum(E,2);
if iflip, E=fliplr(E); end
for i=1:length(Ay)-1
  if n(i)
    t=linspace(Ay(i),Ay(i+1),20)*pi/180;
    r1=ri;
    for j=1:length(Ag)-1
      r2=E(i,j)/n(i) *b(i) +r1;

      x=[r1*cos(t(1)) r2*cos(t) r1*cos(fliplr(t))];
      y=[r1*sin(t(1)) r2*sin(t) r1*sin(fliplr(t))];

      if iflip, jcor=length(Ag)-1-j+1;
      else, jcor=j;
      end

      if E(i,j)>0, handles(end+1)=fill(x,y,cor{jcor}); end
      r1=r2;
    end
  end
end
axis equal
axis off

% uistack has problems in some matlab versions, so:
%uistack(labs,'top')
%uistack(circles,'top')
ch=get(wrAx,'children');
if inorm
  % only bring circles up in inorm case.
  for i=1:length(circles)
    ch(ch==circles(i))=[]; ch=[circles(i); ch];
  end
end
for i=1:length(labs)
  ch(ch==labs(i))=[]; ch=[labs(i); ch];
end
set(wrAx,'children',ch);


% N S E W labels:
bg='none';
args={'BackgroundColor',bg,'FontSize',8};
h(1)=text(-g-ri, 0,'WEST', 'VerticalAlignment','top',   'HorizontalAlignment','left', args{:});
h(2)=text( g+ri, 0,'EAST', 'VerticalAlignment','top',   'HorizontalAlignment','right',args{:});
h(3)=text( 0,-g-ri,'SOUTH','VerticalAlignment','bottom','HorizontalAlignment','left', args{:});
h(4)=text( 0, g+ri,'NORTH','VerticalAlignment','top',   'HorizontalAlignment','left', args{:});
handles=[handles h];

% scale legend:
L=(g*rl-g-ri)/7;
h=(g+ri)/10;
dy=h/3;

x0=g+ri+(g*rl-g-ri)/7;
x1=x0+L;
y0=-g-ri;

if legType==1 % contimuous.
  for j=1:length(Ag)-1
    lab=num2str(Ag(j));
    if j==1 & hasL & IncHiLow
      lab='';
    end
    y1=y0+h;
    handles(end+1)=fill([x0 x1 x1 x0],[y0 y0 y1 y1],cor{j});
    handles(end+1)=text(x1+L/4,y0,lab,'VerticalAlignment','middle','fontsize',8);
    y0=y1;
  end
  if ~ (hasH & IncHiLow)
    handles(end+1)=text(x1+L/4,y0,num2str(Ag(end)),'VerticalAlignment','middle','fontsize',8);
  end
elseif legType==2 % separated boxes.
  for j=1:length(Ag)-1
    lab=[num2str(Ag(j)) ' - ' num2str(Ag(j+1))];
    if j==1 & hasL & IncHiLow
      lab=['<',num2str(Ag(2))];
    end
    if j==length(Ag)-1 & hasH & IncHiLow
      lab=['>=',num2str(Ag(j))];
    end
    y1=y0+h;
    handles(end+1)=fill([x0 x1 x1 x0],[y0+dy y0+dy y1 y1],cor{j});
    handles(end+1)=text(x1+L/4,(y0+dy+y1)/2,lab,'VerticalAlignment','middle','fontsize',8);
    y0=y1;
  end

end

% title and legend label:
x=mean([-g*rs,g*rl]);
y=mean([g+ri,g*rs]);
handles(end+1)=text(x,y,titStr,'HorizontalAlignment','center');

x=x0;
y=y1+dy;
handles(end+1)=text(x,y,legStr,'HorizontalAlignment','left','VerticalAlignment','bottom');

if onAxes
  place_wr(onAxes,wrAx,onAxesX,onAxesY,onAxesR);
end

if nargout>=1
  varargout{1}=handles;
end
if nargout>=2
  varargout{2}=E;
end

function place_wr(ax,ax2,x,y,width)
if nargin < 5
  width=1/5;
end
uax=get(ax,'units');
pax=get(ax,'position');
set(ax,'units',uax)
axXlim=get(ax,'xlim');
axYlim=get(ax,'ylim');

x_ax2=pax(1)+pax(3)*(x-axXlim(1))/diff(axXlim);
y_ax2=pax(2)+pax(4)*(y-axYlim(1))/diff(axYlim);

pax2=get(ax2,'position');
width=pax(3)*width;
height=pax2(4)*width/pax2(3);
pax2=[x_ax2 y_ax2 width height];

if 1
  % place at centre of the wr, not the bottom left corner:
  ax2Xlim=get(ax2,'xlim');
  ax2Ylim=get(ax2,'ylim');
  dx=(0-ax2Xlim(1))/diff(ax2Xlim)*pax2(3);
  dy=(0-ax2Ylim(1))/diff(ax2Ylim)*pax2(4);
  x_ax2=x_ax2-dx;
  y_ax2=y_ax2-dy;
  pax2=[x_ax2 y_ax2 width height];
end
set(ax2,'position',pax2)



function cor = caxcolor(val,cax,cmap)
%CAXCOLOR   Caxis color for value
%   Find the color for a given value in a colormap.
%
%   Syntax:
%     COLOR = CAXCOLOR(VALUE,CAXIS,COLORMAP)
%
%   Inputs:
%      VALUE
%      CAXIS   Default is current caxis
%      COLORMAP   Default is current colormap
%
%   Output:
%      COLOR   RGB color vector
%
%   Example:
%      figure
%      pcolor(peaks)
%      color=caxcolor(0);
%      set(gcf,'color',color)
%
%   MMA 28-5-2007, martinho@fis.ua.pt

% Department of Physics
% University of Aveiro, Portugal

if nargin < 3
  cmap = get(gcf,'colormap');
end
if nargin < 2
  cax = caxis;
end

n=size(cmap,1);
i= (val-cax(1))/diff(cax) * (n-1) +1;
a=i-floor(i);
i=floor(i);

i=min(i,n);
i=max(i,1);

if i==n
  cor=cmap(n,:);
elseif i==1
  cor=cmap(1,:);
else
  cor=cmap(i,:)*(1-a) + cmap(i+1,:)*a;
end
