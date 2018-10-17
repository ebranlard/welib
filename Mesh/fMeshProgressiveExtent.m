function x=fMeshProgressiveExtent(xstart, xend, dx, p)
    % ]xstart xend]

dx=abs(dx);
s=sign(xend-xstart);
dx=dx*s;
p=1+p;

x=[]; % nasty
xlast=xstart;
xi=xlast+dx*p;
dx_last=dx*p;
while abs(xi)<abs(xend);
    x=[x xi];
    xlast=xi;
    xi=xlast+dx_last*p;
    dx_last=dx_last*p;

end


x=sort(x);
