function [S_smooth,f_smooth]=fLogSmooth(vf,vS,n_per_decade)
%  Performs bin averaging per periods

% Dropping zero frequency if present
if(vf(1)==0) 
    vS=vS(2:end);
    vf=vf(2:end);
end

nmax=floor(n_per_decade * log(vf(end)/vf(1))/log(10)+1  );
n_high=0;
f_smooth=zeros(1,nmax);
S_smooth=zeros(1,nmax);
cpt=0;
for i=1:nmax
    n_low=n_high+1;
    if i==nmax
        n_high=length(vf);
    else
        n_high=floor(10^(i/n_per_decade));
    end
    if(n_low>n_high)
        %
    else
        cpt=cpt+1;
        f_loc = sum(vf(n_low:n_high))/(n_high-n_low+1);
        s_loc = sum(vS(n_low:n_high))/(n_high-n_low+1);
        f_smooth(cpt)=f_loc;
        S_smooth(cpt)=s_loc;
    end
end
f_smooth=f_smooth(1:cpt);
S_smooth=S_smooth(1:cpt);


