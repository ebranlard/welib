function damp=getDampingMatrix(Mass,eigenFREQ,delta)
damp=zeros(11,11);
ind=1;
GM=diag(Mass(3:end,3:end));
for i=3:11
     if ind>3
         ind=1;
     end
     damp(i,i)=GM(i-2)*2*pi*eigenFREQ(ind)*delta/pi;
     ind=ind+1;
end
 