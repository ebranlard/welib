function K=getStiffnessMatrix(Mass,eigenFREQ,k)
     K=[k 0
         0 0];
     ind=1;
     GM=diag(Mass(3:end,3:end));
     for i=3:11
         if ind>3
             ind=1;
         end
         K(i,i)=GM(i-2).*(2*pi*eigenFREQ(ind))^2;
         ind=ind+1;
     end
end
 

 
