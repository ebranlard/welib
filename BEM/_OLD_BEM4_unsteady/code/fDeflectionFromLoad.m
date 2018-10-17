function [uy uz]=fDeflectionFromLoad(py,pz,beta,nu,r,EI1,EI2)
    n=length(r);
    %Initialization includes boundary condition at the tip
    My=zeros(1,n);
    Mz=zeros(1,n);
    Ty=zeros(1,n);
    Tz=zeros(1,n);
    M1=zeros(1,n);
    M2=zeros(1,n);
    ky=zeros(1,n);
    kz=zeros(1,n);

    %loop to get T and M from p
    for i=n:-1:2
        Ty(i-1)=Ty(i)+0.5*(py(i-1)+py(i))*(r(i)-r(i-1));
        Tz(i-1)=Tz(i)+0.5*(pz(i-1)+pz(i))*(r(i)-r(i-1));

        My(i-1)=My(i) - Tz(i)*(r(i)-r(i-1)) - ( 1/6*pz(i-1)+1/3*pz(i) )*( r(i)-r(i-1) )^2;
        Mz(i-1)=Mz(i) + Ty(i)*(r(i)-r(i-1)) + ( 1/6*py(i-1)+1/3*py(i) )*( r(i)-r(i-1) )^2;
    end

    %going from principal directions to y and z
    M1=My.*cos(beta+nu) - Mz.*sin(beta+nu);
    M2=My.*sin(beta+nu) + Mz.*cos(beta+nu);
    k1=M1./EI1;
    k2=M2./EI2;
    kz=-k1.*sin(beta+nu)+k2.*cos(beta+nu);
    ky=k1.*cos(beta+nu)+k2.*sin(beta+nu);

    %Initialization includes boundary condition at the hub
    Thetay=zeros(1,n);
    Thetaz=zeros(1,n);
    uy=zeros(1,n);
    uz=zeros(1,n);

    for i=1:(n-1)
        Thetay(i+1)=Thetay(i) + 0.5*(ky(i+1)+ky(i))*(r(i+1)-r(i));
        Thetaz(i+1)=Thetaz(i) + 0.5*(kz(i+1)+kz(i))*(r(i+1)-r(i));

        uy(i+1)=uy(i) + Thetaz(i)*(r(i+1)-r(i)) + (1/6*kz(i+1)+1/3*kz(i))*(r(i+1)-r(i))^2;
        uz(i+1)=uz(i) - Thetay(i)*(r(i+1)-r(i)) - (1/6*ky(i+1)+1/3*ky(i))*(r(i+1)-r(i))^2;
    end
end