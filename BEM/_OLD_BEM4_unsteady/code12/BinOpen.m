%% 
disp('Loading Stochastic wind...')
file='../data/windPWTA1012';
fid =fopen(file,'rb');
U=fread(fid,'float32');
fclose(fid);

file='../data/windPWTA1022';
fid =fopen(file,'rb');
V=fread(fid,'float32');
fclose(fid);

file='../data/windPWTA1032';
fid =fopen(file,'rb');
W=fread(fid,'float32');
fclose(fid);
%%
N1=1024;N2=128;N3=128;
MU=zeros(N1,N2,N3);
MV=zeros(N1,N2,N3);
MW=zeros(N1,N2,N3);
ii=1;
for i = 1:N1
    for j= 1:N2
        for k=1:N3
            MU(i,j,k)=U(ii);
            MV(i,j,k)=V(ii);
            MW(i,j,k)=W(ii);
            ii=ii+1;
        end
    end
end
disp('Done.')
%%
U=0;
V=0;
W=0;
%%
Z=linspace(0,70,N2);
Y=linspace(0,70,N2);
X=linspace(0,500,N1);
% %%
% MeshYZ=meshgrid(Y,Z);
% figure(1)
% for i=1:256
%     clf;
%     surface(Y,Z,squeeze(MU(i,:,:)),squeeze(MU(i,:,:)))
%     zlim([-3 3])
%     view(-35,45)
%     pause(0.5)
% end
% %%  
% MeshXY=meshgrid(X,Y);
% figure(1)
% for i=1:N2
%     clf;
%     surface(Y(1:20),X,squeeze(MU(:,1:20,i)),squeeze(MU(:,1:20,i)),'LineStyle','none')
%     zlim([-3 3])
%     view(-35,45)
%     pause(0.5)
% end

%
