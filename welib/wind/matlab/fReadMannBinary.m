function [u]=fReadMannBinary(filename,Nx,Ny,Nz)
N=Nx*Ny*Nz;
try
    fid = fopen(filename,'rb');
    A = fread(fid,N,'single');
catch e
    disp(['error reading: ' filename]);
    rethrow(e);
end
% disp(filename)
% disp(A(1:5))
u=permute(reshape(A,[Nz,Ny,Nx]),[3 2 1]);
