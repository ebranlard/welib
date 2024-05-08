""" 
2D and 3D source point

Reference: 
  [1] Branlard (2017) Wind turbine aerodynamics and vorticity-based methods
"""
import numpy as np


# --------------------------------------------------------------------------------}
# --- 3D 
# --------------------------------------------------------------------------------{
def sp_u(X, Y, Z, Ps, Sigma=1, regParam=0, regMethod=None): 
    """ 
    3D point source
    """
    pass
#     DX = X - Pv[0]
#     DY = Y - Pv[1]
#     r2 = DX ** 2 + DY ** 2
#     tX = -DY
#     tY =  DX
# function [ ui vi wi grad] = fUi_PointSourceN(Xcp,Ycp,Zcp,Sigmas,Xs,Ys,Zs,bComputeGrad)
# % low level implementation
# % takes flat vectors and returns flat vectors
# 
# % [ ui vi wi] = fUi_PointSourceN(1:4,4:7,7:10,1:6,0:5,0:5,0:5) 
# 
# 
# ncp=length(Xcp); %number of control points
# ns=length(Xs); %number of Sources
# 
# 
# Xs=[Xs(:) Ys(:) Zs(:)];              % Source coordinates
# SIG=[Sigmas(:) Sigmas(:) Sigmas(:)]; % Sources intensities repeated for vectorial multiplication
# 
# Ui=zeros(ncp,3); % columns for each component of velocity
# if bComputeGrad 
#     grad=zeros(ncp,9);
# else
#     grad=0;
# end
# 
# %% The vectorial way
# for icp=1:ncp
#     x0=[Xcp(icp) Ycp(icp) Zcp(icp)];
#     X0=repmat(x0,ns,1); % repmat to make use of vectorial computation
#     X0XS=X0-Xs;
#     nX0XS=sum((X0XS.*X0XS)').^(3/2); %line vector
#     Ui(icp,:)=sum(SIG.*(X0XS)./repmat(nX0XS',1,3),1)/(4*pi);  % specify the dimension of summation otherwise for one particle it won't work
# 
#     if bComputeGrad 
#         % below it's for a blob, for the point source, it will be something similar
# %         C=cross(OM,XCPXP);
# %         Cx=C(:,1);
# %         Cy=C(:,2);
# %         Cz=C(:,3);
# %         grad(icp,-3+4) = 1/(4*pi)* sum(                         -3*XCPXP(:,1) .* Cx ./ nXCPXP'.^5); % 11
# %         grad(icp,-3+5) = 1/(4*pi)* sum( -OM(:,3) ./ nXCPXP'.^3  -3*XCPXP(:,2) .* Cx ./ nXCPXP'.^5); % 12
# %         grad(icp,-3+6) = 1/(4*pi)* sum( +OM(:,2) ./ nXCPXP'.^3  -3*XCPXP(:,3) .* Cx ./ nXCPXP'.^5); % 13
# % 
# %         grad(icp,-3+7) = 1/(4*pi)* sum( +OM(:,3) ./ nXCPXP'.^3  -3*XCPXP(:,1) .* Cy ./ nXCPXP'.^5); % 21
# %         grad(icp,-3+8) = 1/(4*pi)* sum(                         -3*XCPXP(:,2) .* Cy ./ nXCPXP'.^5); % 22
# %         grad(icp,-3+9) = 1/(4*pi)* sum( -OM(:,1) ./ nXCPXP'.^3  -3*XCPXP(:,3) .* Cy ./ nXCPXP'.^5); % 23
# % 
# %         grad(icp,-3+10)= 1/(4*pi)* sum( -OM(:,2) ./ nXCPXP'.^3  -3*XCPXP(:,1) .* Cz ./ nXCPXP'.^5); % 21
# %         grad(icp,-3+11)= 1/(4*pi)* sum( +OM(:,1) ./ nXCPXP'.^3  -3*XCPXP(:,2) .* Cz ./ nXCPXP'.^5); % 23
# %         grad(icp,-3+12)= 1/(4*pi)* sum(                         -3*XCPXP(:,3) .* Cz ./ nXCPXP'.^5); % 22
#     end
# end
# ui=Ui(:,1)';
# vi=Ui(:,2)';
# wi=Ui(:,3)';
# 
# 
# --------------------------------------------------------------------------------}
# --- 2d point source 
# --------------------------------------------------------------------------------{
def sp2d_u(X, Y, Ps, Sigma=1, U0=0, regParam=0, regMethod=None, grad=False): 
    """ 
    Induced velocity by one 2D source point on multiple Control Points (CP)
    INPUTS:
    - CP: position of control point
    - Ps: position of source

    EXAMPLE:
        # [Ui]= sp2d_u(1, 0, [0, 0], 2*np.pi)
        # [Ui]= sp2d_u(0, 1, [0, 0], 2*np.pi)

    """
    DX = X - Ps[0]
    DY = Y - Ps[1]
    r2 = DX ** 2 + DY ** 2
    rX = DX
    rY = DY
    #e_t = np.array([-DP[1], DP[0]])
    with np.errstate(divide='ignore', invalid='ignore'):
        if regMethod is None:
            U = Sigma/(2*np.pi) * rX / r2   # [1] Eq. 32.3
            V = Sigma/(2*np.pi) * rY / r2 
        else:                  
            U = Sigma/(2*np.pi) * rX / r2  * (1 - np.exp(- r2 / regParam ** 2))
            V = Sigma/(2*np.pi) * rY / r2  * (1 - np.exp(- r2 / regParam ** 2))
    U += U0
    return U, V

# % TODO: below it's optimized for many sources and few control points, formulas can be nicer for many CPs and loop on sources
# SIG=repmat(Sigmas(:)',1,ndim); % Sources intensities repeated for vectorial multiplication
# % loop on control point
# for icp=1:ncp
#     x0=CP(icp,:);
#     X0=repmat(x0,ns,1); % repmat to make use of vectorial computation
#     X0XS=X0-XS;  % direction from source to CP
#     n2X0XS=sum((X0XS.*X0XS)'); %line vector of the squared norm
#     % ui= sigma /(2*pi nr) r/nr = sigma * r/(2*pi nr^2))
#     Ui(icp,:)=sum(SIG.*(X0XS)./repmat(n2X0XS',1,ndim),1)/(2*pi);  % specify the dimension of summation otherwise for one particle it won't work
# end


def sp2d_psi(x, y, Ps, Sigma):
    """ Streamfunction from one 2D source point """
    z = x-Ps[0] + (y-Ps[1]) * 1j
    # streamfunction (imaginary part of the complex potential)
    with np.errstate(divide='ignore', invalid='ignore'):
        F   = Sigma/(2*np.pi)*np.log(z)
    psi = np.imag(F)
    psi = Sigma/(2*np.pi) * np.arctan2(y,x)
    return psi

def sp2d_phi(x, y, Ps, Sigma):
    """ Potential from one 2D source point """
    z = x-Ps[0] + (y-Ps[1]) * 1j
    # streamfunction (imaginary part of the complex potential)
    with np.errstate(divide='ignore', invalid='ignore'):
        F   = np.log(z)
    phi = np.real(F)
    # phi = (Sigma/4/pi) * np.log(x**2+y**2);
    return phi

