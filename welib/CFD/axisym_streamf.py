""" 
Tools to solve the streamfunction equation in axisymmetric coordinates


TODO:
    analytical solution Circular Couette flow
      u_tehta = Ar + B/r
      A = (Omega_in R_in^2 - Omega_out R_out^2)/(R_in^2 - R_out^2)
      B = Rin^2 R_out^2 (Omega_in-Omega_out)/(R_out^2 - R_in^2)

    lid driven cavity flow

"""
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as sciint


def velocity_from_streamfunction(psi, r, dr, dz):
    """
    ur=-1/r dpsi/dz
    uz= 1/r dpsi/dr
    """
    nz,nr= psi.shape
    ur = np.zeros_like(psi)
    uz = np.zeros_like(psi)
    for i in range (nz):
        #for j in range (1,nj-1):
        for j in range (1,nr):
            
            #skip over walls, otherwise differencing on neighbors will be off
            #if (node_type[i,j]==WALL): 
            #    continue

            #ur = -1/r dpsi/dz
            if (i==0):
                #ur[i,j] = -(psi[i+1,j]-psi[i,j])/(dz*r[j])
                # Second order forward
                ur[i,j] = -(-psi[i+2,j]+4*psi[i+1,j]-3*psi[i,j])/(2*dz*r[j])
            elif (i==nz-1):
                #ur[i,j] = -(psi[i,j]-psi[i-1,j])/(dz*r[j])
                # Second order backward
                ur[i,j] = -(3*psi[i,j]-4*psi[i-1,j]+psi[i-2,j])/(2*dz*r[j])
            else:
                ur[i,j] = -(psi[i+1,j]-psi[i-1,j])/(2*dz*r[j])

            #uz = 1/r dpsi/dr            
            if (j==0):
                # r=0
                pass
                #uz[i,j] =  (psi[i,j+1]-psi[i,j])/(dr*r[j])
            elif (j==nr-1):
                #uz[i,j] =  (psi[i,j]-psi[i,j-1])/(dr*r[j])
                # Second order backward
                uz[i,j] =  (3*psi[i,j]-4*psi[i,j-1]+psi[i,j-2])/(2*dr*r[j])
            else:
                uz[i,j] =  (psi[i,j+1] - psi[i,j-1])/(2*dr*r[j])
       
    # axis (j=0): u velocity on the axis from q=2*pi*psi
    uz[:,0] = 2*psi[:,1]/((r[1])**2)
    ur[:,0] = 0
    
    #similar approach to get u velocity on nj-1
    #uz[:,nr-1] = (psi[:,nr-1] - psi[:,nr-2])/dr
    #u[node_type==WALL] = 0

    return ur, uz

BC_Dir_0     = 0000
BC_Dir_Uz0   = 2000
BC_Dir_uz    = 2100
BC_Neu_ur0   = 5000
BC_Neu_ur    = 5100
BC_Neu_uz0   = 5200
BC_Neu_uz    = 5300

def streamfunction(om, r, psi, ur, uz, dr, dz, 
        BC_Walls={'south':BC_Dir_0, 'west':BC_Dir_Uz0, 'north':BC_Neu_uz, 'east':BC_Neu_ur0},
        nIter=300000, rTol=1e-12, verbose=False):
    """ 
    
    """
    psi2 = np.copy(psi)
    nz,nr= psi.shape
    
    idz2 = 1/(dz*dz)
    idr2 = 1/(dr*dr)
    for it in range(nIter):   
        psi2[1:nz-1,1:nr-1] = (om[1:nz-1,1:nr-1]*r[1:nr-1]  
                             + idz2*(psi[0:nz-2,1:nr-1]+psi[2:nz,1:nr-1])
                             + idr2*(psi[1:nz-1,0:nr-2]+psi[1:nz-1,2:nr])
                             - 1/(2*dr*r[1:nr-1])*(psi[1:nz-1,2:nr]-psi[1:nz-1,0:nr-2])
                             )/(2*(idz2+idr2))        
                            
        #replace values on boundary nodes with previous values
        #psi2[node_type>0] = psi[node_type>0]
        # --- Inlet (i=0, j)
        if BC_Walls['west']==BC_Dir_Uz0:
            # - Dirichlet - uz known and constant
            psi2[0,:] = 0.5*uz[0,0]*r**2
        elif BC_Walls['west']==BC_Dir_uz:
            # - Dirichlet - uz known
            psi2[0,1:] = sciint.cumtrapz(r*uz[0,:], r)
        elif BC_Walls['west']==BC_Neu_ur0:
            # - Neumann - No radial flow: ur=0,  dpsi/dz = -r ur = 0 - First order
            #psi2[0,:] = psi2[1,:] 
            # - Neumann - No radial flow: ur=0,  dpsi/dz = -r ur = 0 - Second order
            psi2[0,:] = 1/3*(-psi2[2,:]+4*psi2[1,:])
        elif BC_Walls['west']==BC_Neu_ur:
            # - Neumann - ur known,  dpsi/dz = -ur r  - First order
            #psi2[0,:] = psi2[1,:] + dz*ur[0,:]*r[:]
            # - Neumann - ur known,  dpsi/dz = -ur r  - Second order forward
            psi2[0,:] = 1/3*(-psi2[2,:]+4*psi2[1,:] + 2*dz*ur[0,:]*r[:])
        else:
            raise NotImplementedError('BC_walls west', BC_Walls)
        

        # --- Outlet (i=nz-1, j)
        if BC_Walls['east']==BC_Dir_Uz0:
            # - Dirichlet - uz known and constant
            psi2[nz-1,:] = 0.5*uz[0,0]*r**2
        elif BC_Walls['east']==BC_Dir_uz:
            # - Dirichlet - uz known
            psi2[n1-1,1:] = sciint.cumtrapz(r*uz[nz-1,:], r)
        elif BC_Walls['east']==BC_Neu_ur0:
            # - Neumann - no radial flow ur=0 - First order
            #psi2[nz-1,:] = psi2[nz-2,:]
            # - Neumann - no radial flow ur=0 - Second order
            psi2[nz-1,:] = 1/3*(-psi2[nz-3,:] + 4*psi2[nz-2,:])
        elif BC_Walls['east']==BC_Neu_ur:
            # - Neumann - ur known:  dpsi/dz = -ur r  - First order
            #psi2[nz-1,:] = psi2[nz-2,:] - dz*ur[nz-1,:]*r[:]
            # - Neumann - ur known:  dpsi/dz = -ur r  - Second order
            psi2[nz-1,:] = 1/3*(-psi2[nz-3,:] + 4*psi2[nz-2,:] - 2*dz*ur[nz-1,:]*r[:])
        #for j in range(nr):
        #    if (node_type[nz-1,j]!=WALL):
        #        psi2[nz-1,j] = psi2[nz-2,j] - dz*(v[nz-1,j])*r[j]
        else:
            raise NotImplementedError('BC_walls east', BC_Walls)

        # --- Top boundary (i, j=nr-1)
        if BC_Walls['north']==BC_Dir_Uz0:
            # - Dirichlet - freestream constant velocity u0
            psi2[:,nr-1] = 0.5*uz[0,0]*r[nr-1]**2
        elif BC_Walls['north']==BC_Dir_uz:
            # - Dirichlet - freestream varying velocity
            psi2[:,nr-1] = np.trapz(r*uz[0,:],r)
        elif BC_Walls['north']==BC_Neu_uz:
            # - Neumann - uz known:  dpsi/dr = r uz - First order
            #psi2[:,nr-1] = psi2[:,nr-2] + dr*uz[:,nr-1]*r[nr-1]
            # - Neumann - uz known:  dpsi/dr = r uz - Second order
            psi2[:,nr-1] = 1/3*(-psi2[:,nr-3] + 4*psi2[:,nr-2] + 2*dr*uz[:,nr-1]*r[nr-1])
        elif BC_Walls['north']==BC_Neu_uz0:
            # - Neumann - dpsi/dr = ruz=0 - First order
            #psi2[:,nr-1] = psi2[:,nr-2]
            # - Neumann - dpsi/dr = ruz=0 - Second order
            psi2[:,nr-1] = 1/3*(-psi2[:,nr-3] + 4*psi2[:,nr-2])
        else:
            raise NotImplementedError('BC_walls south', BC_Walls)

        # --- Bottom boundary (i, j=0)
        if BC_Walls['south']==BC_Dir_0:
            # - Dirichlet - psi known to be 0
            psi2[:,0] = 0
        else:
            raise NotImplementedError('BC_walls', BC_Walls)
 
        # copy back solution        
        psi = np.copy(psi2) 
            
        # --- Convergence check
        if (it%25==0):
            R = np.zeros_like(psi)
            R[1:nz-1,1:nr-1] = (om[1:nz-1,1:nr-1]*r[1:nr-1] + 
                    idz2*(psi[0:nz-2,1:nr-1]-2*psi[1:nz-1,1:nr-1]+psi[2:nz,1:nr-1])+
                    idr2*(psi[1:nz-1,0:nr-2]-2*psi[1:nz-1,1:nr-1]+psi[1:nz-1,2:nr])-
                    1/(2*dr*r[1:nr-1])*(psi[1:nz-1,2:nr]-psi[1:nz-1,0:nr-2]))
            #R[node_type>0] = 0
            norm     = np.linalg.norm(R)
            norm_avg = norm/(nz*nr)
            if (norm_avg<rTol): 
                if verbose:
                    print('>>> Converged, number of iterations',it)
                return psi2
            if (it%1000==0) and verbose:
                print('    it: {:07d} norm: {:6.2e}   {:6.2e}'.format(it, norm, norm_avg))
    print('[FAIL] Psi failed to converge, norm = {:6.2e}'.format(norm))
    return psi



#---------------- VORTICITY -------------------------------
def applyVorticityBoundaries(w,psi,ur,uz,dr,dz):
    ni,nj = w.shape
    dz2 = dz*dz
    dr2 = dr*dr    
    #apply boundaries
    # Generation of vorticity on walls
    #for i in range(ni):
    #    for j in range(nj):
    #        count = 0
    #        ww = 0
    #        #left wall
    #        if (i<ni-1 and node_type[i,j]==WALL and node_type[i+1,j]==OPEN):
    #            ww += 2*(psi[i,j]-psi[i+1,j])/(r[j]*dz2) - 2*v[i,j]/dz
    #            count += 1
    #            
    #        #right wall
    #        if (i>0 and node_type[i,j]==WALL and node_type[i-1,j]==OPEN):
    #            ww += 2*(psi[i,j]-psi[i-1,j])/(r[j]*dz2) + 2*v[i,j]/dz
    #            count += 1
    #        
    #        #top wall
    #        if (j>0 and node_type[i,j]==WALL and node_type[i,j-1]==OPEN):
    #            ww += 2*(psi[i,j]-psi[i,j-1])/(r[j]*dr2) - 2*uz[i,j]/dr + uz[i,j]/r[j]
    #            count +=1 
    #        
    #        #set values
    #        if (count>0):
    #            w[i,j] = ww/count
                
    #outlet on right side, dw/dz = 0
    w[ni-1,:] = w[ni-2,:]
    
    #outlet on rmax, dw/dr = 0
    #for i in range(ni-outlet_dz,ni):
    #    for i in range(ni-outlet_dz,ni):
    w[:,nj-1]=w[:,nj-2]
    
    #inlet on left side
    #for j in range(1,inlet_nn):
    #for j in range(1,inlet_nn):
    for j in range(1,nj):
        if (j<nj-1):
            duz_dr = (uz[0,j+1]-uz[0,j-1])/(2*dr)
        else:
            duz_dr = (uz[0,j]-uz[0,j-1])/dr                    
        #w[i,j] = 2*(psi[i,j]-psi[i+1,j])/dz2 - 2*v[i,j]/dz + duz_dr                 
        w[0,j] = -duz_dr
    
    # --- Axis
    w[:,0] = 0

#computes RHS for vorticity equation
def R(w, ur, uz, r, dr, dz, nu):
    ni,nj = w.shape
    dz2 = dz*dz
    dr2 = dr*dr    
    #make copy so we use consistent data
    res = np.zeros_like(w)
    for i in range(1,ni-1):
        for j in range(1,nj-1):
            #if (node_type[i,j]>0): continue
                
            #viscous term, d^2w/dz^2+d^2w/dr^2+(1/r)dw/dr
            A = nu*(
                (w[i-1][j]-2*w[i][j]+w[i+1][j])/dz2 + 
                (w[i][j-1]-2*w[i][j]+w[i][j+1])/dr2 +
                (w[i][j+1]-w[i][j-1])/(2*dr*r[j]))
            
            #convective term uz*dw/dz    
            B = uz[i][j]*(w[i+1][j]-w[i-1][j])/(2*dz)
            
            #convective term ur*dw/dr
            C = ur[i][j]*(w[i][j+1]-w[i][j-1])/(2*dr)

            res[i][j] = A - B - C
    return res
            
#advances vorticity equation using RK4
def advanceRK4(dt, w, psi, ur, uz, r, z, dr, dz, nu):
    applyVorticityBoundaries(w,psi,ur,uz,dr,dz)
        
    #compute the four terms of RK4
    Rk = R(w, ur, uz, r, dr, dz, nu)
    w1 = w + 0.5*dt*Rk
    
    R1 = R(w1, ur, uz, r, dr, dz, nu)
    w2 = w + 0.5*dt*R1
    
    R2 = R(w2, ur, uz, r, dr, dz, nu)
    w3 = w + dt*R2
    
    R3 = R(w3, ur, uz, r, dr, dz, nu)
    w_new = w + (dt/6.0)*(Rk + 2*R1 + 2*R2 +R3) 
              
    #return new value
    return w_new


if __name__ == '__main__':
    pass
    #psi = streamfunction(om, r, psi, ur0, uz0, dr, dz)
    #ur,uz = velocity_from_streamfunction(psi, r, dr, dz)


