import numpy as np



def tapered_cylinder_geom(R1, R2, H):
    """ conical taper geometry calculations (volume and center of volume)
    INPUTS: 
     - R1: base radius
     - R2: "top" radius
     - H: cylinder height
    OUTPUTS: 
     - V   : volume of tapered section
     - h_c : center of volume offset from base
    """
    m = (R2-R1)/H
    if R1 == R2: # cylinder
       V   = abs(np.pi*R1*R1*H)
       h_c = H/2.0
    elif R1==0: 
        V   = abs(1.0/3.0*np.pi*R2*R2*H)
        h_c = 3.0/4.0*H
    else:
       V   = abs(np.pi/3.0/m*(R2**3 - R1**3))
       h_c = H*(R1**2 + 2*R1*R2 + 3*R2**2)/4.0/(R1**2 + R1*R2 + R2**2) # ( coneV*1./4.*coneH - coneVtip*(1./4.*(coneH-H) + H) )/ taperV
    return V, h_c  


def tapered_cylinder_inertia(R1, R2, H, rho):
    """ Moments of inertias for a tapered cylinder
    INPUTS:
      - R1
      - R2
      - H
      - rho ! density of material
    OUTPUTS:
      - Il
      - Ir radial moment of inertia about node 1
    """
    m = (R2-R1)/H
    if R1 == R2: # cylinder
        Ir = abs(1.0/12.0* rho*np.pi*R1*R1*H *(3.0*R1*R1 + 4.0*H*H)) # radial inertia about node 1 
        Il = abs(0.5* rho*np.pi*R1*R1*H *R1*R1)    
    elif R1==0:
        Ir = abs(rho*np.pi*(1.0/20.0/m + 1.0/5.0/m**3) * R2**5)      
        Il = abs(1.0/10.0*rho*np.pi/m*R2**5)            
    else:
       h_c = H*(R1**2 + 2*R1*R2 + 3*R2**2)/4.0/(R1**2 + R1*R2 + R2**2) 
       #l_c = R1/M + (R2-R1)/m *(R1**2 + 2*R1*R2 + 3*R2**2)/4/(R1**2 + R1*R2 + R2**2) 
       Ir_tip = abs(np.pi/20.0 *rho/m*(1.0 + 4.0/m**2) * (R2**5 - R1**5))                    # radial moment of inertia about tip of cone
       Ir = abs(Ir_tip - rho/3.0/m*np.pi*(R2**3-R1**3) * (R1/m + 2.0*h_c)*R1/m )  # radial moment of inertia about node 1
       Il = abs(1.0/10.0/m*rho*np.pi*(R2**5 - R1**5))  
    return Ir, Il




def tapered_cylinder_prop_MG(R1, R2, Rmg1, Rmg2, L, rho):
    """ 
    Properties (volume, inertia) of a tapered cylinder with marine growth
    INPUTS:
        - R1
        - R2
        - Rmg1
        - Rmg2
        - L
        - rho:   density of material
    OUTPUTS
        - Vinner : volume from inner radius
        - Vouter : volume from outer radius
        - m_mg   : mass of marine growth
        - h_c    : center of mass offset from first node
        - Ilmg   : moment of inertia about axis
        - Irmg   : moment of inertia about radial axis from first node
    """
    # get V and CV for element
    Vinner, cVinner = tapered_cylinder_geom(R1, R2, L)
    # get V and CV for marine growth displacement
    Vouter, cVouter = tapered_cylinder_geom(Rmg1, Rmg2, L)
    # get mass and CV specific to marine growth thickness
    m_mg = (Vouter - Vinner)*rho
    if m_mg==0:
        h_c = 0.0
    else:
        h_c = (cVouter*Vouter - Vinner*cVinner)/(Vouter - Vinner)
    # get two moments of inertia for marine growth as if solid...
    Ilouter, Irouter = tapered_cylinder_inertia(Rmg1, Rmg2, L, rho) # inertias for marine growth if solid
    Ilinner, Irinner = tapered_cylinder_inertia(R1  , R2  , L, rho) # inertias for element if filled with marine growth
    # subtract to get moments of inertia of marine growth shell
    Ilmg = Ilouter - Ilinner
    Irmg = Irouter - Irinner
    return Vinner, Vouter, m_mg, h_c, Ilmg, Irmg
