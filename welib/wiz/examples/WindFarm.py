# --- General
import unittest
import numpy as np
import matplotlib.pyplot as plt
# --- Local
from welib.wiz.WindTurbine import WindTurbine
from welib.wiz.WindFarm import WindFarm

def main():
    # --- Parameters for this sript
    Model    = 'VC'         # Model used for flow computation ['VC','VCFF','VD','SS']
    Ground   = True         # Add ground effect
    no_wake  = False        # Removes wake - set to True when coupled to a wake model
    R        = 65           # Rotor radius [m]
    h_hub    = 1.5*R        # Hub height [m] 
    CT0      = 0.80         # Turbine CT [-]
    U0       = 10           # Free stream velocity [m/s]
    R        = 65           # Wind turbine radius
    wd       = 0*np.pi/180  # Wind direction [rad]
    ye       = 0*np.pi/180  # Yaw error [rad], note: turbine yaw = WD- YE

    Layout = np.array([     
            [0,-3*R, h_hub],
            [0, 3*R, h_hub],
            ])
    #
    U0_g = [U0*np.cos(wd),U0*np.sin(wd),0] 

    # --- Creating wind farm (list of WT)
    WF = WindFarm(name='2x1')
    for r_hub in Layout:
        WF.append( WindTurbine(R=R, e_shaft_yaw0=[1,0,0],e_vert=[0,0,1], r_hub=r_hub) )
    print(WF)

    # --- Setting WT specific values
    for WT in WF:
        WT.update_yaw_pos(wd-ye)   # Wind turbine yaw in [rad]
        WT.update_wind(U0_g)       # Free stream vector, local to that turbine [m/s]
        WT.update_loading(Ct=CT0)  # Thrust coefficient, based on turbine local free stream [-]

    # --- Flow field (NOTE: set U0 to None to remove main flow) 
    nx=50
    ny=51
    x = np.linspace(-4*R,4*R,nx)
    y = np.linspace(-6*R,6*R,ny)
    [X,Y]=np.meshgrid(x,y)
    Z=Y*0+h_hub
    ux, uy, uz = WF.velocity_field(X, Y, Z, no_wake=no_wake, U0=U0_g, ground=Ground, Model=Model)

    fig=plt.figure()
    ax=fig.add_subplot(111)
    Speed=np.sqrt(ux**2+uy**2)
    Speed[Speed>U0*1.5]=U0*1.5
    im=ax.contourf(X/R,Y/R,Speed)#,levels=30,vmin=0,vmax=1.0)
    for WT in WF:
        Rotor=WT.rotor_disk_points()
        ax.plot(Rotor[0,:]/R,Rotor[1,:]/R,'k--')
    cb=fig.colorbar(im)
    sp=ax.streamplot(x/R,y/R,ux,uy,color='k',linewidth=0.7,density=2)
    ax.set_xlabel('x/R [-]')
    ax.set_ylabel('y/R [-]')
    ax.set_aspect('equal')
    plt.show()

if __name__ == "__main__":
    main()
#     unittest.main()
