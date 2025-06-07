"""




References:
    [1] E. Branlard, A. Meyer Forsting (2015), Using a cylindrical vortex model to assess the induction zone in front of aligned and yawedd rotors
"""
# --- General
import unittest
import matplotlib.pyplot as plt
import numpy as np
# --- Local
from welib.vortilib.elements.VortexCylinder import vc_tang_u
from welib.wiz.WindFarm       import windfarm_gridlayout_CTconst


# --------------------------------------------------------------------------------}
# --- Helper functions 
# --------------------------------------------------------------------------------{
def fColrs(v,m):
    return np.array([0.55,0.55,0.55])*(v-1)/(m-1)

def hline(ax,y,*args,**kwargs):
    ax.plot(ax.get_xlim(),[y,y],*args,**kwargs)
def vline(ax,y,*args,**kwargs):
    ax.plot([x,x],ax.get_ylim(),*args,**kwargs)

# --------------------------------------------------------------------------------}
# --- Influence of layout and CT 
# --------------------------------------------------------------------------------{
def layout_ct(test=False):
    U0      = 7    
    R       = 50   
    zTarget = -5*R # 2.5D                                 ; 

    # Sensitivity study 1 - Squared Wind Farm, Influence of spacing and number of turbines
    if test:
        vCT0     = [0.95, 0.1]
        vnWT     = [5]
        vSpacing = np.linspace(1,12,5)
    else:
        vCT0     = [0.95, 0.4, 0.1]
        vnWT     = [31, 11, 5]
        vSpacing = np.linspace(1,12,20)

    vColrs=[]
    vColrs.append(fColrs(1,4))
    vColrs.append(fColrs(3,4))
    vColrs.append(fColrs(4,4))

    vLW=[2,1,1];
    vSty=['-','--',':'];
    legds=[];
    fig = plt.figure()
    ax  = fig.add_subplot(111)
    #  Setting up grid
    ax.set_xlim([0,12])
    ax.set_ylim([0,1.5])
    hline(ax,1   ,color='k',linestyle='-',linewidth=0.5)
    hline(ax,0.5 ,color='k',linestyle='-',linewidth=0.5)
    hline(ax,0.40,color='k',linestyle=':',linewidth=0.5)
    hline(ax,0.30,color='k',linestyle=':',linewidth=0.5)
    hline(ax,0.20,color='k',linestyle=':',linewidth=0.5)
    hline(ax,0.10,color='k',linestyle=':',linewidth=0.5)
    ax.grid(axis='x',linestyle='-',linewidth=0.5,color='k')
    # --- Loop on CT and spacing 
    for iCT,CT0 in enumerate(vCT0):
        if not test:
            print('CT0 = %.2f'%CT0)
        for iN,nWT in enumerate(vnWT):
            nxWT=nWT
            nzWT=nWT
            if not test:
                print('Turbine grid %d x %d '%(nxWT,nzWT));
            u25D_th = 1-0.5*(1-np.sqrt(1-CT0))*(1+zTarget/np.sqrt(R**2+zTarget**2)) ;
            vuz_wf=np.zeros(vSpacing.shape);
            for iS,Spacing in enumerate(vSpacing):
                if not test:
                    print('Spacing %.1f'%Spacing)
                Spacing=Spacing*2*R
                ux,uy,uz = windfarm_gridlayout_CTconst(0,0,zTarget,R,CT0,U0,nxWT,Spacing,nzWT,Spacing)
                vuz_wf[iS]=(uz+U0)/U0;
            if iN>1 and iCT>1:
                label=None
            else:
                label='C_T=%.2f - %dx%d'%(CT0,nxWT,nzWT);
            ax.plot(vSpacing,(1-vuz_wf/u25D_th)*100,color=vColrs[iCT],linestyle=vSty[iN],linewidth=vLW[iN],label=label);
    ax.set_yticks([0,0.5,1,1.5])
    ax.set_yticklabels(['0.0%','0.5%','1.0%','1.5%'])
    ax.set_xlabel('Turbine spacing [D]')
    ax.set_ylabel('(1-U_z/U_{z,th}) at 2.5D [%]')
    ax.legend()
    ax.set_title('WindFarmInfluenceofSpacing')

def layout(test=False):
    # Sensitivity study 2 - Wind Row, Influence of spacing and number of turbines
    U0      = 7    
    R       = 50   
    zTarget = -5*R # 2.5D                                 ; 
    CT0 = 0.4
    if test:
        vnWT=[1,5];
        vSpacing=np.linspace(1,12,3)
    else:
        vnWT=[1,31];
        vSpacing=np.linspace(1,12,20)

    vColrs=[]
    vColrs.append(fColrs(1,2))
    vColrs.append(fColrs(1,2))
    vLW=[2,2,2,2];
    vSty=['-',':','--','-'];
    fig = plt.figure()
    ax  = fig.add_subplot(111)
    #  Setting up grid
    ax.text(0.2,0.93,'CT=%.2f'%CT0)
    ax.set_xlim([0,12])
    ax.set_ylim([0,1.0])
    hline(ax,1   ,color='k',linestyle='-',linewidth=0.5)
    hline(ax,0.5 ,color='k',linestyle='-',linewidth=0.5)
    hline(ax,0.40,color='k',linestyle=':',linewidth=0.5)
    hline(ax,0.30,color='k',linestyle=':',linewidth=0.5)
    hline(ax,0.20,color='k',linestyle=':',linewidth=0.5)
    hline(ax,0.10,color='k',linestyle=':',linewidth=0.5)
    ax.grid(axis='x',linestyle='-',linewidth=0.5,color='k')
    k=0;
    for iN,nWT in enumerate(vnWT):
        for iin,nnWT in enumerate(vnWT):
            nzWT=nWT
            nxWT=nnWT
            if not test:
                print('Turbine grid %d x %d '%(nxWT,nzWT))
            u25D_th = 1-0.5*(1-np.sqrt(1-CT0))*(1+zTarget/np.sqrt(R**2+zTarget**2)) ;
            vuz_wf=np.zeros(vSpacing.shape)
            for iS,Spacing in enumerate(vSpacing):
                if not test:
                    print('Spacing %.1f'%Spacing);
                Spacing=Spacing*2*R
                ux,uy,uz = windfarm_gridlayout_CTconst(0,0,zTarget,R,CT0,U0,nxWT,Spacing,nzWT,Spacing)
                vuz_wf[iS]=(uz+U0)/U0;
            k=k+1;
            if k>1:
                label='nFront x nSide = %dx%d'%(nxWT,nzWT)
                ax.plot(vSpacing,(1-vuz_wf/u25D_th)*100,color=vColrs[0],linestyle=vSty[k-1],linewidth=vLW[k-1],label=label)
    ax.set_yticks([0,0.5,1])
    ax.set_yticklabels(['0.0%','0.5%','1.0%'])
    ax.set_xlabel('Turbine spacing [D]')
    ax.set_ylabel('(1-U_z/U_{z,th}) at 2.5D [%]')
    ax.legend()
    ax.set_title('WindFarmLayoutInfluenceofSpacing')
    #     text(0.2,0.93,sprintf('CT=%.2f',CT0)) # TODO


def main(test=False):
    layout_ct(test)
    layout(test)

class Test(unittest.TestCase):
    def test_Article_WindFarm_NoGround(self):
        main(test=True)

if __name__ == "__main__":
    main(test=False)
    plt.show()


