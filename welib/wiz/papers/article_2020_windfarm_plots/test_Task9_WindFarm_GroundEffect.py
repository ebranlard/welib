"""
Reproduce Figure 13 from [2]

References:
    [1] E. Branlard, A. Meyer Forsting (2015), Using a cylindrical vortex model to assess the induction zone in front of aligned and yawedd rotors
    [2] Branlard, Forsting (2020) Assessing the blockage effect of wind turbines and wind farms
using an analytical vortex model

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
# --- Influence of ground effect
# --------------------------------------------------------------------------------{
def ground_effect(test=False):
    """ Ground effect 1 turbine and 11x11 farm"""
    ## Parameters
    U0 = 7
    R = 50
    zTarget = - 5 * R
    ## Ground effect
    if not test:
        print('Ground effect - Influence of hub-height')
    if test:
        vCT0 = np.array([0.95])
        vHubHeight = np.linspace(1,3,3) # IN RADIUS
        nHigh=5
    else:
        vCT0 = np.array([0.95,0.4,0.1])
        vHubHeight = np.linspace(1,3,10) # IN RADIUS
        nHigh=11

    vColrs=[]
    vColrs.append(fColrs(1,4))
    vColrs.append(fColrs(3,4))
    vColrs.append(fColrs(4,4))

    vLW     = np.array([2,2,2])
    vSty    = np.array(['--','-',':'])
    legds   = np.array([])
    Spacing = 6*2*R

    fig = plt.figure()
    ax  = fig.add_subplot(111)

    nxWT    = 1
    nzWT    = 1
    for iCT,CT0 in enumerate(vCT0):
        vuz_wf           = np.zeros(vHubHeight.shape)
        vuz_wf_NO_MIRROR = np.zeros(vHubHeight.shape)
        u25D_th = 1 - 0.5 * (1 - np.sqrt(1 - CT0)) * (1 + zTarget / np.sqrt(R ** 2 + zTarget ** 2))
        for iHH,HH in enumerate(vHubHeight):
            HH=HH*R
            if not test:
                print('nx =%2d - CT0 = %.2f - HubHeight %.1f' % (nxWT,CT0,HH))
            # With ground
            _,_,uz_g = windfarm_gridlayout_CTconst(0,HH,zTarget,R,CT0,U0,nxWT,Spacing,nzWT,Spacing,hub_height=HH,mirror=True)
            vuz_wf[iHH] = (uz_g + U0)/U0
            # Without ground
            _,_,uz_ng = windfarm_gridlayout_CTconst(0,HH,zTarget,R,CT0,U0,nxWT,Spacing,nzWT,Spacing,hub_height=HH,mirror=False)
            vuz_wf_NO_MIRROR[iHH] = (uz_ng + U0)/U0
        DeltaUz = vuz_wf_NO_MIRROR - vuz_wf
        UzNewRef = u25D_th - DeltaUz
#         UzNewRef = vuz_wf
        ax.plot(vHubHeight,(1 - UzNewRef/u25D_th) * 100,color=vColrs[iCT],linestyle=vSty[0],linewidth=vLW[iCT],label='C_T=%.2f - %dx%d'%(CT0,1,1))

    nxWT = nHigh
    nzWT = nHigh
    for iCT,CT0 in enumerate(vCT0):
        vuz_wf           = np.zeros(vHubHeight.shape)
        vuz_wf_NO_MIRROR = np.zeros(vHubHeight.shape)
        u25D_th = 1 - 0.5 * (1 - np.sqrt(1 - CT0)) * (1 + zTarget / np.sqrt(R ** 2 + zTarget ** 2))
        for iHH,HH in enumerate(vHubHeight):
            HH=HH*R
            if not test:
                print('nx =%2d - CT0 = %.2f - HubHeight %.1f' % (nxWT,CT0,HH))
            # With ground
            _,_,uz_g = windfarm_gridlayout_CTconst(0,HH,zTarget,R,CT0,U0,nxWT,Spacing,nzWT,Spacing,hub_height=HH,mirror=True)
            vuz_wf[iHH] = (uz_g + U0)/U0
            # Without ground
            _,_,uz_ng = windfarm_gridlayout_CTconst(0,HH,zTarget,R,CT0,U0,nxWT,Spacing,nzWT,Spacing,hub_height=HH,mirror=False)
            vuz_wf_NO_MIRROR[iHH] = (uz_ng + U0)/U0
        DeltaUz = vuz_wf_NO_MIRROR - vuz_wf
        UzNewRef = u25D_th - DeltaUz
#         UzNewRef = vuz_wf
        ax.plot(vHubHeight,(1- UzNewRef/u25D_th)*100,color=vColrs[iCT],linestyle=vSty[1],linewidth=vLW[iCT],label='C_T=%.2f - %dx%d'%(CT0,nxWT,nzWT))
    ax.text(1.05,1.35,'Spacing=%.0fD'%(Spacing/(2*R)))
    ax.set_xlabel('Hub-Height [R]')
    ax.set_ylabel('(U_{z no ground} - U_{z,ground})/U_{z,th} at 2.5D [%]')
    ax.legend()
    ax.set_title('WindFarmGroundEffect')
    ax.set_xlim([1,3])
    ax.set_ylim([0,1.4])
    # ---

def layout_ct_with_ground(test=False):
    ## Sensitivity study 4 - Squared Wind Farm with ground effect, Influence of spacing and number of turbines
    if not test:
        print('Ground effect - Influence of layout and ct')
    ## Parameters
    U0 = 7
    R = 50
    zTarget = - 5 * R
    if test:
        vCT0 = np.array([0.95])
        vnWT = np.array([5])
        vSpacing = np.linspace(1,12,3)
    else:
        vCT0 = np.array([0.95,0.4,0.1])
        vnWT = np.array([31,11,5])
        vSpacing = np.linspace(1,12,20)

    vColrs=[]
    vColrs.append(fColrs(1,4))
    vColrs.append(fColrs(3,4))
    vColrs.append(fColrs(4,4))
    vLW = np.array([2,1,1])
    vSty = np.array(['-','--',':'])
    HubHeight = 1.5*R

    fig = plt.figure()
    ax  = fig.add_subplot(111)
    ax.plot(np.array([0,12]),np.array([0.5,0.5]),'k-')
    for iCT,CT0 in enumerate(vCT0):
        for iN,nWT in enumerate(vnWT):
            nxWT = nWT
            nzWT = nWT
            u25D_th = 1 - 0.5 * (1 - np.sqrt(1 - CT0)) * (1 + zTarget / np.sqrt(R ** 2 + zTarget ** 2))
            vuz_wf = np.zeros(vSpacing.shape)
            for iS,Spacing in enumerate(vSpacing):
                if not test:
                    print('CT {} - nWT {} - Spacing {}'.format(CT0,nWT,Spacing))
                Spacing = Spacing*2*R
                _,_,uz_g = windfarm_gridlayout_CTconst(0,HubHeight,zTarget,R,CT0,U0,nxWT,Spacing,nzWT,Spacing,hub_height=HubHeight,mirror=True)
                vuz_wf[iS] = (uz_g + U0)/U0
            if iN > 1 and iCT > 1:
                label = ''
            else:
                label = 'C_T=%.2f - %dx%d'%(CT0,nxWT,nzWT)
            ax.plot(vSpacing, (1-vuz_wf/u25D_th)*100, color=vColrs[iCT],linestyle=vSty[iN],linewidth=vLW[iN],label=label)

    ax.text(0.2,3.7,'Hub Height=%.1fR'%(HubHeight/R))
    ax.set_yticks([0,0.5,1,2,3,4])
    ax.set_yticklabels(['0.0%','0.5%','1.0%','2.0%','3.0%','4.0%'])
    ax.set_ylim([0,4.0])
    ax.set_xlim([0,12 ])
    ax.set_xlabel('Turbine spacing [D]')
    ax.set_ylabel('(1-U_z/U_{z,th}) at 2.5D [%]')
    ax.grid(axis='both',linestyle='--',linewidth=0.5,color='k')
    ax.legend()
    ax.set_title('WindFarmLayoutGround')


def main(test=False):
    layout_ct_with_ground(test)
    ground_effect(test)

class Test(unittest.TestCase):
    def test_Article_WindFarm_Ground(self):
        main(test=True)
        plt.close('all')

if __name__ == "__main__":
    main(test=False)
    plt.show()


