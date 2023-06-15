# --- For cmd.py
import os
import subprocess
import multiprocessing

import collections
import glob
import pandas as pd
import numpy as np
import distutils.dir_util
import shutil 
import stat
import re


def spanwiseBD(tsAvg,vr,R,postprofile=None, IR=None):
    # --- Extract radial data
    Columns=[]
    for sB in ['B1','B2','B3']:
        Columns.append(extractSpanTSReg(tsAvg,'^'+sB+'N(\d)TDxr_\[m\]',sB+'TDxr_[m]'))
        Columns.append(extractSpanTSReg(tsAvg,'^'+sB+'N(\d)TDyr_\[m\]',sB+'TDyr_[m]'))
        Columns.append(extractSpanTSReg(tsAvg,'^'+sB+'N(\d)TDzr_\[m\]',sB+'TDzr_[m]'))

    dfRad, nrMax, ValidRow = _HarmonizeSpanwiseData('BeamDyn', Columns, vr, R, IR=IR)

    # --- Export to csv
    if postprofile is not None and dfRad is not None:
        dfRad.to_csv(postprofile,sep='\t',index=False)
    return dfRad

 

def spanwiseED(tsAvg,vr,R,postprofile=None, IR=None):
    nr=len(vr)
    # --- Extract radial data
    Columns=[]
    for sB in ['b1','b2','b3']:
        SB=sB.upper()
        Columns.append(extractSpanTSReg(tsAvg,'^Spn(\d)ALx'+sB+'_\[m/s^2\]',SB+'ALx_[m/s^2]'))
        Columns.append(extractSpanTSReg(tsAvg,'^Spn(\d)ALy'+sB+'_\[m/s^2\]',SB+'ALy_[m/s^2]'))
        Columns.append(extractSpanTSReg(tsAvg,'^Spn(\d)ALz'+sB+'_\[m/s^2\]',SB+'ALz_[m/s^2]'))
        Columns.append(extractSpanTSReg(tsAvg,'^Spn(\d)TDx'+sB+'_\[m\]'    ,SB+'TDx_[m]'))
        Columns.append(extractSpanTSReg(tsAvg,'^Spn(\d)TDy'+sB+'_\[m\]'    ,SB+'TDy_[m]'))
        Columns.append(extractSpanTSReg(tsAvg,'^Spn(\d)TDz'+sB+'_\[m\]'    ,SB+'TDz_[m]'))
        Columns.append(extractSpanTSReg(tsAvg,'^Spn(\d)RDx'+sB+'_\[deg\]'  ,SB+'RDx_[deg]'))
        Columns.append(extractSpanTSReg(tsAvg,'^Spn(\d)RDy'+sB+'_\[deg\]'  ,SB+'RDy_[deg]'))
        Columns.append(extractSpanTSReg(tsAvg,'^Spn(\d)RDz'+sB+'_\[deg\]'  ,SB+'RDz_[deg]'))
        Columns.append(extractSpanTSReg(tsAvg,'^Spn(\d)FLx'+sB+'_\[kN\]'   ,SB+'FLx_[kN]'))
        Columns.append(extractSpanTSReg(tsAvg,'^Spn(\d)FLy'+sB+'_\[kN\]'   ,SB+'FLy_[kN]'))
        Columns.append(extractSpanTSReg(tsAvg,'^Spn(\d)FLz'+sB+'_\[kN\]'   ,SB+'FLz_[kN]'))
        Columns.append(extractSpanTSReg(tsAvg,'^Spn(\d)MLy'+sB+'_\[kN-m\]' ,SB+'MLx_[kN-m]'  ))
        Columns.append(extractSpanTSReg(tsAvg,'^Spn(\d)MLx'+sB+'_\[kN-m\]' ,SB+'MLy_[kN-m]'   ))
        Columns.append(extractSpanTSReg(tsAvg,'^Spn(\d)MLz'+sB+'_\[kN-m\]' ,SB+'MLz_[kN-m]'   ))

    dfRad, nrMax, ValidRow = _HarmonizeSpanwiseData('ElastoDyn', Columns, vr, R, IR=IR)

    # --- Export to csv
    if postprofile is not None and dfRad is not None:
        dfRad.to_csv(postprofile,sep='\t',index=False)
    return dfRad

def spanwiseAD(tsAvg,vr=None,rho=None,R=None,nB=None,chord=None,postprofile=None,IR=None):
    # --- Extract radial data
    Columns=[]
    for sB in ['B1','B2','B3']:
        Columns.append(extractSpanTSReg(tsAvg,'^'+sB+'N(\d*)Alpha_\[deg\]',sB+'Alpha_[deg]'))
        Columns.append(extractSpanTSReg(tsAvg,'^'+sB+'N(\d*)AOA_\[deg\]'  ,sB+'Alpha_[deg]')) # DBGOuts
        Columns.append(extractSpanTSReg(tsAvg,'^'+sB+'N(\d*)AxInd_\[-\]'  ,sB+'AxInd_[-]'  ))
        Columns.append(extractSpanTSReg(tsAvg,'^'+sB+'N(\d*)TnInd_\[-\]'  ,sB+'TnInd_[-]'  ))
        Columns.append(extractSpanTSReg(tsAvg,'^'+sB+'N(\d*)AIn_\[deg\]'  ,sB+'AxInd_[-]'  )) # DBGOuts NOTE BUG
        Columns.append(extractSpanTSReg(tsAvg,'^'+sB+'N(\d*)ApI_\[deg\]'  ,sB+'TnInd_[-]'  )) # DBGOuts NOTE BUG
        Columns.append(extractSpanTSReg(tsAvg,'^'+sB+'N(\d*)AIn_\[-\]'    ,sB+'AxInd_[-]'  )) # DBGOuts
        Columns.append(extractSpanTSReg(tsAvg,'^'+sB+'N(\d*)ApI_\[-\]'    ,sB+'TnInd_[-]'  )) # DBGOuts
        Columns.append(extractSpanTSReg(tsAvg,'^'+sB+'N(\d*)Cl_\[-\]'     ,sB+'Cl_[-]'   ))
        Columns.append(extractSpanTSReg(tsAvg,'^'+sB+'N(\d*)Cd_\[-\]'     ,sB+'Cd_[-]'   ))
        Columns.append(extractSpanTSReg(tsAvg,'^'+sB+'N(\d*)Cm_\[-\]'     ,sB+'Cm_[-]'   ))
        Columns.append(extractSpanTSReg(tsAvg,'^'+sB+'N(\d*)Cx_\[-\]'     ,sB+'Cx_[-]'   ))
        Columns.append(extractSpanTSReg(tsAvg,'^'+sB+'N(\d*)Cy_\[-\]'     ,sB+'Cy_[-]'   ))
        Columns.append(extractSpanTSReg(tsAvg,'^'+sB+'N(\d*)Cn_\[-\]'     ,sB+'Cn_[-]'   ))
        Columns.append(extractSpanTSReg(tsAvg,'^'+sB+'N(\d*)Ct_\[-\]'     ,sB+'Ct_[-]'   ))
        Columns.append(extractSpanTSReg(tsAvg,'^'+sB+'N(\d*)Re_\[-\]'     ,sB+'Re_[-]' ))
        Columns.append(extractSpanTSReg(tsAvg,'^'+sB+'N(\d*)Vrel_\[m/s\]' ,sB+'Vrel_[m/s]' ))
        Columns.append(extractSpanTSReg(tsAvg,'^'+sB+'N(\d*)Theta_\[deg\]',sB+'Theta_[deg]'))
        Columns.append(extractSpanTSReg(tsAvg,'^'+sB+'N(\d*)Phi_\[deg\]'  ,sB+'Phi_[deg]'))
        Columns.append(extractSpanTSReg(tsAvg,'^'+sB+'N(\d*)Twst_\[deg\]' ,sB+'Twst_[deg]')) #DBGOuts
        Columns.append(extractSpanTSReg(tsAvg,'^'+sB+'N(\d*)Curve_\[deg\]',sB+'Curve_[deg]'))
        Columns.append(extractSpanTSReg(tsAvg,'^'+sB+'N(\d*)Vindx_\[m/s\]',sB+'Vindx_[m/s]'))
        Columns.append(extractSpanTSReg(tsAvg,'^'+sB+'N(\d*)Vindy_\[m/s\]',sB+'Vindy_[m/s]'))
        Columns.append(extractSpanTSReg(tsAvg,'^'+sB+'N(\d*)Fx_\[N/m\]'   ,sB+'Fx_[N/m]'   ))
        Columns.append(extractSpanTSReg(tsAvg,'^'+sB+'N(\d*)Fy_\[N/m\]'   ,sB+'Fy_[N/m]'   ))
        Columns.append(extractSpanTSReg(tsAvg,'^'+sB+'N(\d*)Fl_\[N/m\]'   ,sB+'Fl_[N/m]'   ))
        Columns.append(extractSpanTSReg(tsAvg,'^'+sB+'N(\d*)Fd_\[N/m\]'   ,sB+'Fd_[N/m]'   ))
        Columns.append(extractSpanTSReg(tsAvg,'^'+sB+'N(\d*)Fn_\[N/m\]'   ,sB+'Fn_[N/m]'   ))
        Columns.append(extractSpanTSReg(tsAvg,'^'+sB+'N(\d*)Ft_\[N/m\]'   ,sB+'Ft_[N/m]'   ))
        Columns.append(extractSpanTSReg(tsAvg,'^'+sB+'N(\d*)VUndx_\[m/s\]',sB+'VUndx_[m/s]'))
        Columns.append(extractSpanTSReg(tsAvg,'^'+sB+'N(\d*)VUndy_\[m/s\]',sB+'VUndy_[m/s]'))
        Columns.append(extractSpanTSReg(tsAvg,'^'+sB+'N(\d*)VUndz_\[m/s\]',sB+'VUndz_[m/s]'))
        Columns.append(extractSpanTSReg(tsAvg,'^'+sB+'N(\d*)VDisx_\[m/s\]',sB+'VDisx_[m/s]'))
        Columns.append(extractSpanTSReg(tsAvg,'^'+sB+'N(\d*)VDisy_\[m/s\]',sB+'VDisy_[m/s]'))
        Columns.append(extractSpanTSReg(tsAvg,'^'+sB+'N(\d*)VDisz_\[m/s\]',sB+'VDisz_[m/s]'))
        Columns.append(extractSpanTSReg(tsAvg,'^'+sB+'N(\d*)Vx_\[m/s\]'   ,sB+'Vx_[m/s]'))
        Columns.append(extractSpanTSReg(tsAvg,'^'+sB+'N(\d*)Vy_\[m/s\]'   ,sB+'Vy_[m/s]'))
        Columns.append(extractSpanTSReg(tsAvg,'^'+sB+'N(\d*)Vz_\[m/s\]'   ,sB+'Vz_[m/s]'))
        Columns.append(extractSpanTSReg(tsAvg,'^'+sB+'N(\d*)DynP_\[Pa\]'  ,sB+'DynP_[Pa]' ))
        Columns.append(extractSpanTSReg(tsAvg,'^'+sB+'N(\d*)M_\[-\]'      ,sB+'M_[-]' ))
        Columns.append(extractSpanTSReg(tsAvg,'^'+sB+'N(\d*)Mm_\[N-m/m\]' ,sB+'Mm_[N-m/m]'   ))
        Columns.append(extractSpanTSReg(tsAvg,'^'+sB+'N(\d*)Gam_\['       ,sB+'Gam_[m^2/s]')) #DBGOuts

    # --- AD 14
    Columns.append(extractSpanTSReg(tsAvg,'^Alpha(\d*)_\[deg\]'    ,'Alpha_[deg]'  ,  IR=IR))
    Columns.append(extractSpanTSReg(tsAvg,'^DynPres(\d*)_\[Pa\]'   ,'DynPres_[Pa]' ,  IR=IR))
    Columns.append(extractSpanTSReg(tsAvg,'^CLift(\d*)_\[-\]'      ,'CLift_[-]'    ,  IR=IR))
    Columns.append(extractSpanTSReg(tsAvg,'^CDrag(\d*)_\[-\]'      ,'CDrag_[-]'    ,  IR=IR))
    Columns.append(extractSpanTSReg(tsAvg,'^CNorm(\d*)_\[-\]'      ,'CNorm_[-]'    ,  IR=IR))
    Columns.append(extractSpanTSReg(tsAvg,'^CTang(\d*)_\[-\]'      ,'CTang_[-]'    ,  IR=IR))
    Columns.append(extractSpanTSReg(tsAvg,'^CMomt(\d*)_\[-\]'      ,'CMomt_[-]'    ,  IR=IR))
    Columns.append(extractSpanTSReg(tsAvg,'^Pitch(\d*)_\[deg\]'    ,'Pitch_[deg]'  ,  IR=IR))
    Columns.append(extractSpanTSReg(tsAvg,'^AxInd(\d*)_\[-\]'      ,'AxInd_[-]'    ,  IR=IR))
    Columns.append(extractSpanTSReg(tsAvg,'^TanInd(\d*)_\[-\]'     ,'TanInd_[-]'   ,  IR=IR))
    Columns.append(extractSpanTSReg(tsAvg,'^ForcN(\d*)_\[N\]'      ,'ForcN_[N]'    ,  IR=IR))
    Columns.append(extractSpanTSReg(tsAvg,'^ForcT(\d*)_\[N\]'      ,'ForcT_[N]'    ,  IR=IR))
    Columns.append(extractSpanTSReg(tsAvg,'^Pmomt(\d*)_\[N-m\]'    ,'Pmomt_[N-N]'  ,  IR=IR))
    Columns.append(extractSpanTSReg(tsAvg,'^ReNum(\d*)_\[x10^6\]'  ,'ReNum_[x10^6]',  IR=IR))
    Columns.append(extractSpanTSReg(tsAvg,'^Gamma(\d*)_\[m^2/s\]'  ,'Gamma_[m^2/s]',  IR=IR))

    dfRad, nrMax, ValidRow = _HarmonizeSpanwiseData('AeroDyn', Columns, vr, R, IR=IR)
    
    # --- Compute additional values (AD15 only)
    if chord is not None:
        if vr is not None:
            chord =chord[0:nrMax]
        chord    = chord   [ValidRow[0]]
    for sB in ['B1','B2','B3']:
        try:
            vr_bar=vr/R
            Fx = dfRad[sB+'Fx_[N/m]']
            U0 = tsAvg['Wind1VelX_[m/s]']
            Ct=nB*Fx/(0.5 * rho * 2 * U0**2 * np.pi * r)
            Ct[vr<0.01*R] = 0
            dfRad[sB+'Ct_[-]'] = Ct
            CT=2*np.trapz(vr_bar*Ct,vr_bar)
            dfRad[sB+'CtAvg_[-]']= CT*np.ones(r.shape)
        except:
            pass
        try:
            dfRad[sB+'Gamma_[m^2/s]'] = 1/2 * chord*  dfRad[sB+'Vrel_[m/s]'] * dfRad[sB+'Cl_[-]'] 
        except:
            pass
        try: 
            if not sB+'Vindx_[m/s]' in dfRad.columns:
                dfRad[sB+'Vindx_[m/s]']= -dfRad[sB+'AxInd_[-]'].values * dfRad[sB+'Vx_[m/s]'].values 
                dfRad[sB+'Vindy_[m/s]']=  dfRad[sB+'TnInd_[-]'].values * dfRad[sB+'Vy_[m/s]'].values 
        except:
            pass

    # --- Export to csv
    if postprofile is not None and dfRad is not None:
        dfRad.to_csv(postprofile,sep='\t',index=False)
    return dfRad



def spanwisePostProLegacy(FST_In=None,avgMethod='constantwindow',avgParam=5,out_ext='.outb',postprofile=None,df=None):
    """
    Postprocess FAST radial data

    INPUTS:
        - FST_IN: Fast .fst input file
        - avgMethod='periods', avgParam=2:  average over 2 last periods, Needs Azimuth sensors!!!
        - avgMethod='constantwindow', avgParam=5:  average over 5s of simulation
        - postprofile: outputfile to write radial data
    """
    # --- Opens Fast output  and performs averaging
    if df is None:
        df = weio.read(FST_In.replace('.fst',out_ext)).toDataFrame()
        returnDF=True
    else:
        returnDF=False
    # NOTE: spanwise script doest not support duplicate columns
    df = df.loc[:,~df.columns.duplicated()]
    dfAvg = averageDF(df,avgMethod=avgMethod ,avgParam=avgParam) # NOTE: average 5 last seconds

    # --- Extract info (e.g. radial positions) from Fast input file
    # We don't have a .fst input file, so we'll rely on some default values for "r"
    rho         = 1.225
    chord       = None
    # --- Extract radial positions of output channels
    r_AD, r_ED, r_BD, IR_AD, IR_ED, IR_BD, R, r_hub, fst = FASTRadialOutputs(FST_In, OutputCols=df.columns.values)
    if R is None: 
        R=1
    try:
        chord  = fst.AD.Bld1['BldAeroNodes'][:,5] # Full span
    except:
        pass
    try:
        rho = fst.AD['Rho']
    except:
        rho = fst.AD['AirDens']
    print('r_AD:', r_AD)
    print('r_ED:', r_ED)
    print('r_BD:', r_BD)
    #print('I_AD:', IR_AD)
    #print('I_ED:', IR_ED)
    #print('I_BD:', IR_BD)
    # --- Extract radial data and export to csv if needed
    dfRad_AD    = None
    dfRad_ED    = None
    dfRad_BD    = None

    dfRad_AD   = spanwiseAD(dfAvg.iloc[0],  r_AD, rho , R=R, nB=3, chord=chord, postprofile=postprofile, IR=IR_AD)
    if r_ED is not None:
        dfRad_ED = spanwiseED(dfAvg.iloc[0], r_ED, R=R, IR=IR_ED, postprofile=postprofile)
    if r_BD is not None:
        dfRad_BD = spanwiseBD(dfAvg.iloc[0], r_BD, R=R, IR=IR_BD, postprofile=postprofile)

    if returnDF:
        return dfRad_ED , dfRad_AD, dfRad_BD, df
    else:
        return dfRad_ED , dfRad_AD, dfRad_BD
