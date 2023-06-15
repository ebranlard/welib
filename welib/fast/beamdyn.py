import numpy as np
from numpy import cos, sin
import pandas as pd
import os
import welib.weio as weio
from welib.weio.fast_input_file import FASTInputFile

# --------------------------------------------------------------------------------}
# --- Writer function 
# --------------------------------------------------------------------------------{
def write_beamdyn_sections(filename,span,lK,lM,Mu=None,Label=''):
    """ Write a BeamDyn section file, 
    span : list of nSpan span values from 0 to 1
    lK   : list of nSpan 6*6 stiffness matrices
    lM   : list of nSpan 6*6 mass matrices
    Mu   : damping coefficient
    """
    if (Mu is None):
        Mu=[0]*6
        damp_type=0
    elif not hasattr(Mu, '__len__'):
        Mu=np.asarray([Mu]*6)
    Mu = np.asarray(Mu)
    if len(Mu)==6:
        damp_type=1

    # --- Helper functions
    def mat_tostring(M,fmt='.5e'):
        return '\n'.join(['   '+' '.join(['{:.6E}'.format(m) for m in M[i,:]]) for i in range(np.size(M,1))])
    def beamdyn_section_mat_tostring(x,K,M):
        s=''
        s+='{:.6f}\n'.format(x)
        s+=mat_tostring(K)
        #s+=np.array2string(K)
        s+='\n'
        s+='\n'
        s+=mat_tostring(M)
        #s+=np.array2string(M)
        s+='\n'
        s+='\n'
        return s

    # --- Writing
    with open(filename, 'w') as f:
        f.write('------- BEAMDYN V1.00.* INDIVIDUAL BLADE INPUT FILE --------------------------\n')
        f.write('! {} - Written using {} \n'.format(Label,os.path.basename(__file__)))
        f.write('---------------------- BLADE PARAMETERS --------------------------------------\n')
        f.write('{:5d}  station_total    - Number of blade input stations (-)\n'.format(len(span)))
        f.write('{:5d}  damp_type        - Damping type (switch): 0: no damping; 1: viscous damping\n'.format(damp_type))
        f.write('---------------------- DAMPING COEFFICIENT------------------------------------\n')
        f.write('   mu1        mu2        mu3        mu4        mu5        mu6\n')
        f.write('   (s)        (s)        (s)        (s)        (s)        (s)\n')
        f.write(' {} {} {} {} {} {} \n'.format(*[Mu[i] for i in range(6)]))
        f.write('---------------------- DISTRIBUTED PROPERTIES---------------------------------\n')
        for s,K,M in zip(span,lK,lM):
            f.write(beamdyn_section_mat_tostring(s,K,M))

# --------------------------------------------------------------------------------}
# --- Hawc2 to BeamDyn 
# --------------------------------------------------------------------------------{
def mypolyfit(x,y,exponents=[0,1,2,3]):
    X_poly=np.array([])
    for i,c in enumerate(exponents):
        if i==0:
            X_poly = x**c
        else:
            X_poly = np.vstack((X_poly,x**c))
    try:
        coeffs = np.linalg.lstsq(X_poly.T, y, rcond=None)[0]
    except:
        coeffs = np.linalg.lstsq(X_poly.T, y)
    #print('Poly fit coeffs: ' + '+'.join(['{:.5f}^{}'.format(p,c) for p,c in zip(coeffs,exponents)]))
    return np.dot(coeffs, X_poly)

def htcToBeamDyn(HTCFile, bodyname, BDBldFileOut, BDMainFileOut=None, BDMainTemplate=None, Mu = 1.0e-03, 
                   ref_axis='c2def-polyfit', poly_exp=[2,3,4,5], zref=None, Label='', bPlot=False,
                   bNoOffset=False, bNoPreSweep=False, nRootZeros=0, bCGOnMeanLine=False):  # Experimental options
    """
    Writes BeamDyn inputs files from a HAWC2 htc file and the blade body name
    INPUTS:
      - HTCFile: path to a htc file
      - bodyname
    OTHER INPUTS:
       see hawc2tobeamdyn
    """
    htc = weio.hawc2_htc_file.HAWC2HTCFile(HTCFile)
    dfs = htc.toDataFrame()
    H2MeanLine = dfs[bodyname+'_c2']
    H2Structure = dfs[bodyname+'_st']
    H2MeanLine = H2MeanLine[['x_[m]','y_[m]','z_[m]','twist_[deg]']] # Changing order

    return hawc2ToBeamDyn(H2MeanLine, H2Structure, BDBldFileOut, BDMainFileOut=BDMainFileOut, BDMainTemplate=BDMainTemplate, Mu=Mu, 
                   ref_axis=ref_axis, poly_exp=poly_exp, zref=zref, Label=Label, bPlot=bPlot,
                   bNoOffset=bNoOffset, bNoPreSweep=bNoPreSweep, nRootZeros=nRootZeros, bCGOnMeanLine=bCGOnMeanLine)

def hawc2ToBeamDyn(H2MeanLine, H2Structure, BDBldFileOut, BDMainFileOut=None, BDMainTemplate=None, Mu = 1.0e-03, 
                   ref_axis='c2def-polyfit', poly_exp=[2,3,4,5], zref=None, Label='', bPlot=False,
                   bNoOffset=False, bNoPreSweep=False, nRootZeros=0, bCGOnMeanLine=False):  # Experimental options
    """
    Writes BeamDyn inputs files from two csv files derived from "Hawc2" inputs

    INPUTS:
        - H2MeanLine :  dataframe or csv file with one header line, containing c2 def definition, (Hawc2 coordinates)
                           Column order has to be:  ['x_[m]','y_[m]','z_[m]','twist_[deg]']
        - H2Structure: dataframe or csv that contains Hawc2 beam structural properties (typically found in a hawc2 st file)
                           Colums order has to be:  ['r_[m]','m_[kg/m]','x_cg_[m]','y_cg_[m]','ri_x_[m]',... ,'pitch_[deg]','x_e_[m]','y_e_[m]']
        - BDBldFileOut:  filepath of the BeamDyn Blade file to be written

    OPTIONAL INPUTS:
        - BDMainFileOut:  filepath of the BeamDyn main file to be written 
                          The file will contain the mean line definition. Requires BDMainTemplate
        - BDMainTemplate: filepath of a BeamDyn main file, will be used as template for BDMainFileOut
        - Mu : list of the 6 damping coeffficients to be used in the blade file
        - ref_axis: string defining how the main axis of beamdyn will be defined.
                    'c2def-polyfit': the reference axis is Hawc2 c2def, smoothened out using a polyfit (see poly_exp)
                    'straight': the reference axis is straight (prebend and sweep are still included as offsets) 
        - poly_exp: list of exponents used to perform the polyfit of Hawc2 c2def line
        - zref: specifies "z" locations where sections have to be interpolated to. If "None", hawc2 mean line sections are used
        - Label : string used as a label for the blade file
        - bPlot : boolean, if true, a plot is generated

    EXPERIMENTAL INPUTS (not recommended, keep as default):
        - bNoOffset: do not use offsets from mean line 
                    If used with ref_axis='straight', this results in fully straight blade).
        - bCBOnMeanLine: assumes CB on the mean axis
        - bNoPreSweep: remove presweep
        - nRootZeros: number of points that are set to have zero x and y at the root
    """
    import welib.weio.csv_file
    # --- Mean line definition (position and orientation)  - Hawc2 "c2-def file", BeanDyn point "O"
    if isinstance(H2MeanLine, pd.DataFrame):
        c2def = H2MeanLine
    else:
        c2def = weio.csv_file.CSVFile(H2MeanLine).toDataFrame()
    # For the equations below to be generic we force the column names
    c2def.columns.values[0:4]=['x_[m]','y_[m]','z_[m]','twist_[deg]']

    # --- If necessary, interpolate mean line to user defined positions
    c2def_old = c2def.copy()
    if zref is None:
        # we dont interpolate
        zref=c2def['z_[m]'].values
    else:
        z_old = c2def_old['z_[m]'].values
        # safety checks
        zref=np.asarray(zref)
        if z_old[0]!=zref[0]:
            raise Exception('`zref` start value should be {} to match input'.format(z_old[0]))
        if z_old[-1]!=zref[-1]:
            raise Exception('`zref` end value should be {} to match input'.format(z_old[-1]))
        # interpolating to zref values
        c2def     = c2def[0:0] # emptying
        for c in c2def_old.columns:
            c2def[c] = np.interp(zref, z_old, c2def_old[c])

    # --- Hawc2 ref axis (in BeamDyn system)
    x_O_h2 =   c2def['y_[m]'].values  # kp_xr  
    y_O_h2 = - c2def['x_[m]'].values  # kp_yr  
    z_O_h2 =   c2def['z_[m]'].values  # kp_zr
    twist  = - c2def['twist_[deg]'].values # initial_twist [deg] Hawc2 angle is positive around z, unlike the "twist"

    # --- Compute r_ref, curvilinear position of keypoints (for st file)
    dr= np.sqrt((x_O_h2[1:]-x_O_h2[0:-1])**2 +(y_O_h2[1:]-y_O_h2[0:-1])**2 +(z_O_h2[1:]-z_O_h2[0:-1])**2)
    r_ref = np.concatenate(([0],np.cumsum(dr)))


    # --- BeamDyn ref axis
    # Default: taken as c2def 
    x_O = x_O_h2.copy() # kp_xr
    y_O = y_O_h2.copy() # kp_yr
    z_O = z_O_h2.copy() # kp_zr
    x_O[:nRootZeros] = 0
    y_O[:nRootZeros] = 0
    if ref_axis=='c2def':
        # (see above)
        pass
    elif ref_axis=='straight':
        # straight axis, with everything as offsets
        x_O = 0*x_O_h2   # kp_xr
        y_O = 0*y_O_h2   # kp_yr
    elif ref_axis=='y-straight-polyfit':
        # y-axis straight, x-axis poly fitted,  with everything as offsets
        y_O = 0*y_O_h2   # kp_yr
        x_O = mypolyfit(z_O_h2, x_O, poly_exp)  # kp_xr NOTE: we fit x_O (where nRoot was already inforced
        x_O[:nRootZeros] =0 # enforcing zero displacements at root
    elif ref_axis=='c2def-polyfit':
        # Smooth mean line definition (position and orientation)
        x_O = mypolyfit(z_O_h2, x_O, poly_exp)  # kp_xr NOTE: we fit x_O (where nRoot was already inforced
        y_O = mypolyfit(z_O_h2, y_O, poly_exp)  # kp_yr  
        x_O[:nRootZeros] =0 # enforcing zero displacements at root
        y_O[:nRootZeros] =0
    else:
        raise NotImplementedError('ref_axis: {}'.format(ref_axis))

    # difference between input_axis and smooth axis in global (blade root, BeamDyn convention)
    x_off_g = (x_O_h2-x_O)
    y_off_g = (y_O_h2-y_O)

    if bNoOffset:
        x_off_g = 0*x_off_g  # no offsets
        y_off_g = 0*y_off_g

    # transform offset from global to local axis orientation
    theta_z = -twist
    x_off_s =   x_off_g * np.cos(theta_z) + y_off_g * np.sin(theta_z)
    y_off_s =  -x_off_g * np.sin(theta_z) + y_off_g * np.cos(theta_z)




    # --- Cross section properties 
    if isinstance(H2Structure, pd.DataFrame):
        hwc_in = H2Structure
    else:
        hwc_in = weio.csv_file.CSVFile(H2Structure).toDataFrame()
    # For the equations below to be generic we force the column names
    hwc_in.columns=['r_[m]','m_[kg/m]','x_cg_[m]','y_cg_[m]','ri_x_[m]','ri_y_[m]','x_sh_[m]','y_sh_[m]','E_[N/m^2]','G_[N/m^2]','I_x_[m^4]','I_y_[m^4]','I_p_[m^4]','k_x_[-]','k_y_[-]','A_[m^2]','pitch_[deg]','x_e_[m]','y_e_[m]']
    # --- Interpolating to match c2def positions
    hwc = pd.DataFrame(columns=hwc_in.columns)
    r_old = hwc_in['r_[m]'].values
    for c in hwc.columns:
        hwc[c] = np.interp(r_ref, r_old, hwc_in[c])
    if r_old[-1]<r_ref[-1]:
        # NOTE: interp won't do extrap , small hack here...
        hwc['r_[m]'].values[-1] = r_ref[-1]

    # --- Safety check
    if len(hwc)!=len(c2def):
        raise Exception('Interpolation failed, wrong length. Debug me.')
    if any(np.abs(hwc['r_[m]'].values-r_ref)>1e-9):
        raise Exception('Interpolation failed, radial position mismatch. Debug me.')

    # --- Setting Mass and stiffness matrices of each cross section
    lM=[]; lK=[]
    vx_G=[]; vy_G=[];
    vx_S=[]; vy_S=[];
    vx_C=[]; vy_C=[];
    for i,row in hwc.iterrows():
        if i<nRootZeros:
            x_G =0
            y_G =0
            x_S =0
            y_S =0
            x_C =0
            y_C =0
        else:
            if bCGOnMeanLine:
                x_G     = 0
                y_G     = 0
            else:
                x_G     =  row['y_cg_[m]']  +x_off_s[i]
                y_G     = -row['x_cg_[m]']  +y_off_s[i]
            x_S     =  row['y_sh_[m]']  +x_off_s[i]
            y_S     = -row['x_sh_[m]']  +y_off_s[i]
            x_C     =  row['y_e_[m]']   +x_off_s[i]
            y_C     = -row['x_e_[m]']   +y_off_s[i]

        EA      = row['E_[N/m^2]']* row['A_[m^2]']
        GKt     = row['G_[N/m^2]']* row['I_p_[m^4]']
        GA      = row['G_[N/m^2]']* row['A_[m^2]']
        kxs     = row['k_y_[-]']
        kys     = row['k_x_[-]']
        EIxp    = row['E_[N/m^2]']* row['I_y_[m^4]']   # Should be [N.m^2]
        EIyp    = row['E_[N/m^2]']* row['I_x_[m^4]']
        theta_s = row['pitch_[deg]']*np.pi/180 # [rad]
        theta_p = row['pitch_[deg]']*np.pi/180 # NOTE: hawc2 and our convention here for angles are positive around z (whereas beamdyn take them negative around z)
        theta_i = row['pitch_[deg]']*np.pi/180 # [rad]

        m       = row['m_[kg/m]']
        Ixi     = row['ri_y_[m]']**2 * m    # [kg.m]              
        Iyi     = row['ri_x_[m]']**2 * m    # [kg.m]
        I_p     = Ixi+Iyi                   # [kg.m]

        M =  MM(m, Ixi, Iyi, I_p, x_G, y_G, theta_i) # NOTE: theta_i in rad
        K =  KK(EA, EIxp, EIyp, GKt, GA, kxs, kys, x_C, y_C, theta_p, x_S, y_S, theta_s) # Note theta_p/s in rad
        
        lM.append(M)
        lK.append(K)
        vx_G.append(x_G); vy_G.append(y_G)
        vx_S.append(x_S); vy_S.append(y_S)
        vx_C.append(x_C); vy_C.append(y_C)

    # --- Writing BeamDyn blade file
    span=hwc['r_[m]'].values
    s_bar=span/span[-1]
    print('Writing BeamDyn blade file:',BDBldFileOut)
    write_beamdyn_sections(BDBldFileOut,s_bar,lK,lM,Mu,Label=Label)

    # --- db
    #M=np.column_stack((zref, x_off, y_off))
    #np.savetxt(BDBldFileOut.replace('.dat','offsets.txt'), M, delimiter=',',header='z_[m], xoff_[m], yoff_[m]')

    # --- Writing BeamDyn main file based on template file
    if BDMainTemplate is not None and BDMainFileOut is not None:
        BD=FASTInputFile(BDMainTemplate)
        #print(BD.keys())
        BD.comment=Label
        BD['MemberGeom'] = np.column_stack((x_O,y_O,z_O,twist))
        BD['kp_total']   = len(x_O)
        BD['BldFile']    = '"'+os.path.basename(BDBldFileOut)+'"' 
        # TODO TODO 
        BD.data[BD.getID('kp_total')+1]['value']= '1 {}'.format(len(x_O))

        print('Writing BeamDyn file:',BDMainFileOut)
        BD.write(BDMainFileOut)


    # ---
    if bPlot:
        import matplotlib.pyplot as plt
        colrs=plt.rcParams['axes.prop_cycle'].by_key()['color']

        EdgStiff= np.array([K[3,3] for K in lK])
        FlpStiff= np.array([K[4,4] for K in lK])

        EIxp    = hwc['E_[N/m^2]']*hwc['I_y_[m^4]'].values   # Should be [N.m^2]
        EIyp    = hwc['E_[N/m^2]']*hwc['I_x_[m^4]'].values

#         fig=plt.figure()
        fig,axes = plt.subplots(4, 2, sharex=True, figsize=(12.4,09.)) # (6.4,4.8)
        fig.subplots_adjust(left=0.07, right=0.99, top=0.98, bottom=0.07, hspace=0.25, wspace=0.15)
        for ax in axes.ravel():
            ax.tick_params(direction='in')

        # --- Plot mean line from hawc2 and beamdyn
        x_O_h2 =   c2def_old['y_[m]'].values   # kp_xr  
        y_O_h2 =  -c2def_old['x_[m]'].values   # kp_yr  
        z_O_h2 =   c2def_old['z_[m]'].values   # kp_zr
        twist  =  -c2def_old['twist_[deg]'].values # initial_twist [deg]
#         fig,axes = plt.subplots(2, 1, sharex=True, figsize=(6.4,4.8)) # (6.4,4.8)
#         ax=axes[0,0]

        axes[0,0].text(0.5, 1.01, 'Mean line x', horizontalalignment='center', verticalalignment='bottom', transform = axes[0,0].transAxes)
        axes[0,1].text(0.5, 1.01, 'Mean line y', horizontalalignment='center', verticalalignment='bottom', transform = axes[0,1].transAxes)
        axes[0,0].plot(z_O, x_O   , '-' , label = 'BD smooth)')
        axes[0,0].plot(z_O, x_O_h2, '--', label = 'H2 c2def', ms=3, color='k')
        axes[0,0].plot(z_O, x_off_g, ':', label = r'"Delta" to c2def', color=colrs[6])
        axes[0,1].plot(z_O, y_O   , '-' , label = 'BD y (smooth)')
        axes[0,1].plot(z_O, y_O_h2, '--' , label = 'H2 "y"', ms=3, color='k')
        axes[0,1].plot(z_O, y_off_g , ':', label = 'y_off', color=colrs[6])
        if 'Relative_thickness_[%]' and 'Chord_[m]' in c2def.columns.values:
            c = c2def['Chord_[m]']
            t = c2def['Relative_thickness_[%]'] *c/100
            axes[0,0].plot(z_O, x_O_h2+c/2*np.sin(twist*np.pi/180), '-', color=[0.5,0.5,0.5] )
            axes[0,0].plot(z_O, x_O_h2-c/2*np.sin(twist*np.pi/180), '-', color=[0.5,0.5,0.5] )
            axes[0,1].plot(z_O, y_O_h2+c/2*np.cos(twist*np.pi/180), '-', color=[0.5,0.5,0.5] )
            axes[0,1].plot(z_O, y_O_h2-c/2*np.cos(twist*np.pi/180), '-', color=[0.5,0.5,0.5] )

#	Chord_[m]	Relative_thickness_[%]
#         axes[0,0].set_xlabel('z [m]')
#         axes[0,1].set_xlabel('z [m]')
        axes[0,0].set_ylabel('x [m]')
        axes[0,1].set_ylabel('y [m]')
        axes[0,0].legend(loc='upper right', fontsize=8)

        # --- Plot COG, Shear Center
        vx_G=np.asarray(vx_G); vy_G=np.asarray(vy_G)
        vx_S=np.asarray(vx_S); vy_S=np.asarray(vy_S)
        vx_C=np.asarray(vx_C); vy_C=np.asarray(vy_C)
#         fig,axes = plt.subplots(2, 1, sharex=True, figsize=(6.4,4.8)) # (6.4,4.8)
#         fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        axes[1,0].text(0.5, 1.01, 'Abs. position, x', horizontalalignment='center', verticalalignment='bottom', transform = axes[1,0].transAxes)
        axes[1,1].text(0.5, 1.01, 'Abs. position, y', horizontalalignment='center', verticalalignment='bottom', transform = axes[1,1].transAxes)
        axes[1,0].plot(z_O              , x_O     , '-' , label = 'BD meanline')
        axes[1,0].plot(z_O              , x_O    + vx_G           , 'd' , ms=6, color=None, markeredgecolor=colrs[1], markerfacecolor="None", label = 'G (COG)') 
        axes[1,0].plot(z_O              , x_O    + vx_S           , 's' , ms=6, color=None, markeredgecolor=colrs[2], markerfacecolor="None", label = 'S (shear center)') 
        axes[1,0].plot(z_O              , x_O    + vx_C           , 'o' , ms=6, color=None, markeredgecolor=colrs[3], markerfacecolor="None", label = 'C (elastic center)') 
        axes[1,0].plot(hwc['r_[m]'].values, x_O_h2 + hwc['y_cg_[m]'].values, 'd' , ms=1, color=colrs[1] , label='HAWC2')
        axes[1,0].plot(hwc['r_[m]'].values, x_O_h2 + hwc['y_sh_[m]'].values, 's' , ms=1, color=colrs[2] )
        axes[1,0].plot(hwc['r_[m]'].values, x_O_h2 + hwc['y_e_[m]' ].values , 'o' , ms=1, color=colrs[3] )
        axes[1,1].plot(z_O, y_O     , '-' , label = 'BD y (smooth)')
        axes[1,1].plot(z_O              , y_O    + vy_G           , 'd' , ms=6, color=None, markeredgecolor=colrs[1], markerfacecolor="None", label = 'G (COG)') 
        axes[1,1].plot(z_O              , y_O    + vy_S           , 's' , ms=6, color=None, markeredgecolor=colrs[2], markerfacecolor="None", label = 'S (shear center)') 
        axes[1,1].plot(z_O              , y_O    + vy_C           , 'o' , ms=6, color=None, markeredgecolor=colrs[3], markerfacecolor="None", label = 'C (elastic center)') 
        axes[1,1].plot(hwc['r_[m]'].values, y_O_h2 - hwc['x_cg_[m]'].values, 'd' , ms=1, color=colrs[1] )
        axes[1,1].plot(hwc['r_[m]'].values, y_O_h2 - hwc['x_sh_[m]'].values, 's' , ms=1, color=colrs[2] )
        axes[1,1].plot(hwc['r_[m]'].values, y_O_h2 - hwc['x_e_[m]'].values, 'o' , ms=1, color=colrs[3] )
#         axes[1,0].set_xlabel('z [m]')
#         axes[1,1].set_xlabel('z [m]')
        axes[1,0].set_ylabel('x [m]')
        axes[1,1].set_ylabel('y [m]')
        axes[1,0].legend(loc='upper right', fontsize=8)
        if 'Relative_thickness_[%]' and 'Chord_[m]' in c2def.columns.values:
            c = c2def['Chord_[m]'].values
            t = c2def['Relative_thickness_[%]'].values *c/100
            axes[1,0].plot(z_O, x_O_h2+c/2*np.sin(twist*np.pi/180), '-', color=[0.5,0.5,0.5] )
            axes[1,0].plot(z_O, x_O_h2-c/2*np.sin(twist*np.pi/180), '-', color=[0.5,0.5,0.5] )
            axes[1,1].plot(z_O, y_O_h2+c/2*np.cos(twist*np.pi/180), '-', color=[0.5,0.5,0.5] )
            axes[1,1].plot(z_O, y_O_h2-c/2*np.cos(twist*np.pi/180), '-', color=[0.5,0.5,0.5] )
# 

        # --- Positions rel to mean line
        axes[2,0].text(0.5, 1.01, r'Pos. wrt. meanline, x', horizontalalignment='center', verticalalignment='bottom', transform = axes[2,0].transAxes)
        axes[2,1].text(0.5, 1.01, r'Pos. wrt. meanline, y', horizontalalignment='center', verticalalignment='bottom', transform = axes[2,1].transAxes)
        axes[2,0].plot(z_O              , x_O-x_O     , '-' , label = 'BD meanline')
        axes[2,0].plot(z_O              , vx_G        , 'd' , ms=6, color=None, markeredgecolor=colrs[1], markerfacecolor="None", label = 'G (COG)') 
        axes[2,0].plot(z_O              , vx_S        , 's' , ms=6, color=None, markeredgecolor=colrs[2], markerfacecolor="None", label = 'S (shear center)') 
        axes[2,0].plot(z_O              , vx_C        , 'o' , ms=6, color=None, markeredgecolor=colrs[3], markerfacecolor="None", label = 'C (elastic center)') 
        axes[2,0].plot(z_O              , x_off_g     , ':', label = r'"Delta" to c2def', color=colrs[6])
#         axes[1,0].plot(hwc['r_[m]'], x_O_h2 + hwc['y_cg_[m]'], 'o' , ms=1, color=colrs[1] , label='HAWC2')
#         axes[1,0].plot(hwc['r_[m]'], x_O_h2 + hwc['y_sh_[m]'], 'o' , ms=1, color=colrs[2] )
#         axes[1,0].plot(hwc['r_[m]'], x_O_h2 + hwc['y_e_[m]'] , 'd' , ms=1, color=colrs[3] )
        axes[2,1].plot(z_O              , y_O -y_O    , '-' , label = 'BD meanline')
        axes[2,1].plot(z_O              , vy_G        , 'd' , ms=6, color=None, markeredgecolor=colrs[1], markerfacecolor="None", label = 'G (COG)') 
        axes[2,1].plot(z_O              , vy_S        , 's' , ms=6, color=None, markeredgecolor=colrs[2], markerfacecolor="None", label = 'S (shear center)') 
        axes[2,1].plot(z_O              , vy_C        , 'o' , ms=6, color=None, markeredgecolor=colrs[3], markerfacecolor="None", label = 'C (elastic center)') 
        axes[2,1].plot(z_O              , y_off_g     , ':', label = r'"Delta" to c2def', color=colrs[6])
#         axes[1,1].plot(hwc['r_[m]'], y_O_h2 - hwc['x_cg_[m]'], 'o' , ms=1, color=colrs[1] )
#         axes[1,1].plot(hwc['r_[m]'], y_O_h2 - hwc['x_sh_[m]'], 's' , ms=1, color=colrs[2] )
#         axes[1,1].plot(hwc['r_[m]'], y_O_h2 - hwc['x_e_[m]'] , 'd' , ms=1, color=colrs[3] )
        axes[2,0].set_xlabel('z [m]')
        axes[2,1].set_xlabel('z [m]')
        axes[2,0].set_ylabel(r'Delta x [m]')
        axes[2,1].set_ylabel(r'Delta y [m]')
        axes[2,0].legend(loc='upper right', fontsize=8)
# 
        # --- Plot Stiffness
#         ax=fig.add_subplot(111)
        ax=axes[3,0]
        ax.text(0.5, 1.01, 'Stiffnesses', horizontalalignment='center', verticalalignment='bottom', transform = ax.transAxes)
        ax.plot(z_O,EdgStiff,'-' , color=colrs[0], label='Edge Stiffness (K_44)')
        ax.plot(z_O,EIxp.values    ,'--', color=colrs[0], label='EIx "edge" at elastic center')
        ax.plot(z_O,FlpStiff,'-' , color=colrs[1], label='Flap Stiffness (K_55)')
        ax.plot(z_O,EIyp.values    ,'--', color=colrs[1], label='EIy "flap" at elastic center')
        ax.set_xlabel('z [m]')
        ax.set_ylabel('Stiffness [Nm^2]')
        ax.legend(fontsize=8)

        return fig 
        #fig.savefig(BDMainFileOut.replace('.dat','.png'))
        #plt.show()



# --------------------------------------------------------------------------------}
# ---  
# --------------------------------------------------------------------------------{
def beamDynToHawc2(BD_mainfile, BD_bladefile, H2_htcfile=None, H2_stfile=None, bodyname=None, A=None, E=None, G=None, theta_p_in=None, FPM=False):
    """ 
    
     FPM: fully populated matrix, if True, use the FPM format of hawc2
    """
    # --- Read BeamDyn files
    if isinstance(BD_mainfile, str):
        BD_mainfile = weio.read(BD_mainfile)
    if isinstance(BD_bladefile, str):
        BD_bladefile = weio.read(BD_bladefile)
    bdLine = BD_mainfile.toDataFrame()
    bd     = BD_bladefile.toDataFrame()

    # --- Extract relevant info
    prop  = bd['BeamProperties']
    kp_x  = bdLine['kp_xr_[m]'].values
    kp_y  = bdLine['kp_yr_[m]'].values
    kp_z  = bdLine['kp_zr_[m]'].values
    twist = bdLine['initial_twist_[deg]'].values*np.pi/180 # BeamDyn convention
    r_bar = prop['Span'].values

    K = np.zeros((6,6),dtype='object')
    M = np.zeros((6,6),dtype='object')
    for i in np.arange(6):
        for j in np.arange(6):
            K[i,j]=prop['K{}{}'.format(i+1,j+1)].values
            M[i,j]=prop['M{}{}'.format(i+1,j+1)].values

    # Map 6x6 data to "beam" data
    # NOTE: theta_* are in [rad]
    EA, EIx, EIy, kxsGA, kysGA, GKt, x_C, y_C, x_S, y_S, theta_p, theta_s = K66toProps(K, theta_p_in)
    m, Ixi, Iyi, Ip, x_G, y_G, theta_i = M66toProps(M)
#     print('kxGA    {:e}'.format(np.mean(kxsGA)))
#     print('kyGA    {:e}'.format(np.mean(kysGA)))
#     print('EA      {:e}'.format(np.mean(EA)))
#     print('EIx     {:e}'.format(np.mean(EIx)))
#     print('EIy     {:e}'.format(np.mean(EIy)))
#     print('GKt     {:e}'.format(np.mean(GKt)))
#     print('xC    ',np.mean(x_C))
#     print('yC    ',np.mean(y_C))
#     print('xS    ',np.mean(x_S))
#     print('yS    ',np.mean(y_S))
#     print('thetap',np.mean(theta_p))
#     print('thetas',np.mean(theta_s))
#     print('m     ',np.mean(m))
#     print('Ixi   ',np.mean(Ixi))
#     print('Iyi   ',np.mean(Iyi))
#     print('Ip    ',np.mean(Ip))
#     print('x_G   ',np.mean(x_G))
#     print('y_G   ',np.mean(y_G))
#     print('thetai',np.mean(theta_i))

    # Convert to Hawc2 system
    if FPM:
        dfMeanLine , dfStructure = beamDyn2Hawc2FPM_raw(r_bar,
                kp_x, kp_y, kp_z, twist,  # BeamDyn convention, twist around -z [in rad]
                m, Ixi, Iyi, x_G, y_G, theta_i,  # theta_i/p around z (in rad)
                x_C, y_C, theta_p, K)

    else:

        dfMeanLine , dfStructure = beamDyn2Hawc2_raw(r_bar,
                kp_x, kp_y, kp_z, twist, 
                m, Ixi, Iyi, x_G, y_G, theta_i,
                EA, EIx, EIy, GKt, kxsGA, kysGA, x_C, y_C, theta_p, x_S, y_S, theta_s, 
                A=A, E=E, G=G)

    # --- Rewrite st file
    if H2_stfile is not None:
        with open(H2_stfile, 'w') as f:
            f.write('%i ; number of sets, Nset\n' % 1)
            f.write('-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n')
            f.write('#%i ; set number\n' % 1)
            if FPM:
                cols=['r','m_[kg/m]','x_cg_[m]','y_cg_[m]','ri_x_[m]','ri_y_[m]','pitch_[deg]','x_e_[m]','y_e_[m]','K11','K12','K13','K14','K15','K16','K22','K23','K24','K25','K26','K33','K34','K35','K36','K44','K45','K46','K55','K56','K66']
            else:
                cols=['r_[m]','m_[kg/m]','x_cg_[m]','y_cg_[m]','ri_x_[m]','ri_y_[m]', 'x_sh_[m]','y_sh_[m]','E_[N/m^2]','G_[N/m^2]','I_x_[m^4]','I_y_[m^4]','I_p_[m^4]','k_x_[-]','k_y_[-]','A_[m^2]','pitch_[deg]','x_e_[m]','y_e_[m]']
            f.write('\t'.join(['{:20s}'.format(s) for s in cols])+'\n')
            f.write('$%i %i\n' % (1, dfStructure.shape[0]))
            f.write('\n'.join('\t'.join('%19.13e' %x for x in y) for y in dfStructure.values))

    # --- Rewrite htc file
    if H2_htcfile is not None:
        def readToMarker(lines, marker, i, nMax=None, noException=False):
            l_sel=[]
            if nMax is None: nMax=len(lines)
            while i<nMax:
                line=lines[i]
                if line.replace(' ','').lower().find(marker)>=0:
                    break
                l_sel.append(line.strip())
                i+=1
            if line.strip().replace(' ','').lower().find(marker)<0:
                if noException:
                    return None, None, None
                else:
                    raise Exception('Marker not found '+ marker)
            return l_sel, line, i

        with open(H2_htcfile, 'r') as f:
            lines_in = f.readlines()
        lines_out = []
        bodyNotFound=True
        iBodyEnd=0
        nBodies=0
        while bodyNotFound and nBodies<10:
            _, line, iBodyStart = readToMarker(lines_in, 'beginmain_body',iBodyEnd)
            _, line, iBodyEnd = readToMarker(lines_in, 'endmain_body', iBodyStart)
            _, line, iBody = readToMarker(lines_in, 'name'+bodyname, iBodyStart, iBodyEnd, True)
            nBodies+=1
            if line is None:
                iBody=-1
            else:
                #print('Body {} found between lines {} and {} '.format(bodyname, iBodyStart+1, iBodyEnd+1))
                bodyNotFound=False
        if nBodies>=10:
            raise Exception('Body {} not found in file'.format(bodyname))

        _, line, iC2Start = readToMarker(lines_in, 'beginc2_def', iBodyStart, iBodyEnd)
        _, line, iC2End   = readToMarker(lines_in, 'endc2_def'  , iC2Start, iBodyEnd)

        _, line, iTIStart = readToMarker(lines_in, 'begintimoschenko_input', iBodyStart, iBodyEnd)
        _, line, iTIEnd   = readToMarker(lines_in, 'endtimoschenko_input'  , iTIStart, iBodyEnd)


        simdir        = os.path.dirname(H2_htcfile)
        H2_stfile_rel = os.path.relpath(H2_stfile, simdir)

        lines_out  = lines_in[:iTIStart+1]
        lines_out += ['      filename {};\n'.format(H2_stfile_rel)]
        if FPM:
            lines_out += ['      FPM 1 ;\n']
        #    lines_out += ['      FPM 0 ;\n']
        lines_out += ['      set 1 1 ;\n']
        lines_out += lines_in[iTIEnd:iC2Start+1]
        lines_out += ['      nsec {} ;\n'.format(dfMeanLine.shape[0])]
        for i, row in dfMeanLine.iterrows():
            lines_out += ['      sec {:4d}\t{:13.6e}\t{:13.6e}\t{:13.6e}\t{:13.6e};\n'.format(i+1, row['x_[m]'],row['y_[m]'],row['z_[m]'],row['twist_[deg]'])]
        lines_out += lines_in[iC2End:]

        with open(H2_htcfile, 'w') as f:
            f.write(''.join(lines_out))

    return dfMeanLine, dfStructure


def beamDyn2Hawc2FPM_raw(r_bar, kp_x, kp_y, kp_z, twist,
        m, Ixi, Iyi, x_G, y_G, theta_i, 
        x_C, y_C, theta_p,  
        K):
    """
    NOTE: all angles are in radians
    
    """
    import scipy.linalg
    # --- BeamDyn to Hawc2 Structural data
    # Hawc2 = BeamDyn
    x_cg    = -y_G
    y_cg    = x_G
    x_e     = -y_C
    y_e     = x_C
    pitch   = theta_p*180/np.pi # [deg] NOTE: could use theta_p, theta_i or theta_s
    if np.all(np.abs(m)<1e-16):
        ri_y    = m*0
        ri_x    = m*0
    else:
        ri_y    = np.sqrt(Ixi/m)    # [m]
        ri_x    = np.sqrt(Iyi/m)    # [m]
    # Curvilinear position of keypoints (only used to get max radius...)
    dr= np.sqrt((kp_x[1:]-kp_x[0:-1])**2 +(kp_y[1:]-kp_y[0:-1])**2 +(kp_z[1:]-kp_z[0:-1])**2)
    r_p= np.concatenate(([0],np.cumsum(dr)))
    r=r_bar * r_p[-1]

    RotMat=np.array([  # From Hawc2 to BeamDyn
            [0 ,1,0],
            [-1,0,0],
            [0,0,1]])
    RR= scipy.linalg.block_diag(RotMat,RotMat)

    nSpan = len(K[0,0])
    KH2=np.zeros((6,6,nSpan))
    for iSpan in np.arange(nSpan):
        Kbd = np.zeros((6,6))
        for i in np.arange(6):
            for j in np.arange(6):
                Kbd[i,j] = K[i,j][iSpan]
        Kh2 = (RR.T).dot(Kbd).dot(RR)
        for i in np.arange(6):
            for j in np.arange(6):
                KH2[i,j][iSpan]=Kh2[i,j]

    K11 = KH2[0,0]
    K22 = KH2[1,1]
    K33 = KH2[2,2]
    K44 = KH2[3,3]
    K55 = KH2[4,4]
    K66 = KH2[5,5]

    K12 = KH2[0,1]
    K13 = KH2[0,2]
    K14 = KH2[0,3]
    K15 = KH2[0,4]
    K16 = KH2[0,5]
    K23 = KH2[1,2]
    K24 = KH2[1,3]
    K25 = KH2[1,4]
    K26 = KH2[1,5]
    K34 = KH2[2,3]
    K35 = KH2[2,4]
    K36 = KH2[2,5]
    K44 = KH2[3,3]
    K45 = KH2[3,4]
    K46 = KH2[3,5]
    K55 = KH2[4,4]
    K56 = KH2[4,5]


    columns=['r_[m]','m_[kg/m]','x_cg_[m]','y_cg_[m]','ri_x_[m]','ri_y_[m]','pitch_[deg]','x_e_[m]','y_e_[m]','K11','K12','K13','K14','K15','K16','K22','K23','K24','K25','K26','K33','K34','K35','K36','K44','K45','K46','K55','K56','K66']
    data = np.column_stack((r, m, x_cg, y_cg, ri_x, ri_y, pitch, x_e, y_e, K11,K12,K13,K14,K15,K16,K22,K23,K24,K25,K26,K33,K34,K35,K36,K44,K45,K46,K55,K56,K66))
    dfStructure = pd.DataFrame(data=data, columns=columns)

#     #      Siemens z ->  BeamDyn z
#     #      Siemens x -> -BeamDyn y
#     #      Siemens y ->  BeamDyn x
#     KSiemens = np.zeros((nSpan,6,6))
#     for i in np.arange(6):
#         for j in np.arange(6):
#             if j>=i:
#                 key='d{}{}'.format(i+1,j+1)
#             else:
#                 key='d{}{}'.format(j+1,i+1)
#             KSiemens[:,i,j] = sd[key]
# 
#     for i in np.arange(len(bp['span'])):
#         K = bp['K'][i]*0
#         Ks= KSiemens[i]
#         K = RR.dot(Ks).dot(RR.T)
#         bp['K'][i] = K



    # --- BeamDyn to Hawc2 Reference axis
    X_H2     = -kp_y
    Y_H2     = kp_x
    Z_H2     = kp_z
    twist_H2 = - twist*180/np.pi # - [deg]
    columns=['x_[m]', 'y_[m]', 'z_[m]', 'twist_[deg]']
    data = np.column_stack((X_H2, Y_H2, Z_H2, twist_H2))
    dfMeanLine = pd.DataFrame(data=data, columns=columns)

    return dfMeanLine, dfStructure



def beamDyn2Hawc2_raw(r_bar, kp_x, kp_y, kp_z, twist,
        m, Ixi, Iyi, x_G, y_G, theta_i, 
        EA, EIx, EIy, GKt, kxsGA, kysGA, x_C, y_C, theta_p, x_S, y_S, theta_s, 
        A=None, E=None, G=None):
    """ 
    NOTE: all angles are in radians
    """
    # --- BeamDyn to Hawc2 Structural data
    if A is None: A = np.ones(x_G.shape)
    if E is None: E = EA/A
    if G is None: G = E/2/(1+0.3) # Young modulus
    # Hawc2 = BeamDyn
    x_cg    = -y_G
    y_cg    = x_G
    x_sh    = -y_S
    y_sh    = x_S
    x_e     = -y_C
    y_e     = x_C
    I_y     = EIx/E            # [m^4] Hawc2 Iy is wrt to principal bending ye axis
    I_x     = EIy/E            # [m^4] Hawc2 Ix is wrt to principal bending xe axis
    I_p     = GKt/G            # [m^4]
    k_y     = kxsGA/(G*A)
    k_x     = kysGA/(G*A)
    pitch   = theta_p*180/np.pi # [deg] NOTE: could use theta_p, theta_i or theta_s
    if np.all(np.abs(m)<1e-16):
        ri_y    = m*0
        ri_x    = m*0
    else:
        ri_y    = np.sqrt(Ixi/m)    # [m]
        ri_x    = np.sqrt(Iyi/m)    # [m]
    # Curvilinear position of keypoints (only used to get max radius...)
    dr= np.sqrt((kp_x[1:]-kp_x[0:-1])**2 +(kp_y[1:]-kp_y[0:-1])**2 +(kp_z[1:]-kp_z[0:-1])**2)
    r_p= np.concatenate(([0],np.cumsum(dr)))
    r=r_bar * r_p[-1]

    columns=['r_[m]','m_[kg/m]','x_cg_[m]','y_cg_[m]','ri_x_[m]','ri_y_[m]','x_sh_[m]','y_sh_[m]','E_[N/m^2]','G_[N/m^2]','I_x_[m^4]','I_y_[m^4]','I_p_[m^4]','k_x_[-]','k_y_[-]','A_[m^2]','pitch_[deg]','x_e_[m]','y_e_[m]']
    data = np.column_stack((r,m,x_cg,y_cg,ri_x, ri_y, x_sh, y_sh, E, G, I_x, I_y, I_p, k_x, k_y, A, pitch, x_e, y_e))
    dfStructure = pd.DataFrame(data=data, columns=columns)

    # --- BeamDyn to Hawc2 Reference axis
    X_H2     = -kp_y
    Y_H2     = kp_x
    Z_H2     = kp_z
    twist_H2 = - twist*180/np.pi # -[deg]
    columns=['x_[m]', 'y_[m]', 'z_[m]', 'twist_[deg]']
    data = np.column_stack((X_H2, Y_H2, Z_H2, twist_H2))
    dfMeanLine = pd.DataFrame(data=data, columns=columns)

    return dfMeanLine, dfStructure




# --------------------------------------------------------------------------------}
# --- Functions for 6x6 matrices 
# --------------------------------------------------------------------------------{
def M66toProps(M, convention='BeamDyn'):
    """ 
    Convert mass properties of a 6x6 section to beam properties
    This assumes that the axial and bending loads are decoupled.

    INPUTS:
     - M : 6x6 array of mass elements. Each element may be an array (e.g. for all spanwise values)
     OUTPUTS:
      - m: section mass
      - Ixx, Iyy, Ixy: area moment of inertia
      - x_G, y_G
    """
    M11=M[0,0]
    M44=M[3,3]
    M55=M[4,4]
    M66=M[5,5]
    M16=M[0,5]
    M26=M[1,5]
    M45=M[3,4]

    m=M11
    if convention=='BeamDyn':
        if np.all(np.abs(m)<1e-16):
            Ixi = Iyi = Ipp = x_G = y_G = theta_i = m*0
            return m, Ixi, Iyi, Ipp, x_G, y_G, theta_i
        y_G= -M16/m
        x_G=  M26/m
        # sanity
        np.testing.assert_array_almost_equal([M[0,3],M[0,4]],[0*M[0,3],0*M[0,3]])
        np.testing.assert_array_almost_equal([M[1,3],M[1,4]],[0*M[0,3],0*M[0,3]])
        np.testing.assert_array_almost_equal([M[3,5],M[4,5]],[0*M[0,3],0*M[0,3]])
        
        Ixx =  M44-m*y_G**2
        Iyy =  M55-m*x_G**2
        Ixy = -M45-m*x_G*y_G
        Ipp =  M66 -m*(x_G**2 + y_G**2)

        if np.all(np.abs(Ixy)<1e-16):
            # NOTE: Assumes theta_i ==0
            #print('>>> Assume theta_i 0')
            Ixi     = Ixx
            Iyi     = Iyy
            theta_i = Ixx*0
        else:
            #print('>>> Minimize theta_i')
            Ixi = np.zeros(Ixx.shape)
            Iyi = np.zeros(Ixx.shape)
            theta_i = np.zeros(Ixx.shape)
            for i, (hxx, hyy, hxy) in enumerate(zip(Ixx,Iyy,Ixy)):
                Ixi[i],Iyi[i],theta_i[i] = solvexytheta(hxx,hyy,hxy)

                MM2= MM(m[i],Ixi[i],Iyi[i],Ipp[i],x_G[i],y_G[i],theta_i[i])

                np.testing.assert_allclose(MM2[3,3], M[3,3][i], rtol=1e-3)
                np.testing.assert_allclose(MM2[4,4], M[4,4][i], rtol=1e-3)
                np.testing.assert_allclose(MM2[5,5], M[5,5][i], rtol=1e-3)
                np.testing.assert_allclose(MM2[3,4], M[3,4][i], rtol=1e-3)

        np.testing.assert_array_almost_equal(Ipp, Ixx+Iyy, 2)
        np.testing.assert_array_almost_equal(Ipp, Ixi+Iyi, 2)
         
    else:
        raise NotImplementedError()

    return m, Ixi, Iyi, Ipp, x_G, y_G, theta_i

def solvexytheta(Hxx,Hyy,Hxy):
    """ 
     Solve for a system of three unknown, used to get:
       - EIy, EIx and thetap given Hxx,Hyy,Hxy
       - kxs*GA, kys*GA and thetas given Kxx,Kyy,Kxy
       - I_x, I_y and theta_is given Ixx,Iyy,Ixy
    """
    from scipy.optimize import fsolve
    def residual(x):
        EI_x, EI_y, theta_p =x
        res=np.array([
            Hxx - EI_x*np.cos(theta_p)**2 - EI_y*np.sin(theta_p)**2 ,
            Hyy - EI_x*np.sin(theta_p)**2 - EI_y*np.cos(theta_p)**2,
            Hxy - (EI_y-EI_x)*np.sin(theta_p)*np.cos(theta_p)]
                ).astype(float)
        return res
    x0 = [Hxx,Hyy,0]
    x_opt   = fsolve(residual, x0)
    EI_x,EI_y,theta_p = x_opt
    theta_p = np.mod(theta_p,2*np.pi)
    return EI_x, EI_y, theta_p



def K66toProps(K, theta_p_in=None, convention='BeamDyn'):
    """ 
    Convert stiffness properties of a 6x6 section to beam properties
    This assumes that the axial and bending loads are decoupled.

    INPUTS:
     - K : 6x6 array of stiffness elements. Each element may be an array (e.g. for all spanwise values)
    INPUTS OPTIONAL:
     - theta_p_in : angle from section to principal axis [rad], positive around z
     - convention : to change coordinate systems in the future
    OUTPUTS:
     - EA, EIx, EIy: axial and bending stiffnesses
     - kxGA, kyGA, GKt: shear and torsional stiffness
     - xC,yC : centroid
     - xS,yS : shear center
     - theta_p, theta_s: angle to principal axes and shear axes [rad]
    """
    K11=K[0,0]
    K22=K[1,1]
    K33=K[2,2]
    K44=K[3,3]
    K55=K[4,4]
    K66=K[5,5]

    K12=K[0,1]
    K16=K[0,5]
    K26=K[1,5]
    K34=K[2,3]
    K35=K[2,4]
    K45=K[3,4]

    if convention=='BeamDyn':
        # --- EA, EI, centroid, principal axes
        EA =  K33
        yC =  K34/EA
        xC = -K35/EA
        Hxx=  K44-EA*yC**2
        Hyy=  K55-EA*xC**2 # NOTE: xC fixed
        Hxy= -K45-EA*xC*yC # NOTE: sign changed

        if theta_p_in is not None:
            theta_p=theta_p_in
            print('>>> theta_p given')
            C2=np.cos(theta_p)**2
            S2=np.sin(theta_p)**2
            C4=np.cos(theta_p)**4
            S4=np.sin(theta_p)**4
            EIxp = (Hxx*C2 - Hyy*S2)/(C4-S4)
            EIyp = (Hxx*S2 - Hyy*C2)/(S4-C4)
            Hxyb = (EIyp-EIxp)*np.sin(theta_p)*np.cos(theta_p)

            bNZ=np.logical_and(Hxy!=0, Hxyb!=0)
            np.testing.assert_allclose(Hxy[bNZ], Hxyb[bNZ], rtol=1e-3)
            np.testing.assert_allclose(EIxp+EIyp, Hxx+Hyy, rtol=1e-3)

        else:
            if np.all(np.abs(Hxy)<1e-16):
                #print('>>>> assume theta_p=0')
                # NOTE: Assumes theta_p ==0
                EIx = Hxx
                EIy = Hyy
                theta_p=0*EA
            else:
                #print('>>> Minimization for theta_p')
                EIxp= np.zeros(Hxx.shape)
                EIyp= np.zeros(Hxx.shape)
                theta_p = np.zeros(Hxx.shape)
                for i, (hxx, hyy, hxy) in enumerate(zip(Hxx,Hyy,Hxy)):
                    EIxp[i],EIyp[i],theta_p[i] = solvexytheta(hxx,hyy,hxy)

            theta_p[theta_p>np.pi]=theta_p[theta_p>np.pi]-2*np.pi

        # --- Torsion, shear terms, shear center
        Kxx =  K11
        Kxy = -K12
        Kyy =  K22
        yS  = (Kyy*K16+Kxy*K26)/(-Kyy*Kxx + Kxy**2)
        xS  = (Kxy*K16+Kxx*K26)/( Kyy*Kxx - Kxy**2)
        GKt = K66 - Kxx*yS**2 -2*Kxy*xS*yS - Kyy*xS**2
        if np.all(np.abs(Kxy)<1e-16):
            # Assumes theta_s=0
            kxsGA = Kxx # Kxx = kxs*GA
            kysGA = Kyy
            theta_s=0*EA
        else:
            kxsGA = np.zeros(Kxx.shape)
            kysGA = np.zeros(Kxx.shape)
            theta_s = np.zeros(Hxx.shape)
            for i, (kxx, kyy, kxy) in enumerate(zip(Kxx,Kyy,Kxy)):
                kxsGA[i],kysGA[i],theta_s[i] = solvexytheta(kxx,kyy,kxy)

        theta_s[theta_s>np.pi]=theta_s[theta_s>np.pi]-2*np.pi


        # sanity checks
        KK2= KK(EA, EIxp, EIyp, GKt, EA*0+1, kxsGA, kysGA, xC, yC, theta_p, xS, yS, theta_s)
        np.testing.assert_allclose(KK2[0,0], K[0,0], rtol=1e-2)
        np.testing.assert_allclose(KK2[1,1], K[1,1], rtol=1e-2)
        np.testing.assert_allclose(KK2[2,2], K[2,2], rtol=1e-2)
        np.testing.assert_allclose(KK2[3,3], K[3,3], rtol=1e-1)
#         np.testing.assert_allclose(KK2[4,4], K[4,4], rtol=1e-2)
        np.testing.assert_allclose(KK2[5,5], K[5,5], rtol=1e-1)
        np.testing.assert_allclose(KK2[2,3], K[2,3], rtol=1e-2)
        np.testing.assert_allclose(KK2[2,4], K[2,4], rtol=1e-2)

        np.testing.assert_allclose(K16, -Kxx*yS-Kxy*xS)
#         np.testing.assert_allclose(KK2[0,5], K[0,5],rtol=1e-3)
#         np.testing.assert_allclose(KK2[1,5], K[1,5],rtol=5e-2) # Kxy harder to get
        #np.testing.assert_allclose(KK2[3,4], K[3,4]) # <<< hard to match

    else:
        raise NotImplementedError()

    return EA, EIxp, EIyp, kxsGA, kysGA, GKt, xC, yC, xS, yS, theta_p, theta_s

# --------------------------------------------------------------------------------}
# --- Bauchau 
# --------------------------------------------------------------------------------{
def K_sheartorsion_xbeam(J,K22,K23,K33,x2,x3):
    """ Returns Shear-torsion stiffness matrix.  See Eq.(13) of DyMore manual """
    return np.array( [
        [J + K22*x3**2 + 2*K23*x2*x3 + K33*x2**2, -K22*x3 - K23*x2, K23*x3 + K33*x2],
        [-K22*x3 - K23*x2, K22, -K23],
        [K23*x3 + K33*x2, -K23, K33]])

def K_axialbending_xbeam(S,H22,H23,H33,x2,x3):
    """ Returns Axial-Bending stiffness matrix. See Eq.(20) of DyMore manual """
    return np.array([
        [S, S*x3, -S*x2],
        [S*x3, H22 + S*x3**2, -H23 - S*x2*x3],
        [-S*x2, -H23 - S*x2*x3, H33 + S*x2**2]])

# --------------------------------------------------------------------------------}
# --- BeamDyn 
# --------------------------------------------------------------------------------{
def K_axialbending(EA, EI_x, EI_y, x_C=0, y_C=0, theta_p=0):
    """
    Axial bending problem. See KK for notations.
    """
    H_xx = EI_x*cos(theta_p)**2 + EI_y*sin(theta_p)**2 
    H_yy = EI_x*sin(theta_p)**2 + EI_y*cos(theta_p)**2
    H_xy = (EI_y-EI_x)*sin(theta_p)*cos(theta_p)
    return np.array([
        [EA      , EA*y_C             , -EA*x_C            ] ,
        [EA*y_C  , H_xx + EA*y_C**2   , -H_xy - EA*x_C*y_C ] ,
        [-EA*x_C , -H_xy - EA*x_C*y_C , H_yy + EA*x_C**2   ] 
        ])

def K_sheartorsion(GKt, GA, kxs, kys, x_S=0, y_S=0, theta_s=0):
    """
    Shear torsion problem. See KK for notations.
    """
    K_xx = GA * ( kxs*cos(theta_s)**2 + kys*sin(theta_s)**2   ) 
    K_yy = GA * ( kxs*sin(theta_s)**2 + kys*cos(theta_s)**2   )
    K_xy = GA * ( (kys-kxs)*sin(theta_s)*cos(theta_s)         )
    return np.array([
        [K_xx                 , -K_xy               , -K_xx*y_S - K_xy*x_S                             ] ,
        [-K_xy                , K_yy                , K_xy*y_S + K_yy*x_S                              ] ,
        [-K_xx*y_S - K_xy*x_S , K_xy*y_S + K_yy*x_S , GKt + K_xx*y_S**2 + 2*K_xy*x_S*y_S + K_yy*x_S**2 ]
        ])

def KK(EA, EI_x, EI_y, GKt, GA, kxs, kys, x_C=0, y_C=0, theta_p=0, x_S=0, y_S=0, theta_s=0):
    """ 
    Returns 6x6 stiffness matrix at the cross section origin O, based on inputs at centroid and shear center.
    INPUTS:
        - EA, EI_x, EI_y: diagonal terms for the axial bending expressed at the centroid and in the principal axis frame
        - GKt, GA*kxs, GA*kys: diagonal terms for the shear/torsion expressed at the shear center and in the princial shear direction frame
        - kxs, kys: dimensionless shear parameters
        - x_C, y_C: coordinates of the centroid (elastic center/ neutral axis), expressed FROM the origin of the cross section O
        - x_S, y_S:       "            shear center            "                  "                                             
        - theta_p : angle (around z) FROM the reference axes to the principal axes [rad]
        - theta_s :       "            "             "              principal shear axes [rad]
    """
    H_xx = EI_x*cos(theta_p)**2 + EI_y*sin(theta_p)**2 
    H_yy = EI_x*sin(theta_p)**2 + EI_y*cos(theta_p)**2
    H_xy = (EI_y-EI_x)*sin(theta_p)*cos(theta_p)
    K_xx = GA * ( kxs*cos(theta_s)**2 + kys*sin(theta_s)**2   ) 
    K_yy = GA * ( kxs*sin(theta_s)**2 + kys*cos(theta_s)**2   )
    K_xy = GA * ( (kys-kxs)*sin(theta_s)*cos(theta_s)         )
    return np.array([
        [K_xx                 , -K_xy               , 0*EA    , 0*EA               , 0*EA               , -K_xx*y_S - K_xy*x_S                             ]    , 
        [-K_xy                , K_yy                , 0*EA    , 0*EA               , 0*EA               , K_xy*y_S + K_yy*x_S                              ]    , 
        [0*EA                 , 0*EA                , EA      , EA*y_C             , -EA*x_C            , 0*EA                                                ] , 
        [0*EA                 , 0*EA                , EA*y_C  , H_xx + EA*y_C**2   , -H_xy - EA*x_C*y_C , 0*EA                                                ] , 
        [0*EA                 , 0*EA                , -EA*x_C , -H_xy - EA*x_C*y_C , H_yy + EA*x_C**2   , 0*EA                                                ] , 
        [-K_xx*y_S - K_xy*x_S , K_xy*y_S + K_yy*x_S , 0*EA    , 0*EA               , 0*EA               , GKt + K_xx*y_S**2 + 2*K_xy*x_S*y_S + K_yy*x_S**2 ]
        ])

def MM(m,I_x,I_y,I_p,x_G=0,y_G=0,theta_i=0):
    """ 
    Returns the mass matrix at a given point O and with respect to given orientation axes based 
    on the values at the center of gravity and in the inertia axis frame.
    The convention is such that:
      - x_G,y_G      : the distaances FROM point O to point G
      - theta_i      : angle (around z) FROM the reference axes to the inertial axes
      - I_x, I_y, I_p: "diagonal" inertias for the body expressed in the inertial frame and at point G
    """
    Ixx = I_x*cos(theta_i)**2 + I_y*sin(theta_i)**2 
    Iyy = I_x*sin(theta_i)**2 + I_y*cos(theta_i)**2
    Ixy = (I_y-I_x)*sin(theta_i)*cos(theta_i)

    return np.array([
        [m      , 0*m     , 0*m      , 0*m                , 0*m                , -m*y_G]                    , 
        [0*m      , m     , 0*m      , 0*m                , 0*m                , m*x_G]                     , 
        [0*m      , 0*m     , m      , m*y_G            , -m*x_G           , 0*m]                         , 
        [0*m      , 0*m     , m*y_G  , Ixx + m*y_G**2   , -Ixy - m*x_G*y_G , 0*m]                         , 
        [0*m      , 0*m     , -m*x_G , -Ixy - m*x_G*y_G , Iyy + m*x_G**2   , 0*m]                         , 
        [-m*y_G , m*x_G , 0*m      , 0*m                , 0*m                , I_p + m*x_G**2 + m*y_G**2]
        ])


if __name__=='__main__':
    np.set_printoptions(linewidth=300)

	# ---  Hawc2 to BeamDyn
    H2MeanLine     = '../../data/Hawc2/Blade_Planform_Hawc2.csv' # csv file with c2def columns: ['x_[m]','y_[m]','z_[m]','twist_[deg]']
    H2Structure    = '../../data/Hawc2/Blade_Structural_Hawc2.csv' # csv file with columns ['r_[m]','m_[kg/m]','x_cg_[m]','y_cg_[m]','ri_x_[m]',... ,'x_e_[m]','y_e_[m]']
    BDMainTemplate = '../../data/Hawc2/BD.dat'  # template file to write main BD file

    BDOut          = 'BD_Smooth.dat'
    BDBldOut       = 'BD_Blade_Smooth.dat' 
    Mu             = [0.001]*6 # damping
    hawc2ToBeamDyn(H2MeanLine, H2Structure, BDBldOut, BDOut, BDMainTemplate, Mu=Mu, poly_exp=[2,3,4,5], bPlot=False)



    # --- BeamDyn 2 Hawc 2
    BD_mainfile  = 'solid_beam_BeamDyn.dat'
    BD_bladefile = '../solid_beam_BeamDyn_Blade.dat'
    H2_htcfile_old  = './_template.htc'
    H2_htcfile_new  = './solid_beam_hawc2.htc'
    H2_stfile    = './solid_beam_st.dat'

    from shutil import copyfile
    copyfile(H2_htcfile_old, H2_htcfile_new)

    beamDyn2Hawc2(BD_mainfile, BD_bladefile, H2_htcfile_new, H2_stfile, 'beam_1', A=None, E=None, G=None)


