import numpy as np
from numpy import cos, sin
import pandas as pd
import os


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

def hawc2tobeamdyn(H2MeanLineFile, H2StructureFile, BDBldFileOut, BDMainFileOut=None, BDMainTemplate=None, Mu = 1.0e-03, 
                   ref_axis='c2def-polyfit', poly_exp=[2,3,4,5], zref=None, Label='', bPlot=False,
                   bNoOffset=False, bNoPreSweep=False, nRootZeros=0, bCGOnMeanLine=False):  # Experimental options
    """
    Writes BeamDyn inputs files from two csv files derived from "Hawc2" inputs

    INPUTS:
        - H2MeanLineFile:  csv file with one header line, containing c2 def definition, (Hawc2 coordinates)
                           Column order has to be:  ['X_[m]','Y_[m]','Z_[m]','Twist_theta_z_[deg]']
        - H2StructureFile: contains Hawc2 beam structural properties (typically found in a hawc2 st file)
                           Colums order has to be:  ['Radius_[m]','m_[kg/m]','x_cg_[m]','y_cg_[m]','ri_x_[m]',... ,'pitch_[deg]','x_e_[m]','y_e_[m]']
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
    # --- Mean line definition (position and orientation)  - Hawc2 "c2-def file", BeanDyn point "O"
    c2def = pd.read_csv(H2MeanLineFile,delim_whitespace=True)
    # For the equations below to be generic we force the column names
    c2def.columns.values[0:4]=['X_[m]','Y_[m]','Z_[m]','Twist_theta_z_[deg]']

    # --- If necessary, interpolate mean line to user defined positions
    c2def_old = c2def.copy()
    if zref is None:
        # we dont interpolate
        zref=c2def['Z_[m]']
    else:
        z_old = c2def_old['Z_[m]'].values
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
    x_O_h2 =   c2def['Y_[m]'].values  # kp_xr  
    y_O_h2 = - c2def['X_[m]'].values  # kp_yr  
    z_O_h2 =   c2def['Z_[m]'].values  # kp_zr
    twist  = - c2def['Twist_theta_z_[deg]'].values # initial_twist [deg]

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

    x_off = (x_O_h2-x_O) # difference between input_axis and smooth axis
    y_off = (y_O_h2-y_O) # difference between input_axis and smooth axis

    if bNoOffset:
        x_off = 0*x_off  # no offsets
        y_off = 0*y_off


    # --- Cross section properties 
    hwc_in = pd.read_csv(H2StructureFile,delim_whitespace=True)
    # For the equations below to be generic we force the column names
    hwc_in.columns=['Radius_[m]','m_[kg/m]','x_cg_[m]','y_cg_[m]','ri_x_[m]','ri_y_[m]','x_sh_[m]','y_sh_[m]','E_[N/m^2]','G_[N/m^2]','I_x_[m^4]','I_y_[m^4]','I_p_[m^4]','k_x_[-]','k_y_[-]','A_[m^2]','pitch_[deg]','x_e_[m]','y_e_[m]']
    # --- Interpolating to match c2def positions
    hwc = pd.DataFrame(columns=hwc_in.columns)
    z_old = hwc_in['Radius_[m]']
    for c in hwc.columns:
        hwc[c] = np.interp(zref, z_old, hwc_in[c])

    # --- Safety check
    if len(hwc)!=len(c2def):
        raise Exception('Interpolation failed, wrong length. Debug me.')
    if any(hwc['Radius_[m]']!=z_O):
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
                x_G     =  row['y_cg_[m]']  +x_off[i]
                y_G     = -row['x_cg_[m]']  +y_off[i]
            x_S     =  row['y_sh_[m]']  +x_off[i]
            y_S     = -row['x_sh_[m]']  +y_off[i]
            x_C     =  row['y_e_[m]']   +x_off[i]
            y_C     = -row['x_e_[m]']   +y_off[i]

        EA      = row['E_[N/m^2]']* row['A_[m^2]']
        GKt     = row['G_[N/m^2]']* row['I_p_[m^4]']
        GA      = row['G_[N/m^2]']* row['A_[m^2]']
        kxs     = row['k_y_[-]']
        kys     = row['k_x_[-]']
        EIxp    = row['E_[N/m^2]']* row['I_y_[m^4]']   # Should be [N.m^2]
        EIyp    = row['E_[N/m^2]']* row['I_x_[m^4]']
        theta_s = row['pitch_[deg]']*np.pi/180
        theta_p = row['pitch_[deg]']*np.pi/180
        theta_i = row['pitch_[deg]']*np.pi/180

        m       = row['m_[kg/m]']
        Ixi     = row['ri_y_[m]']**2 * m    # [kg.m]              
        Iyi     = row['ri_x_[m]']**2 * m    # [kg.m]
        I_p     = Ixi+Iyi                   # [kg.m]

        #print(EIxp,EIyp)

        M =  MM(m, Ixi, Iyi, I_p, x_G, y_G, theta_i)
        K =  KK(EA, EIxp, EIyp, GKt, GA, kxs, kys, x_C, y_C, theta_p, x_S, y_S, theta_s)
        
        lM.append(M)
        lK.append(K)
        vx_G.append(x_G); vy_G.append(y_G)
        vx_S.append(x_S); vy_S.append(y_S)
        vx_C.append(x_C); vy_C.append(y_C)



    # --- Writing BeamDyn blade file
    span=hwc['Radius_[m]'].values
    s_bar=span/span[-1]
    print('Writing BeamDyn file:',BDBldFileOut)
    write_beamdyn_sections(BDBldFileOut,s_bar,lK,lM,Mu,Label=Label)

    # --- db
    #M=np.column_stack((zref, x_off, y_off))
    #np.savetxt(BDBldFileOut.replace('.dat','offsets.txt'), M, delimiter=',',header='z_[m], xoff_[m], yoff_[m]')

    # --- Writing BeamDyn main file based on template file
    if BDMainTemplate is not None and BDMainFileOut is not None:
        try:
            import welib.weio as weio
        except:
            try:
                import weio
            except:
                print('[WARN] Not writing BD main file without weio library (https://github.com/ebranlard/weio)')
                return

        BD=weio.FASTInFile(BDMainTemplate)
        #print(BD.keys())
        BD.data[1]['value']=Label
        BD['MemberGeom'] = np.column_stack((x_O,y_O,z_O,twist))
        BD['kp_total']   = len(x_O)
        BD['BldFile']    = '"'+os.path.basename(BDBldFileOut)+'"' 
        BD.data[BD.getID('kp_total')+1]['value']= '1 {}'.format(len(x_O))

        print('Writing BeamDyn file:',BDMainFileOut)
        BD.write(BDMainFileOut)


    # ---
    if bPlot:
        import matplotlib.pyplot as plt
        colrs=plt.rcParams['axes.prop_cycle'].by_key()['color']

        EdgStiff= np.array([K[3,3] for K in lK])
        FlpStiff= np.array([K[4,4] for K in lK])

        EIxp    = hwc['E_[N/m^2]']*hwc['I_y_[m^4]']   # Should be [N.m^2]
        EIyp    = hwc['E_[N/m^2]']*hwc['I_x_[m^4]']

#         fig=plt.figure()
        fig,axes = plt.subplots(4, 2, sharex=True, figsize=(12.4,09.)) # (6.4,4.8)
        fig.subplots_adjust(left=0.07, right=0.99, top=0.98, bottom=0.07, hspace=0.25, wspace=0.15)
        for ax in axes.ravel():
            ax.tick_params(direction='in')

        # --- Plot mean line from hawc2 and beamdyn
        x_O_h2 =   c2def_old['Y_[m]']   # kp_xr  
        y_O_h2 =  -c2def_old['X_[m]']   # kp_yr  
        z_O_h2 =   c2def_old['Z_[m]']   # kp_zr
        twist  =  -c2def_old['Twist_theta_z_[deg]'] # initial_twist [deg]
#         fig,axes = plt.subplots(2, 1, sharex=True, figsize=(6.4,4.8)) # (6.4,4.8)
#         ax=axes[0,0]

        axes[0,0].text(0.5, 1.01, 'Mean line x', horizontalalignment='center', verticalalignment='bottom', transform = axes[0,0].transAxes)
        axes[0,1].text(0.5, 1.01, 'Mean line y', horizontalalignment='center', verticalalignment='bottom', transform = axes[0,1].transAxes)
        axes[0,0].plot(z_O, x_O   , '-' , label = 'BD smooth)')
        axes[0,0].plot(z_O, x_O_h2, '--', label = 'H2 c2def', ms=3, color='k')
        axes[0,0].plot(z_O, x_off , ':', label = r'"$\Delta$" to c2def', color=colrs[6])
        axes[0,1].plot(z_O, y_O   , '-' , label = 'BD y (smooth)')
        axes[0,1].plot(z_O, y_O_h2, '--' , label = 'H2 "y"', ms=3, color='k')
        axes[0,1].plot(z_O, y_off , ':', label = 'y_off', color=colrs[6])
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
        axes[1,0].plot(hwc['Radius_[m]'], x_O_h2 + hwc['y_cg_[m]'], 'd' , ms=1, color=colrs[1] , label='HAWC2')
        axes[1,0].plot(hwc['Radius_[m]'], x_O_h2 + hwc['y_sh_[m]'], 's' , ms=1, color=colrs[2] )
        axes[1,0].plot(hwc['Radius_[m]'], x_O_h2 + hwc['y_e_[m]'] , 'o' , ms=1, color=colrs[3] )
        axes[1,1].plot(z_O, y_O     , '-' , label = 'BD y (smooth)')
        axes[1,1].plot(z_O              , y_O    + vy_G           , 'd' , ms=6, color=None, markeredgecolor=colrs[1], markerfacecolor="None", label = 'G (COG)') 
        axes[1,1].plot(z_O              , y_O    + vy_S           , 's' , ms=6, color=None, markeredgecolor=colrs[2], markerfacecolor="None", label = 'S (shear center)') 
        axes[1,1].plot(z_O              , y_O    + vy_C           , 'o' , ms=6, color=None, markeredgecolor=colrs[3], markerfacecolor="None", label = 'C (elastic center)') 
        axes[1,1].plot(hwc['Radius_[m]'], y_O_h2 - hwc['x_cg_[m]'], 'd' , ms=1, color=colrs[1] )
        axes[1,1].plot(hwc['Radius_[m]'], y_O_h2 - hwc['x_sh_[m]'], 's' , ms=1, color=colrs[2] )
        axes[1,1].plot(hwc['Radius_[m]'], y_O_h2 - hwc['x_e_[m]'] , 'o' , ms=1, color=colrs[3] )
#         axes[1,0].set_xlabel('z [m]')
#         axes[1,1].set_xlabel('z [m]')
        axes[1,0].set_ylabel('x [m]')
        axes[1,1].set_ylabel('y [m]')
        axes[1,0].legend(loc='upper right', fontsize=8)
        if 'Relative_thickness_[%]' and 'Chord_[m]' in c2def.columns.values:
            c = c2def['Chord_[m]']
            t = c2def['Relative_thickness_[%]'] *c/100
            axes[1,0].plot(z_O, x_O_h2+c/2*np.sin(twist*np.pi/180), '-', color=[0.5,0.5,0.5] )
            axes[1,0].plot(z_O, x_O_h2-c/2*np.sin(twist*np.pi/180), '-', color=[0.5,0.5,0.5] )
            axes[1,1].plot(z_O, y_O_h2+c/2*np.cos(twist*np.pi/180), '-', color=[0.5,0.5,0.5] )
            axes[1,1].plot(z_O, y_O_h2-c/2*np.cos(twist*np.pi/180), '-', color=[0.5,0.5,0.5] )


        # --- Positions rel to mean line
        axes[2,0].text(0.5, 1.01, r'Pos. wrt. meanline, x', horizontalalignment='center', verticalalignment='bottom', transform = axes[2,0].transAxes)
        axes[2,1].text(0.5, 1.01, r'Pos. wrt. meanline, y', horizontalalignment='center', verticalalignment='bottom', transform = axes[2,1].transAxes)
        axes[2,0].plot(z_O              , x_O-x_O     , '-' , label = 'BD meanline')
        axes[2,0].plot(z_O              , vx_G        , 'd' , ms=6, color=None, markeredgecolor=colrs[1], markerfacecolor="None", label = 'G (COG)') 
        axes[2,0].plot(z_O              , vx_S        , 's' , ms=6, color=None, markeredgecolor=colrs[2], markerfacecolor="None", label = 'S (shear center)') 
        axes[2,0].plot(z_O              , vx_C        , 'o' , ms=6, color=None, markeredgecolor=colrs[3], markerfacecolor="None", label = 'C (elastic center)') 
        axes[2,0].plot(z_O              , x_off       , ':', label = r'"$\Delta$" to c2def', color=colrs[6])
#         axes[1,0].plot(hwc['Radius_[m]'], x_O_h2 + hwc['y_cg_[m]'], 'o' , ms=1, color=colrs[1] , label='HAWC2')
#         axes[1,0].plot(hwc['Radius_[m]'], x_O_h2 + hwc['y_sh_[m]'], 'o' , ms=1, color=colrs[2] )
#         axes[1,0].plot(hwc['Radius_[m]'], x_O_h2 + hwc['y_e_[m]'] , 'd' , ms=1, color=colrs[3] )
        axes[2,1].plot(z_O              , y_O -y_O    , '-' , label = 'BD meanline')
        axes[2,1].plot(z_O              , vy_G        , 'd' , ms=6, color=None, markeredgecolor=colrs[1], markerfacecolor="None", label = 'G (COG)') 
        axes[2,1].plot(z_O              , vy_S        , 's' , ms=6, color=None, markeredgecolor=colrs[2], markerfacecolor="None", label = 'S (shear center)') 
        axes[2,1].plot(z_O              , vy_C        , 'o' , ms=6, color=None, markeredgecolor=colrs[3], markerfacecolor="None", label = 'C (elastic center)') 
        axes[2,1].plot(z_O              , y_off       , ':', label = r'"$\Delta$" to c2def', color=colrs[6])
#         axes[1,1].plot(hwc['Radius_[m]'], y_O_h2 - hwc['x_cg_[m]'], 'o' , ms=1, color=colrs[1] )
#         axes[1,1].plot(hwc['Radius_[m]'], y_O_h2 - hwc['x_sh_[m]'], 's' , ms=1, color=colrs[2] )
#         axes[1,1].plot(hwc['Radius_[m]'], y_O_h2 - hwc['x_e_[m]'] , 'd' , ms=1, color=colrs[3] )
        axes[2,0].set_xlabel('z [m]')
        axes[2,1].set_xlabel('z [m]')
        axes[2,0].set_ylabel(r'$\Delta x$ [m]')
        axes[2,1].set_ylabel(r'$\Delta y$ [m]')
        axes[2,0].legend(loc='upper right', fontsize=8)

        # --- Plot Stiffness
#         ax=fig.add_subplot(111)
        ax=axes[3,0]
        ax.text(0.5, 1.01, 'Stiffnesses', horizontalalignment='center', verticalalignment='bottom', transform = ax.transAxes)
        ax.plot(z_O,EdgStiff,'-' , color=colrs[0], label='Edge Stiffness (K_44)')
        ax.plot(z_O,EIxp    ,'--', color=colrs[0], label='EIx "edge" at elastic center')
        ax.plot(z_O,FlpStiff,'-' , color=colrs[1], label='Flap Stiffness (K_55)')
        ax.plot(z_O,EIyp    ,'--', color=colrs[1], label='EIy "flap" at elastic center')
        ax.set_xlabel('z [m]')
        ax.set_ylabel('Stiffness [Nm^2]')
        ax.legend(fontsize=8)

        fig.savefig(BDMainFileOut.replace('.dat','.png'))

        plt.show()
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
        - theta_p : angle (around z) FROM the reference axes to the principal axes
        - theta_s :       "            "             "              principal shear axes
    """
    H_xx = EI_x*cos(theta_p)**2 + EI_y*sin(theta_p)**2 
    H_yy = EI_x*sin(theta_p)**2 + EI_y*cos(theta_p)**2
    H_xy = (EI_y-EI_x)*sin(theta_p)*cos(theta_p)
    K_xx = GA * ( kxs*cos(theta_s)**2 + kys*sin(theta_s)**2   ) 
    K_yy = GA * ( kxs*sin(theta_s)**2 + kys*cos(theta_s)**2   )
    K_xy = GA * ( (kys-kxs)*sin(theta_s)*cos(theta_s)         )
    return np.array([
        [K_xx                 , -K_xy               , 0       , 0                  , 0                  , -K_xx*y_S - K_xy*x_S                             ] ,
        [-K_xy                , K_yy                , 0       , 0                  , 0                  , K_xy*y_S + K_yy*x_S                              ] ,
        [0                    , 0                   , EA      , EA*y_C             , -EA*x_C            , 0                                                ] ,
        [0                    , 0                   , EA*y_C  , H_xx + EA*y_C**2   , -H_xy - EA*x_C*y_C , 0                                                ] ,
        [0                    , 0                   , -EA*x_C , -H_xy - EA*x_C*y_C , H_yy + EA*x_C**2   , 0                                                ] ,
        [-K_xx*y_S - K_xy*x_S , K_xy*y_S + K_yy*x_S , 0       , 0                  , 0                  , GKt + K_xx*y_S**2 + 2*K_xy*x_S*y_S + K_yy*x_S**2 ] 
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
        [m      , 0     , 0      , 0                , 0                , -m*y_G]                    , 
        [0      , m     , 0      , 0                , 0                , m*x_G]                     , 
        [0      , 0     , m      , m*y_G            , -m*x_G           , 0]                         , 
        [0      , 0     , m*y_G  , Ixx + m*y_G**2   , -Ixy - m*x_G*y_G , 0]                         , 
        [0      , 0     , -m*x_G , -Ixy - m*x_G*y_G , Iyy + m*x_G**2   , 0]                         , 
        [-m*y_G , m*x_G , 0      , 0                , 0                , I_p + m*x_G**2 + m*y_G**2]
        ])


if __name__=='__main__':
    np.set_printoptions(linewidth=300)

    H2MeanLine     = '../../_data/Hawc2/Blade_Planform_Hawc2.csv' # csv file with c2def columns: ['X_[m]','Y_[m]','Z_[m]','Twist_theta_z_[deg]']
    H2Structure    = '../../_data/Hawc2/Blade_Structural_Hawc2.csv' # csv file with columns ['Radius_[m]','m_[kg/m]','x_cg_[m]','y_cg_[m]','ri_x_[m]',... ,'x_e_[m]','y_e_[m]']
    BDMainTemplate = '../../_data/Hawc2/BD.dat'  # template file to write main BD file

    BDOut          = 'BD_Smooth.dat'
    BDBldOut       = 'BD_Blade_Smooth.dat' 
    Mu             = [0.001]*6 # damping
    hawc2tobeamdyn(H2MeanLine, H2Structure, BDBldOut, BDOut, BDMainTemplate, Mu=Mu, poly_exp=[2,3,4,5], bPlot=False)


