import os
import numpy as np
import pandas as pd

def getCase(case, solver=None, mid_out=None, **kwargs):

    scriptDir = os.path.dirname(__file__)

    def _opts(default_opts, kwargs):
        opts = default_opts.copy()
        for k,v in kwargs.items():
            if k in opts:
                opts[k] = v
            else:
                raise Exception('Unsuported key {} for case {}'.format(k, case))
        return opts


    # --------------------------------------------------------------------------------}
    # --- Non lifting cases 
    # --------------------------------------------------------------------------------{
    if case=='cylinder':
        default_opts={'m':30, 'U0':1, 'R':2 }
        out = _opts(default_opts, kwargs)

        out['theta']     = -np.linspace(0, 2*np.pi, out['m']+1) # NOTE: counterclockwise
        out['XP']        = out['R']*np.cos(out['theta'])
        out['YP']        = out['R']*np.sin(out['theta'])
        dtheta           = out['theta'][1]-out['theta'][0]
        out['theta_mid'] = out['theta'][:-1]+dtheta/2
        out['Ut']        = 2.0*out['U0']*np.sin(out['theta_mid'])     # gamma_exact = Ut at CP
        out['Cp']        = 1-out['Ut']**2/out['U0']**2
        out['Uxy']       = (out['U0'], 0)
        out['maxVal']    = 2*out['U0'] # For velocity field
        
    elif case=='ellipse':
        default_opts={'m':130, 'U0':1, 'R':1, 'ratio':0.5}
        out = _opts(default_opts, kwargs)
        out['theta']     = -np.linspace(0,2*np.pi, out['m']+1) # NOTE: counterclockwise
        abyr      = np.sqrt((1.-out['ratio'])/(1.+out['ratio'])) # ! see Lewis p 50
        out['XP']        = out['R']             * np.cos(out['theta'])
        out['YP']        = out['R']*out['ratio']* np.sin(out['theta'])
        dtheta           = out['theta'][1]-out['theta'][0]
        out['theta_mid'] = out['theta'][:-1]+dtheta/2

        # Note for ellipst, we can't use "theta" from the CP
        out['Ut']        = (2.0*out['U0']*np.sin(out['theta_mid']))/np.sqrt(1.0+ (abyr**4) - 2.*(abyr**2)*np.cos(2.*out['theta_mid']) ) # See Lewis p50
        out['Cp']        = 1-out['Ut']**2/out['U0']**2
        out['Uxy']       = (out['U0'], 0)
        out['maxVal']    =  1.5 *out['U0'] # For velocity field

    # --------------------------------------------------------------------------------}
    # --- Lifting or nonlfting cases
    # --------------------------------------------------------------------------------{
    elif case in ['VonDeVooren_lift', 'VonDeVooren']: # Sharp
        # Test case 2 - Van de Vooren (Katz Plotkin)
        # NUMBER OF AIRFOIL PANELS, M    :   90
        # THE ANGLE OF ATTACK IN DEGREES :   5
        # THICKNESS COEFF. Eps (<1)      :   0.075
        # T.E. ANGLE COEFF. K (1-2)      :   1.90555555555
        alpha=0
        if 'lift' in case :
            alpha=5
        default_opts={'alpha':alpha*np.pi/180, 'U0':1}
        out = _opts(default_opts, kwargs)
        airfoil_file = os.path.join(scriptDir,'data/VonDeVooren_esp0.075_k1.906_AFOIL2.csv')
        df = pd.read_csv(airfoil_file)
        out['XP'], out['YP'] = df['x'].values, df['y'].values
        out['Uxy'] = (out['U0']*np.cos(out['alpha']), out['U0']*np.sin(out['alpha'])) # Freestream velocity vector [m/s]

        if alpha==5:
            # ref
            if solver is None:
                df_ref = pd.read_csv(os.path.join(scriptDir, 'data/VonDeVooren_Cp_theory.csv')) # KatzPlotkin example - VanDeVooren
            elif solver=='LVP':
                df_ref  = pd.read_csv(os.path.join(scriptDir, 'data/VonDeVooren_Cp_lvortex.csv')) # KatzPlotkin example - VanDeVooren
            else:
                raise NotImplementedError(solver)
            out['x'] = df_ref['x']
            out['Cp'] = df_ref['Cp']
            # Mid
            if mid_out: 
                x = df_ref['x'] .values
                Cp = df_ref['Cp'] .values
                out['x'] = (x[1:] + x[:-1])/2
                out['Cp'] = (Cp[1:] + Cp[:-1])/2

    elif case in ['NACA2412',  'NACA2412_lift']: # Slighlty blunt
        # Test case 3 - XFoil NACA 2412 PPAR N 170 P 4 T 1 R 1
        alpha=0
        if 'lift' in case :
            alpha=5
        default_opts={'alpha':alpha*np.pi/180, 'U0':1}
        out = _opts(default_opts, kwargs)
        airfoil_file = os.path.join(scriptDir,'data/NACA2412.txt')
        df = pd.read_csv(airfoil_file)
        mid = (df['x'].values[1:] + df['x'].values[:-1]) / 2
        out['x_ref'] = mid

        out['XP'], out['YP'] = df['x'].values, df['y'].values
        out['Uxy'] = (out['U0']*np.cos(out['alpha']), out['U0']*np.sin(out['alpha'])) # Freestream velocity vector [m/s]
        
        if alpha==5:
            import pickle
            with open(os.path.join(scriptDir, 'tests/LSN_LV1_NACA2412.pkl'), 'rb') as f:
                data = pickle.load(f) 
                for k,v in data.items():
                    out[k] = v
            out['Cp'] = data['Cp']

    else:
        raise NotImplementedError(case)

    if 'Cp' not in out:
        out['Cp'] = None

    return out 

