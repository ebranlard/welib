import sys, os, yaml
import numpy as np

from yaml.loader import SafeLoader 

import dynamic_stall_mhh as mhh

MyDir=os.path.dirname(__file__)

##########################
# Paths and additional imports from sviv_2d repo
##########################

sviv_2d_path = '../../../../sviv_2d/' # default cloned in same repo

pff_path = os.path.join(MyDir, sviv_2d_path, 'python/PFF')
sys.path.append(pff_path)
import peak_filter_fit as pff
    


def pff_summary(t, x, mode_shapes, nom_freq, dict, aoa,
                half_bandwidth_frac=0.05, tstart=10, remove_end=7, reportnum=30):
    """
    Function to conduct PFF analysis on nc file data

    Inputs:
      t - time series
      x - physical coordinates, displacements
      mode_shapes - 3x3 matrix of mode shapes. Columns are modes=[Flap, Edge, Twist], rows are
                    flap, edge, twist contributions of the given mode. 

      tstart - start of time analysis, probably want to eliminate the first several cycles
      remove_end - remove a number of points off the end of the pff analysis, 
                   for multiharmonic signals this is needed to eliminate filtering end effects
      reportnum - number of values from PFF to report
      

    Outputs: 
      - Various values to save summarizing PFF results
    """

    create_append_dict(dict, 'mode_shapes', mode_shapes.reshape(-1).tolist())
    create_append_dict(dict, 'aoa', float(aoa))

    if np.isnan(x).any():
        create_append_dict(dict,  'flap_mode_amp', (np.ones(reportnum)*np.nan).tolist())
        create_append_dict(dict,  'edge_mode_amp', (np.ones(reportnum)*np.nan).tolist())
        create_append_dict(dict, 'twist_mode_amp', (np.ones(reportnum)*np.nan).tolist())

        create_append_dict(dict,  'flap_mode_damp', (np.ones(reportnum)*np.nan).tolist())
        create_append_dict(dict,  'edge_mode_damp', (np.ones(reportnum)*np.nan).tolist())
        create_append_dict(dict, 'twist_mode_damp', (np.ones(reportnum)*np.nan).tolist())

        create_append_dict(dict,  'flap_mode_freq', (np.ones(reportnum)*np.nan).tolist())
        create_append_dict(dict,  'edge_mode_freq', (np.ones(reportnum)*np.nan).tolist())
        create_append_dict(dict, 'twist_mode_freq', (np.ones(reportnum)*np.nan).tolist())
        
    else:    
    
        # Convert to the modal domain
        modal_q = (np.linalg.inv(mode_shapes) @ x.T).T # use inv of modes to convert to modal domain.
    
        freq_rad_s_q, damp_frac_crit_q, report_t_q, report_amp_q, intermediate_data_q = \
               pff.pff_analysis(t, modal_q, nom_freq, tstart, half_bandwidth_frac, remove_end=remove_end)
    
        create_append_dict(dict, 'flap_mode_amp', report_amp_q[0][-reportnum:].tolist())
        create_append_dict(dict, 'edge_mode_amp', report_amp_q[1][-reportnum:].tolist())
        create_append_dict(dict, 'twist_mode_amp', report_amp_q[2][-reportnum:].tolist())

        create_append_dict(dict,  'flap_mode_damp', damp_frac_crit_q[0][-reportnum:].tolist())
        create_append_dict(dict,  'edge_mode_damp', damp_frac_crit_q[1][-reportnum:].tolist())
        create_append_dict(dict, 'twist_mode_damp', damp_frac_crit_q[2][-reportnum:].tolist())

        create_append_dict(dict,  'flap_mode_freq_rad_s', freq_rad_s_q[0][-reportnum:].tolist())
        create_append_dict(dict,  'edge_mode_freq_rad_s', freq_rad_s_q[1][-reportnum:].tolist())
        create_append_dict(dict, 'twist_mode_freq_rad_s', freq_rad_s_q[2][-reportnum:].tolist())



def create_append_dict(dict, key, val):
    """
    Quick function to either append data to list or start a list for the key
    """

    if key in dict:
        dict[key] += [val]
    else:
        dict[key] = [val]


if __name__ == '__main__':

    velocities = np.array([20, 40, 50])

    aoa = 50 # deg
    t_max = 100 # sec
    t_ramp = 0.0 # sec

    # signal processing parameters
    t_start = 2.0*t_ramp+20 # time to start signal processing
    nominal_freq = (0.5065+0.6935)/2.0
    half_bandwidth_frac = 0.05

    output = './initial_sweep_aoa50.yaml'

    ##########################
    # Create output dictionary
    ##########################
    dict = {}

    ##########################
    # Start Loop over AOA here
    ##########################

    sys.path.append(os.path.join(sviv_2d_path, 'python/ConstructModel'))
    import create_model_funs as create

    struct_file = os.path.join(MyDir, './IEA15_aoa{}_3dof.yaml'.format(aoa))
    bd_yaml = os.path.join(sviv_2d_path, 'python/bd_driver.BD.sum.yaml')
    create.construct_IEA15MW_chord(bd_yaml, struct_file, aoa, node_interest=7)

    # import matplotlib
    # matplotlib.use('MacOSX') # Have to reset the backend after constructing model

    # Mode shapes for processing
    iea15mw_sviv2d = yaml.load(open(struct_file),Loader=yaml.UnsafeLoader)
    mode_shapes = np.array(iea15mw_sviv2d['phi_matrix']).reshape(3,3)

    ##########################
    # Start Loop over velocity here
    ##########################

    for vel in velocities:

        t, xytheta, xytheta_dot = mhh.sviv_2d(vel, struct_file=struct_file, t_max=t_max, 
                                          t_ramp=t_ramp, return_data=True)

        if t.shape[0] == xytheta.shape[1]:
            # Do PFF analysis and augment data dictionary
            print('Need to do analysis here')

            pff_summary(t, xytheta.T, mode_shapes, nominal_freq, dict, aoa,
                       half_bandwidth_frac=half_bandwidth_frac, tstart=t_start)

    ##########################
    # save output dictionary 
    ##########################
    print(dict)

    with open(output, 'w') as outfile:
        yaml.dump(dict, outfile)
