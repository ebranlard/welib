import numpy as np

from welib.weio.fast_input_file import FASTInputFile



def resampleBlade(r_new, fileIn, fileOut, resamplePolars=False):

    bld = FASTInputFile(fileIn)
    M = bld['BldAeroNodes']
    r_old = M[:,0]

    M_new= np.zeros((len(r_new), M.shape[1]))
    for ic in np.arange(1,M.shape[1]):
        M_new[:,ic] = np.interp(r_new, r_old, M[:,ic])
    M_new[:,0] = r_new
    M_new[:,6] = np.around(M_new[:,6]).astype(int)
    bld['BldAeroNodes'] = M_new
    #bld['NumBlNds']     = M_new.shape[0]

    bld.write(fileOut)



