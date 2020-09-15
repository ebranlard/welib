import os
import glob
import pandas as pd
import numpy as np

try:
    from welib.fastlib import fastlib
except:
    try:
        import fastlib
    except:
        raise Exception('This script needs welib from https://github.com/ebranlard/welib')

try:
    import weio
except:
    raise Exception('This script needs weio from https://github.com/ebranlard/weio')




def two_files():
    # --- Main Parameters
    fstfiles=['NM90_OF2_1.fst','NM90_OF2_2.fst']
    outputDir='./'
    print(fstfiles)
    for FST_In in fstfiles:
        suffix  = os.path.splitext(os.path.basename(FST_In))[0]
        outfile = os.path.join(outputDir,'Spanwise_Aero_'+suffix+'.csv')
        dfRad=fastlib.spanwisePostPro(FST_In, avgMethod='constantwindow',avgParam=5,out_ext='.outb',outfile=outfile)


def extractSpanTSNew(ts, col_pattern, colname, IR=None):
    """ Helper function to extract spanwise results, like B1N1Cl B1N2Cl etc. 

    Example
        col_pattern: 'B1N(\d*)Cl_\[-\]'
        colname    : 'B1Cl_[-]'
    """
    # Extracting columns matching pattern
    cols, sIdx = find_matching_pattern(ts.keys(), col_pattern)
    if len(cols) ==0:
        return (None,None)

    # Sorting by ID
    cols = np.asarray(cols)
    Idx  = np.array([int(s) for s in sIdx])
    Isort = np.argsort(Idx)
    print(Isort)
    Idx  = Idx[Isort]
    cols = cols[Isort]
    print(Idx)
    print(cols)

    nrMax =  np.max(Idx)
    Values = np.zeros((nrMax,1))
    Values[:] = np.nan
#     if IR is None:
#         cols   = [col_pattern.format(ir+1) for ir in range(nr)]
#     else:
#         cols   = [col_pattern.format(ir) for ir in IR]
    for idx,col in zip(Idx,cols):
        Values[idx-1]=ts[col]
    print(Values)
    nMissing = np.sum(np.isnan(Values))
    if nMissing==nrMax:
        return (None,None)
    if len(cols)<nrMax:
        print(Values)
        print('[WARN] Not all values found for {}, missing {}/{}'.format(colname,nMissing,nrMax))
    if len(cols)>nrMax:
        print('[WARN] More values found for {}, found {}/{}'.format(colname,len(cols),nrMax))
    return (colname,Values)


if __name__=='__main__':

    import welib.fastlib.fastlib as fastlib
    import weio
    import fnmatch
    import re


    avgMethod='constantwindow'
    avgParam=5
#     filename='test.outb'
    filename='Main_HelicalWake.outb'
    rho=1.225
    R=63
#     r_FST_aero=[  TODO ]
    r_bar_FST_aero=None
    R=None

    df    = weio.read(filename).toDataFrame()
    dfAvg = fastlib.averageDF(df, avgMethod=avgMethod ,avgParam=avgParam) 
    print(dfAvg.iloc[0].keys())
    dfAeroRad = fastlib.spanwiseAD(dfAvg.iloc[0], r_bar_FST_aero, rho , R, nB=3)
    print(dfAeroRad)

#     r=vr_bar*R
#     nr=len(vr_bar)
#     Columns     = [('r/R_[-]', vr_bar)]
    # --- Extract radial data
    ts=dfAvg.iloc[0]
#     pattern='B1N(\d*)AOA_\[deg\]'
#     print(find_matching_pattern(ts.keys(), pattern))
#     print(reg_pattern.search('B1N02AOA_[deg]').groups(1))
#     extractSpanTSNew(ts,'B1N(\d*)AOA_\[deg\]','B1Alpha_[deg]')
#     import pdb; pdb.set_trace()
