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

def radial_avg():

    import fastlib
    import weio


    avgMethod='constantwindow'
    avgParam=5
    filename='test.outb'
    rho=1.225
    R=63
    r_FST_aero=[  TODO ]

    df    = weio.read(filename).toDataFrame()
    dfAvg = fastlib.averageDF(df, avgMethod=avgMethod ,avgParam=avgParam) 

    dfAeroRad   = fastlib.spanwiseAD(dfAvg.iloc[0], r_FST_aero/R, rho , R, nB=3)



if __name__=='__main__':
