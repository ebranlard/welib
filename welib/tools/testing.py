import numpy as np

from welib.tools.stats import comparison_stats, allclose_errors
from welib.tools.strings import FAIL

def compareSigWithStats(y1, y2, t=None, sig=None, epsTolP=0.01, rtol=1e-3, atol=None, atolP=0.01, printStats=False, test=True, factor=1):
    """ 
    Provide either:
       y1, y2, t
    or 
       y1=df1
       y2=df2
       sig

    - y1: actual
    - y2: desired
    - epsTolP: Tolerance for relative error, default epsTolP=0.01 = 1% 
    - aTopP:  Tolerance for absolute error in percent of range, default atolP=0.01 = 1% 
    """
    if t is None:
        t1, y1 = y1['Time_[s]'].values ,y1[sig].values # Actual
        t2, y2 = y2['Time_[s]'].values ,y2[sig].values # Desired
        if sig is None:
            sig = 'Unknown' 
    else:
        t1=t
        t2=t
        pass
    y1 = y1.copy() /factor
    y2 = y2.copy() /factor

    stats, sStats =  comparison_stats(t2, y2, t1, y1, stats='eps', method='mean', latex=False)

    atol_e, rtol_e = allclose_errors(y1, y2)

    if printStats:
        print(f'{sig:20s}: atol:{atol_e:.4f} - rtol:{rtol_e:.4f} - {sStats}')

    try:
        # np.testing.assert_array_less(1-np.abs(stats['sigRatio']), 0.08)
        np.testing.assert_array_less(np.abs(stats['eps']), epsTolP*100)
    except:
        if test:
            raise
        else:
            FAIL(f'{sig:20s}: atol:{atol_e:.4f} - rtol:{rtol_e:.4f} - {sStats}')
        # np.testing.assert_array_less(1-np.abs(stats['R2']), 0.08)

        if atol is not None:
            pass
        else:
            signal_range = np.ptp(y2)  # Peak-to-peak range
            atol = atopP * signal_range  # % of total signal range

        np.testing.assert_allclose(y1, y2, rtol=rtol, atol=atol)
    return t1, y1, t2, y2
