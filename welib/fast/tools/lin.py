"""
Tools to manipulate linear models from OpenFAST lin files
"""
import numpy as np
def unit(s):
    iu=s.rfind('[')
    if iu>1:
        return s[iu+1:].replace(']','')
    else:
        return ''
def no_unit(s):
    s=s.replace('_[',' [')
    iu=s.rfind(' [')
    if iu>1:
        return s[:iu]
    else:
        return s

def matToSIunits(Mat, name='', verbose=False):
    """ 
    Scale a matrix (pandas dataframe) such that is has standard units

    """
    # TODO columns
    for irow,row in enumerate(Mat.index.values):
        # Scaling based on units
        u  = unit(row).lower()
        nu = no_unit(row)
        if u in['deg','deg/s','deg/s^2']:
            if verbose:
                print('Mat {} - scaling deg2rad   for row {}'.format(name,row))
            Mat.iloc[irow,:] /=180/np.pi # deg 2 rad
            Mat.index.values[irow]=Mat.index.values[irow].replace('rad','deg')
        elif u=='rpm':
            if verbose:
                print('Mat {} - scaling rpm2rad/s for row {}'.format(name,row))
            Mat.iloc[irow,:] /=60/(2*np.pi) # rpm 2 rad/s
            Mat.index.values[irow]=nu+'_[rad/s]'
        elif u=='knm':
            if verbose:
                print('Mat {} - scaling kNm to Nm for row {}'.format(name,row))
            Mat.iloc[irow,:] /=1000 # to Nm
            Mat.index.values[irow]=nu+'_[Nm]'
    return Mat


def subMat(df, rows=None, cols=None, check=True):
    """ Extract relevant part from a dataframe, perform a safety check """
    if rows is None:
        rows=df.index
    if cols is None:
        cols=df.columns
    if check:
        missingRows = [l for l in rows if l not in df.index]
        missingCols = [c for c in cols  if c not in df.columns]
        if len(missingRows)>0:
            raise Exception('The following rows are missing from outputs: {}'.format(missingRows))
        if len(missingCols)>0:
            raise Exception('The following columns are missing from inputs: {}'.format(missingCols))
    return df[cols].loc[rows].copy()

def matLabelNoUnit(df):
    """ remove unit from index and columns of matrix"""
    rows = [no_unit(row) for row in df.index.values]
    cols = [no_unit(col) for col in df.columns.values]
    df.index =  rows
    df.columns = cols
    return df

def matLabelReplace(df, s1, s2):
    """ remove unit from index and columns of matrix"""
    rows = [row.replace(s1,s2) for row in df.index.values]
    cols = [col.replace(s1,s2) for col in df.columns.values]
    df.index =  rows
    df.columns = cols
    return df


def matSimpleStateLabels(df, inplace=False):
    RenameMap={
    'PtfmSurge_[m]'       : 'x',
    'PtfmSway_[m]'        : 'y',
    'PtfmHeave_[m]'       : 'z',
    'PtfmRoll_[rad]'      : 'phi_x',
    'PtfmPitch_[rad]'     : 'phi_y',
    'PtfmYaw_[rad]'       : 'phi_z',
    'psi_rot_[rad]'       : 'psi',
    'qt1FA_[m]'           : 'q_FA1',
    'd_PtfmSurge_[m/s]'   : 'dx',
    'd_PtfmSway_[m/s]'    : 'dy',
    'd_PtfmHeave_[m/s]'   : 'dz',
    'd_PtfmRoll_[rad/s]'  : 'dphi_x',
    'd_PtfmPitch_[rad/s]' : 'dphi_y',
    'd_PtfmYaw_[rad/s]'   : 'dphi_z',
    'd_psi_rot_[rad/s]'   : 'dpsi',
    'd_qt1FA_[m/s]'       : 'dq_FA1',
    }
    return df.rename(columns=RenameMap, index=RenameMap, inplace=inplace)

