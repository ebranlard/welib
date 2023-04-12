"""
Tools to manipulate linear models from OpenFAST lin files
"""
import numpy as np
import pandas as pd

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

def renameList(l, colMap, verbose=False):
    keys = colMap.keys()
    for i,s in enumerate(l):
        if s in keys:
            l[i] = colMap[s] 
        else:
            if verbose:
                print('Label {} not renamed'.format(s))
    return l


def dfToSIunits(Mat, name='', verbose=False):
    def SIscaling(name):
        u  = unit(name).lower()
        nu = no_unit(name)
        if u in['deg','deg/s','deg/s^2']:
            scaling = np.pi/180
            replace = ('deg', 'rad')
        elif u=='rpm':
            scaling = 2*np.pi/60
            replace =  (u, 'rad/s')
        elif u in ['kn']:
            scaling = 1000 # to N
            replace =  (u, 'N')
        elif u in ['kw']:
            scaling = 1000 
            replace = (u, 'W')
        elif u in ['knm', 'kn-m', 'kn*m']:
            scaling = 1000 
            replace = (u, 'Nm')
        else:
            return None, None, None
        newname = nu + '_['+u.replace(replace[0],replace[1])+']'
        return newname, scaling, replace
    newCols = list(Mat.columns.values)
    for icol,col in enumerate(Mat.columns.values):
        # Scaling based on units
        newname, scaling, replace = SIscaling(col)
        if newname is None:
            continue # We skip
        Mat.iloc[:,icol] *= scaling
        newCols[icol] = newname
        if verbose:
            print('Mat {} - scaling {} for col {} > {}'.format(name, col, replace, newname))
    Mat.columns = newCols
    return Mat

def matToSIunits(Mat, name='', verbose=False, row=True, col=True):
    """ 
    Scale a matrix (pandas dataframe) such that is has standard units

    The understanding is that "Columns" are Input and Rows are Outputs
    Because of that, columns are divided by the scaling, rows are multiplied by the scaling



        scalings['rpm']    =  (np.pi/30,'rad/s') 
        scalings['rad' ]   =   (180/np.pi,'deg')
        scalings['kn']     =   (1e3, 'N')
        scalings['kw']     =   (1e3, 'W')
        scalings['knm']    =   (1e3, 'Nm')
        scalings['kn-m']   =   (1e3, 'Nm')
        scalings['kn*m']   =   (1e3, 'Nm')

    """

    def SIscaling(name):
        u  = unit(name).lower()
        nu = no_unit(name)
        if u in['deg','deg/s','deg/s^2']:
            scaling = np.pi/180
            replace = ('deg', 'rad')
        elif u=='rpm':
            scaling = np.pi/30
            replace =  (u, 'rad/s')
        elif u in ['kn']:
            scaling = 1000 # to N
            replace =  (u, 'N')
        elif u in ['kw']:
            scaling = 1000 
            replace = (u, 'W')
        elif u in ['knm', 'kn-m', 'kn*m']:
            scaling = 1000 
            replace = (u, 'Nm')
        else:
            return None, None, None
        newname = nu + '_['+u.replace(replace[0],replace[1])+']'

        return newname, scaling, replace


    if row:
        newIndex = list(Mat.index.values)
        for irow,row in enumerate(Mat.index.values):
            # Scaling based on units
            newname, scaling, replace = SIscaling(row)
            if newname is None:
                continue # We skip
            Mat.iloc[irow,:] *= scaling # NOTE: for row scaling we multiply
            newIndex[irow] = newname
            if verbose:
                print('Mat {} - scaling {} for row {} > {}'.format(name, row, replace, newname))
        Mat.index = newIndex
    if col:
        newCols = list(Mat.columns.values)
        for icol,col in enumerate(Mat.columns.values):
            # Scaling based on units
            newname, scaling, replace = SIscaling(col)
            if newname is None:
                continue # We skip
            Mat.iloc[:,icol] /= scaling # NOTE: for column scaling, we divide!
            newCols[icol] = newname
            if verbose:
                print('Mat {} - scaling {} for col {} > {}'.format(name, col, replace, newname))
        Mat.columns = newCols
    return Mat


def subMat(df, rows=None, cols=None, check=True, name='matrix', removeDuplicates=True):
    """ Extract relevant part from a dataframe, perform a safety check """
    if rows is None:
        rows=df.index
    if cols is None:
        cols=df.columns
    missingRows = [l for l in rows if l not in df.index]
    missingCols = [c for c in cols  if c not in df.columns]
    if check:
        if len(missingRows)>0:
            raise Exception('The following rows are missing from {}: {}'.format(name,missingRows))
        if len(missingCols)>0:
            raise Exception('The following columns are missing from {}: {}'.format(name,missingCols))
    else:
        if len(missingRows)>0:
            print('[WARN] The following rows are missing from {}: {}'.format(name,missingRows))
        if len(missingCols)>0:
            print('[WARN] The following columns are missing from {}: {}'.format(name,missingCols))

    # Create an emtpy dataframe with specifications
    M = np.zeros((len(rows),len(cols)))*np.nan
    df_out = pd.DataFrame(data=M, index=rows, columns=cols)

    # Col/Row that are indeed present in input 
    cols_ = [s for s in cols if s in df.columns]
    rows_ = [s for s in rows if s in df.index]

    # Copy existing 
    df_out.loc[rows_,cols_] = df.loc[rows_, cols_]
    #df = df[cols].loc[rows].copy()
    #if removeDuplicates:
    #    bColDup = df.columns.duplicated()
    #    bIndDup = df.index.duplicated()
    #    df_out = df_out.loc[~bIndDup,~bColDup]
    return df_out

def subSeries(df, rows=None, check=True, name='series', removeDuplicates=True):
    """ Extract relevant part from a dataframe, perform a safety check """
    if rows is None:
        rows=df.index
    missingRows = [l for l in rows if l not in df.index]
    if check:
        if len(missingRows)>0:
            raise Exception('The following rows are missing from {}: {}'.format(name,missingRows))
    else:
        if len(missingRows)>0:
            print('[WARN] The following rows are missing from {}: {}'.format(name,missingRows))
    # Create an emtpy dataframe with specifications
    M = np.zeros(len(rows))*np.nan
    df_out = pd.Series(data=M, index=rows)

    # Col/Row that are indeed present in input 
    rows_ = [s for s in rows if s in df.index]

    # Copy existing 
    df_out.loc[rows_] = df.loc[rows_]
    return df_out


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

