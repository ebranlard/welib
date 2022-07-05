import pandas as pd
import numpy as np

def pd_interp1(x_new, xLabel, df):
    """ Interpolate a panda dataframe based on a set of new value
    This function assumes that the dataframe is a simple 2d-table
    """
    from welib.tools.signal_analysis import multiInterp
    x_old = df[xLabel].values
    data_new=multiInterp(x_new, x_old, df.values.T)
    return pd.DataFrame(data=data_new.T, columns=df.columns.values)
    #nRow,nCol = df.shape
    #nRow = len(xnew)
    #data = np.zeros((nRow,nCol))
    #xref =df[xLabel].values.astype(float)
    #for col,i in zip(df.columns.values,range(nCol)):
    #    yref = df[col].values
    #    if yref.dtype!=float:
    #        raise Exception('Wrong type for yref, consider using astype(float)')
    #    data[:,i] = np.interp(xnew, xref, yref)
    #return pd.DataFrame(data=data, columns = df.columns)

def create_dummy_dataframe(size):
    return pd.DataFrame(data={'col1': np.linspace(0,1,size), 'col2': np.random.normal(0,1,size)})
