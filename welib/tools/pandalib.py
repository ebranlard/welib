import pandas as pd
import numpy as np



def pd_interp1(xLabel,df,xnew):
    # Interpolate a panda dataframe based on a set of new value
    nRow,nCol = df.shape
    nRow = len(xnew)
    data = np.zeros((nRow,nCol))
    xref =df[xLabel].values.astype(float)
    for col,i in zip(df.columns.values,range(nCol)):
        yref = df[col].values
        if yref.dtype!=float:
            raise Exception('Wrong type for yref, consider using astype(float)')
        data[:,i] = np.interp(xnew, xref, yref)
    dfi = pd.DataFrame(data=data,columns = df.columns)
    return dfi

def create_dummy_dataframe(size):
    return pd.DataFrame(data={'col1': np.linspace(0,1,size), 'col2': np.random.normal(0,1,size)})
