""" 
Tools to compare two python object/arrays/values etc.

See function compare in this module

"""

import pandas as pd
import numpy as np


def ok(s):
    print('[ OK ] '+s)
    return True

def fail(s):
    print('[FAIL] '+s)
    return False


def compare(o1,o2,Columns=None,tol=1e-8,n1='o1',n2='o2'):
    """ 
    o1 and o2 are two python objects to be compared
    """
    Columns = [] if Columns is None else Columns

    b=False
    if type(o1)!=type(o2):
        fail('{} {} Type mismatch, First is {}, second is {}'.format(n1,n2,type(o1),type(o2)));
        
    if isinstance(o1,pd.core.series.Series):
        b=compare_pandas_series(o1,o2,Columns=[],tol=tol,n1=n1,n2=n1)
    elif isinstance(o1,pd.core.frame.DataFrame):
        b=compare_pandas_df(o1,o2,Columns=[],tol=tol,n1=n1,n2=n1)
    elif isinstance(o1,np.ndarray):
        b=compare_ndarray(o1,o2,tol=tol,n1=n1,n2=n1)
    else:
        raise Exception('Type not implemented yet in compare function {}'.format(type(o1)))
    return b


def compare_len(o1,o2,n1='o1',n2='o2'):
    if len(o1)!=len(o2):
        Msg='Different length {}:{} {}:{}'.format(n1,len(o1),n2,len(o2))
        return fail(Msg)
    else:
        return True


def compare_ndarray(o1,o2,tol=1e-8,n1='o1',n2='o2'):
    b=False
    
    if not compare_len(o1,o2,n1,n2):
        return b
    else:
        AbsErr   = abs(o1-o2);
        # --- RelError TODO TODO Switch on method
        sref=(abs(o1)+abs(o2))/2;
        #RelErr2  = abs(o1-o2)/o2;
        # --- Handling division by zero
        bZero=sref==0;
        sref[bZero]=1
        RelErr   = AbsErr/sref;
        # --- Handling signals very close to zero
        myEps=1e-8
        bSmall = sref<myEps
        RelErr[bSmall]=AbsErr[bSmall]/(myEps)
        # --- Handling difference below machine precision
        eps=np.finfo(o1[0]).eps
        RelErr[AbsErr<myEps]=0;
        #   
        MaxRelErr   = np.max(abs(RelErr))
        if  MaxRelErr > tol:
            b=fail(n1+' tolerance not matched')
        else:
            b=ok(n1+'\t\t tol. matched ({:e} > {:e})'.format(tol,MaxRelErr))

def compare_pandas_df(o1,o2,Columns=None,tol=1e-8,n1='o1',n2='o2'):
    Columns = [] if Columns is None else Columns
    b=False
    if not compare_len(o1,o2,n1,n2):
        return b
    elif len(o1.columns)!=len(o2.columns):
        Msg='Different length {}:{} {}:{}'.format(n1,len(o1),n2,len(o2))
        fail(Msg)
        return b
    else:
        for j,col in enumerate(o1.columns.values):
            v1=o1[col].values
            v2=o2[col].values
            b=compare_ndarray(v1,v2,tol=tol,n1='{}[''{}'']'.format(n1,col),n2='{}[''[]'']'.format(n1,col))

def compare_pandas_series(o1,o2,Columns=None,tol=1e-8,n1='o1',n2='o2'):
    Columns = [] if Columns is None else Columns
    b=False
    #b=compare_ndarray(o1.values,o2.values,tol=tol,n1='{}[{}]'.format(n1,col),n2='{}[[]]'.format(n1,col))
    return False

if __name__ == '__main__':
    a=np.linspace(0,1,10)
    b=np.linspace(0,1,10)
    compare(a,b)

    df1=pd.DataFrame.from_items([
        ("Time", np.linspace(0.1,1,10)+1e-8),
        ("Value", np.random.normal(0.1, 1, 10)),
    ])
    df2=df1.copy()
    df2['Time']=df2['Time']+2e-2
    compare(df1,df2)
