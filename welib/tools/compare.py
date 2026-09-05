""" 
Tool to compare two python object/arrays/values etc.
"""

import pandas as pd
import numpy as np
from welib.tools.strings import printMat, printVec
from welib.tools.strings import FAIL as FAIL_, OK as OK_


def FAIL(*args, **kwargs):
    FAIL_(*args, **kwargs)
    return False


def OK(*args, verbose=True, **kwargs):
    if verbose:
        OK_(*args, **kwargs)
    return True


def get_user_attributes(obj):
    """Extract non-private, non-callable attributes from a custom class instance."""
    if hasattr(obj, '__dict__'):
        return {k: v for k, v in obj.__dict__.items() if not k.startswith('_')}
    
    attrs = {}
    for k in dir(obj):
        if not k.startswith('_'):
            try:
                val = getattr(obj, k)
                if not callable(val):
                    attrs[k] = val
            except Exception:
                pass
    return attrs


def compare(o1, o2, Columns=None, tol=1e-8, n1='o1', n2='o2', verbose=True):
    """ 
    o1 and o2 are two python objects to be compared
    """
    # Temporary Hack
    if 'Children' in n1:
        return True
    elif 'Connections' in n1:
        return True

    Columns = [] if Columns is None else Columns

    all_ok = True

    # Type check
    if type(o1) != type(o2):
        FAIL(f'{n1} vs {n2}: Type mismatch, First is {type(o1).__name__}, second is {type(o2).__name__}')
        all_ok = False

    # Pandas Series
    if isinstance(o1, pd.Series) or isinstance(o2, pd.Series):
        if isinstance(o1, pd.Series) and isinstance(o2, pd.Series):
            return compare_pandas_series(o1, o2, Columns=Columns, tol=tol, n1=n1, n2=n2, verbose=verbose) and all_ok
        return False

    # Pandas DataFrame
    elif isinstance(o1, pd.DataFrame) or isinstance(o2, pd.DataFrame):
        if isinstance(o1, pd.DataFrame) and isinstance(o2, pd.DataFrame):
            return compare_pandas_df(o1, o2, Columns=Columns, tol=tol, n1=n1, n2=n2, verbose=verbose) and all_ok
        return False

    # NumPy Array
    elif isinstance(o1, np.ndarray) or isinstance(o2, np.ndarray):
        if isinstance(o1, np.ndarray) and isinstance(o2, np.ndarray):
            return compare_ndarray(o1, o2, tol=tol, n1=n1, n2=n2, verbose=verbose) and all_ok
        return False

    # Dictionaries
    elif isinstance(o1, dict) or isinstance(o2, dict):
        if isinstance(o1, dict) and isinstance(o2, dict):
            return compare_dict(o1, o2, Columns=Columns, tol=tol, n1=n1, n2=n2, verbose=verbose) and all_ok
        else:
            # If one is a dict and the other is an object, compare dict against object attributes
            d1 = o1 if isinstance(o1, dict) else get_user_attributes(o1)
            d2 = o2 if isinstance(o2, dict) else get_user_attributes(o2)
            return compare_dict(d1, d2, Columns=Columns, tol=tol, n1=n1, n2=n2, verbose=verbose) and False

    # Sequences
    elif isinstance(o1, (list, tuple)) or isinstance(o2, (list, tuple)):
        if isinstance(o1, (list, tuple)) and isinstance(o2, (list, tuple)):
            return compare_sequence(o1, o2, Columns=Columns, tol=tol, n1=n1, n2=n2, verbose=verbose) and all_ok
        return False

    # Scalars
    elif isinstance(o1, (int, float, np.floating, np.integer)) and isinstance(o2, (int, float, np.floating, np.integer)):
        return compare_scalar(o1, o2, tol=tol, n1=n1, n2=n2, verbose=verbose) and all_ok

    # Strings, Booleans, None
    elif isinstance(o1, (str, bool, type(None))) or isinstance(o2, (str, bool, type(None))):
        if o1 != o2:
            return FAIL(f"{n1} != {n2} ({o1} vs {o2})")
        return OK(f"{n1} == {n2}", verbose=verbose) and all_ok

    # Custom Class Instance Objects
    else:
        return compare_object(o1, o2, Columns=Columns, tol=tol, n1=n1, n2=n2, verbose=verbose) and all_ok


def compare_len(o1, o2, n1='o1', n2='o2', verbose=True):
    if len(o1) != len(o2):
        Msg = f'Different length {n1}:{len(o1)} {n2}:{len(o2)}'
        return FAIL(Msg)
    return True


def compare_scalar(o1, o2, tol=1e-8, n1='o1', n2='o2', verbose=True):
    abs_err = abs(o1 - o2)
    sref = (abs(o1) + abs(o2)) / 2.0
    rel_err = abs_err / 1e-8 if sref < 1e-8 else abs_err / sref

    if abs_err < 1e-8:
        rel_err = 0.0

    if rel_err > tol:
        return FAIL(f"{n1} scalar mismatch: {o1} vs {o2} (rel_err={rel_err:.2e} > tol={tol:.2e})")
    else:
        return OK(f"{n1}\t\t scalar matched ({rel_err:.2e} <= {tol:.2e})", verbose=verbose)


def compare_ndarray(o1, o2, tol=1e-8, n1='o1', n2='o2', verbose=True):
    if not compare_len(o1, o2, n1, n2, verbose=verbose):
        return False
    if o1.shape != o2.shape:
        return FAIL(f"Shape mismatch {n1}:{o1.shape} {n2}:{o2.shape}")
    if o1.size == 0:
        return OK(f"{n1}\t\t (empty)", verbose=verbose)
    
    # Fallback to element-wise recursive comparison if either array is dtype=object
    if o1.dtype == object or o2.dtype == object:
        all_ok = True
        for idx in np.ndindex(o1.shape):
            idx_str = "".join([f"[{i}]" for i in idx])
            v1 = o1[idx]
            v2 = o2[idx]
            b = compare( v1, v2, tol=tol, n1=n1+idx_str, n2=n2+idx_str, verbose=verbose,)
            if not b:
                all_ok = False
        return all_ok

    # Standard numerical comparison for numeric dtypes
    try:
        AbsErr = abs(o1 - o2)
        sref = (abs(o1) + abs(o2)) / 2.0

        bZero = sref == 0
        sref[bZero] = 1.0
        RelErr = AbsErr / sref

        myEps = 1e-8
        bSmall = sref < myEps
        RelErr[bSmall] = AbsErr[bSmall] / myEps
        RelErr[AbsErr < myEps] = 0.0

        MaxRelErr = np.max(abs(RelErr))
        if MaxRelErr > tol:
            return FAIL(
                f"{n1} tolerance not matched ({MaxRelErr:.2e} > {tol:.2e})"
            )
        else:
            return OK(
                f"{n1}\t\t tol. matched ({MaxRelErr:.2e} <= {tol:.2e})",
                verbose=verbose,
            )
    except TypeError:
        # Safety fallback if vector math fails unexpectedly
        all_ok = True
        for idx in np.ndindex(o1.shape):
            idx_str = "".join([f"[{i}]" for i in idx])
            b = compare(
                o1[idx],
                o2[idx],
                tol=tol,
                n1=n1+idx_str,
                n2=n2+idx_str,
                verbose=verbose,
            )
            if not b:
                all_ok = False
        return all_ok


def compare_pandas_df(o1, o2, Columns=None, tol=1e-8, n1='o1', n2='o2', verbose=True):
    Columns = [] if Columns is None else Columns
    all_ok = True
    
    if not compare_len(o1, o2, n1, n2, verbose=verbose):
        return False
        
    if len(o1.columns) != len(o2.columns):
        FAIL(f'Different number of columns {n1}:{len(o1.columns)} {n2}:{len(o2.columns)}')
        all_ok = False

    cols_to_check = Columns if len(Columns) > 0 else o1.columns.values
    for col in cols_to_check:
        if col not in o2.columns:
            FAIL(f"Column '{col}' in {n1} missing from {n2}")
            all_ok = False
            continue
        v1 = o1[col].values
        v2 = o2[col].values
        b = compare_ndarray(v1, v2, tol=tol, n1=n1+f"['{col}']", n2=n2+f"['{col}']", verbose=verbose)
        if not b:
            all_ok = False
            
    return all_ok


def compare_pandas_series(o1, o2, Columns=None, tol=1e-8, n1='o1', n2='o2', verbose=True):
    if not compare_len(o1, o2, n1, n2, verbose=verbose):
        return False
    return compare_ndarray(o1.values, o2.values, tol=tol, n1=n1, n2=n2, verbose=verbose)


def compare_sequence(o1, o2, Columns=None, tol=1e-8, n1='o1', n2='o2', verbose=True):
    if not compare_len(o1, o2, n1, n2, verbose=verbose):
        return False
    
    all_ok = True
    for i, (v1, v2) in enumerate(zip(o1, o2)):
        b = compare(v1, v2, Columns=Columns, tol=tol, n1=f"{n1}[{i}]", n2=f"{n2}[{i}]", verbose=verbose)
        if not b:
            all_ok = False
            
    return all_ok


def compare_dict(o1, o2, Columns=None, tol=1e-8, n1='o1', n2='o2', verbose=True):
    all_ok = True

    keys1 = set(o1.keys())
    keys2 = set(o2.keys())
    
    missing_in_2 = keys1 - keys2
    missing_in_1 = keys2 - keys1
    
    for k in sorted(missing_in_2):
        FAIL(f"Key '{k}' present in {n1} but missing in {n2}")
        all_ok = False
        
    for k in sorted(missing_in_1):
        FAIL(f"Key '{k}' present in {n2} but missing in {n1}")
        all_ok = False

    common_keys = sorted(keys1.intersection(keys2))
    for k in common_keys:
        b = compare(o1[k], o2[k], Columns=Columns, tol=tol, n1=f"{n1}['{k}']", n2=f"{n2}['{k}']", verbose=verbose)
        if not b:
            all_ok = False

    return all_ok


def compare_object(o1, o2, Columns=None, tol=1e-8, n1='o1', n2='o2', verbose=True):
    all_ok = True

    attrs1 = get_user_attributes(o1)
    attrs2 = get_user_attributes(o2)

    keys1 = set(attrs1.keys())
    keys2 = set(attrs2.keys())

    missing_in_2 = keys1 - keys2
    missing_in_1 = keys2 - keys1

    for k in sorted(missing_in_2):
        FAIL(f"Attribute '{k}' present in {n1} but missing in {n2}")
        all_ok = False

    for k in sorted(missing_in_1):
        FAIL(f"Attribute '{k}' present in {n2} but missing in {n1}")
        all_ok = False

    common_attrs = sorted(keys1.intersection(keys2))
    for k in common_attrs:
        b = compare(attrs1[k], attrs2[k], Columns=Columns, tol=tol, n1=f"{n1}.{k}", n2=f"{n2}.{k}", verbose=verbose)
        if not b:
            all_ok = False

    return all_ok


if __name__ == '__main__':
    class StructA:
        def __init__(self):
            self.x = 10.0
            self.arr = np.linspace(0, 1, 5)
            self.meta = "A"

    class StructB:
        def __init__(self):
            self.x = 10.0001
            self.arr = np.linspace(0, 1, 5) + 1e-9
            self.extra = True

    obj1 = StructA()
    obj2 = StructB()

    print("--- Comparing Custom Classes (Verbose = False: Only Errors/Differences) ---")
    compare(obj1, obj2, tol=1e-5, verbose=False)

    print("\n--- Comparing Custom Classes (Verbose = True: All Matches) ---")
    compare(obj1, obj2, tol=1e-3, verbose=True)
