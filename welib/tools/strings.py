import sys
import numpy as np
from html import escape

# ---------- optional libs ----------
try:
    from termcolor import cprint as _tc_cprint
    _HAS_TERMCOLOR = True
except Exception:
    _HAS_TERMCOLOR = False

try:
    from IPython import get_ipython
    from IPython.display import display, HTML
    _IPY = get_ipython()
    # ZMQInteractiveShell => Jupyter Notebook / Lab
    _IN_JUPYTER = _IPY is not None and _IPY.__class__.__name__ == "ZMQInteractiveShell"
except Exception:
    _IN_JUPYTER = False

# --- HTML
_HTML_COLOR = {
    'red': '#d32f2f', 'yellow': '#f7b500', 'green': '#388e3c',
    'blue': '#1976d2', 'magenta': '#8e24aa', 'cyan': '#0097a7', None: 'inherit'
}

# --- ASCII Codes
_ANSI_COLOR = {
    'red':    '\033[91m',
    'yellow': '\033[93m',
    'green':  '\033[92m',
    'blue':   '\033[94m',
    'magenta':'\033[95m',
    'cyan':   '\033[96m',
    None:     ''
}
_ATTRS_ANSI  = {
    'bold':      '\033[1m',
    'underline': '\033[4m'
}
_RESET = '\033[0m'



def cprint_local(msg, color=None, attrs=None, file=sys.stdout, end='\n'):
    color_code = _COLOR.get(color, '')
    attr_code  = ''.join(_ATTR.get(a, '') for a in (attrs or []))
    try:
        print(f"{color_code}{attr_code}{msg}{_RESET}", file=file)
    except Exception:
        # Absolute last resort (no colors, never crash)
        print(msg, file=file, end=end)



def cprint(msg, color=None, attrs=None, file=sys.stdout, end='\n'):
    """Robust colored / bold print. In Jupyter: render HTML for reliable styling.
       In normal terminals: use termcolor if present, else ANSI escapes, else plain print.
       `file` follows print() semantics; when in a notebook and file is stdout/stderr
       the function uses rich HTML output (display)."""
    attrs = attrs or []

    # 1) Jupyter: render HTML so bold + color always show in output cells
    if _IN_JUPYTER and file in (sys.stdout, sys.stderr):
        try:
            color_css = _HTML_COLOR.get(color, color or 'inherit')
            style = ''
            if color_css:
                style += f'color:{color_css};'
            if 'bold' in attrs:
                style += 'font-weight:700;'
            if 'underline' in attrs:
                style += 'text-decoration:underline;'
            safe = escape(msg)
            html = (f"<pre style='margin:0;padding:0;font-family:monospace;"
                    f"white-space:pre-wrap;'><span style=\"{style}\">{safe}</span></pre>")
            display(HTML(html))
            return
        except Exception:
            # fall through to other backends if display fails
            pass

    # 2) termcolor if available (works well in many terminals)
    if _HAS_TERMCOLOR:
        try:
            _tc_cprint(msg, color=color, attrs=attrs, file=file)
            return
        except Exception:
            pass

    # 3) ANSI fallback
    try:
        color_code = _ANSI_COLOR.get(color, '')
        attr_code = ''.join(_ATTRS_ANSI.get(a, '') for a in attrs)
        trailing = _RESET if (color_code or attr_code) else ''
        print(f"{color_code}{attr_code}{msg}{trailing}", file=file, end=end)
    except Exception:
        # 4) last resort: plain text
        try:
            print(msg, file=file, end=end)
        except Exception:
            # silence any error (we never want the logger itself to crash)
            pass


# -------------------------------------------------------------------------
# --- Convenient functions
# -------------------------------------------------------------------------
def print_bold(msg, **kwargs):
    cprint(msg, attrs=['bold'], **kwargs)

def FAIL(msg, label='[FAIL] ', **kwargs):
    msg = ('\n'+ ' ' * len(label)).join( (label+msg).split('\n') ) # Indending new lines
    cprint(msg, color='red', attrs=['bold'], file=sys.stderr, **kwargs)

def WARN(msg, label='[WARN] ', **kwargs):
    msg = ('\n'+ ' ' * len(label)).join( (label+msg).split('\n') ) # Indending new lines
    cprint(msg, color='yellow', attrs=['bold'], **kwargs)

def OK(msg, label='[ OK ] ', **kwargs):
    msg = ('\n'+ ' ' * len(label)).join( (label+msg).split('\n') ) # Indending new lines
    cprint(msg, color='green', attrs=['bold'], **kwargs)

def INFO(msg, label='[INFO] ', **kwargs):
    msg = ('\n'+ ' ' * len(label)).join( (label+msg).split('\n') ) # Indending new lines
    cprint(msg, **kwargs)


# --------------------------------------------------------------------------------
# --- Pretty prints
# --------------------------------------------------------------------------------
def pretty_num(x, digits=None, nchar=None, align='right', xmin=1e-16, center0=True):
    """ 
    Printing number with "pretty" formatting, either:
      - fixed number of decimals by setting digits
    OR
      - fixed number of characters by setting nchar

    """
    if nchar is not None and digits is not None:
        method='fixed_number_of_char_and_digits'

    elif nchar is None:
        nchar=7+digits
        method='fixed_number_of_digits'
    else:
        if digits is None:
            digits=int(nchar/2)
        method='fixed_number_of_char'
        if nchar<8:
            raise Exception('nchar needs to be at least 7 to accomodate exp notation')

    try:
        x = float(x)
    except:
        s=str(x)
        if align=='right':
            return s.rjust(nchar)
        else:
            return s.ljust(nchar)
    
    if np.abs(x)<xmin:
        x=0

    if x==0 and center0:
        s= ''.join([' ']*(nchar-digits-2))+ '0'+''.join([' ']*(digits+1))
    elif method=='fixed_number_of_digits':
        # --- Fixed number of digits
        if type(x)==int:
            s = f"{x:d}"
            #raise NotImplementedError()
        elif digits==6:
            if abs(x)<1000000 and abs(x)>1e-7:
                s= "{:.6f}".format(x)
            else:
               s= "{:.6e}".format(x)
        elif digits==5:
            if abs(x)<100000 and abs(x)>1e-6:
                s= "{:.5f}".format(x)
            else:
               s= "{:.5e}".format(x)
        elif digits==4:
            if abs(x)<10000 and abs(x)>1e-5:
                s= "{:.4f}".format(x)
            else:
               s= "{:.4e}".format(x)
        elif digits==3:
            if abs(x)<10000 and abs(x)>1e-4:
                s= "{:.3f}".format(x)
            else:
               s= "{:.3e}".format(x)
        elif digits==2:
            if abs(x)<100000 and abs(x)>1e-3:
                s= "{:.2f}".format(x)
            else:
               s= "{:.2e}".format(x)
        elif digits==1:
            if abs(x)<100000 and abs(x)>1e-2:
                s= "{:.1f}".format(x)
            else:
               s= "{:.1e}".format(x)
        elif digits==0:
            if abs(x)<1000000 and abs(x)>1e-1:
                s= "{:.0f}".format(x)
            else:
               s= "{:.0e}".format(x)
        else:
            raise NotImplementedError('digits',digits)
    elif method=='fixed_number_of_char':
        xlow  = 10**(-(nchar-2))
        xhigh = 10**( (nchar-1))
        if type(x)==int:
            raise NotImplementedError()
        if abs(x)<xhigh and abs(x)>xlow:
            n = int(np.log10(abs(x)))
            if n<0:
                sfmt='{:'+str(nchar)+'.'+str(nchar-3)+'f'+'}'
            elif nchar-3-n<0:
                sfmt='{:'+str(nchar-1)+'.0'+'f'+'}'
            elif nchar-3-n==0:
                sfmt='{:'+str(nchar-1)+'.0'+'f'+'}.'
            else:
                sfmt='{:'+str(nchar)+'.'+str(nchar-3-n)+'f'+'}'
        else:
            sfmt='{:'+str(nchar)+'.'+str(nchar-7)+'e'+'}' # Need 7 char for exp
        s = sfmt.format(x)
        #print(xlow, xhigh, sfmt, len(s), '>'+s+'<')
    elif method=='fixed_number_of_char_and_digits':
        xlow  = 10**(-(nchar-2))
        xhigh = 10**( (nchar-1))
        s = f"{x:.{digits+1}g}"  # general format with significant digits
        if len(s) > nchar:
            # fallback: scientific notation
            s = f"{x:.{digits+1}e}"
        # truncate or pad to exactly nchar characters
        if len(s) > nchar:
            s = s[:nchar]
    else:
        raise NotImplementedError(method)

    if align=='right':
        return s.rjust(nchar)
    else:
        return s.ljust(nchar)

def prettyMat(M, var=None, digits=2, nchar=None, sindent='   ', align='right', center0=True, newline=True, openChar='[',closeChar=']', sepChar=' ', xmin=1e-16):
    """ 
    return a matrix as a string, with misc output options
    INPUTS:
      - M: array of float/int
      - var: string
    """
    s=''
    if var is not None:
        if not isinstance(var, str):
            raise Exception()
        s=var+':'
        if newline:
            s+='\n'
    # Corner cases, being nice to user..
    if isinstance(M, str):
        s+=M
        return s
    if not hasattr(M,'__len__'):
        s+=pretty_num(M, digits=digits, nchar=nchar, align=align, center0=center0, xmin=xmin)
        return s

    M=np.atleast_2d(M)
    s+=sindent
    for iline,line in enumerate(M):
        s+= openChar+sepChar.join([pretty_num(v, digits=digits, nchar=nchar, align=align, center0=center0, xmin=xmin)  for v in line ])+closeChar
        if iline<M.shape[0]-1:
            s+='\n'+sindent
    return s

def printMat(M, var=None, **kwargs):
    """ 
    print a matrix with formatting options
    see prettyMat for input arguments and documentation
    
    example: 
        printMat(M, 'M')
          or 
        printMat('M', M)

        printMat(M, 'M', digits=1, align='right')
    """
    # Being nice if the user calls it by swapping the two arguments
    var, M = _swapArgs(var, M)
    M=np.asarray(M)
    print(prettyMat(M, var=var, **kwargs))

def printVec(M, var=None, newline=False, **kwargs):
    # Being nice if the user calls it by swapping the two arguments
    var, M = _swapArgs(var, M)
    M=np.asarray(M)
    M=np.atleast_2d(M)
    print(prettyMat(M, var=var, newline=newline, **kwargs))

def printDict(d, var=None, newline=False, digits=2, xmin=1e-16, **kwargs):
    s=''
    if var is not None:
        if not isinstance(var, str):
            raise Exception()
        s=var+':'
        if newline:
            s+='\n'
        print(s)
    sindent='  '
    for k,v in d.items():
        var='{:s}{:20s}'.format(sindent, k)
        if isinstance(v, str):
            print('{}:{}'.format(var, v))
        elif isinstance(v, int):
            print('{}:{:d}'.format(var, v))
        elif isinstance(v, float):
            print('{}:{}'.format(var, pretty_num(v, digits=digits, xmin=xmin, **kwargs)))
        elif isinstance(v, np.ndarray):
            if len(v.shape)==1:
                printVec(v, var, sindent=sindent, digits=digits, xmin=xmin, **kwargs)
            else:
                printMat(v, var, sindent=sindent+'   ', digits=digits, xmin=xmin, **kwargs)
        else:
            print('>>> printDict TYPE', type(v))
#             sindentloc = print('{}{20s}:{}'.format(sindent, k, v)


def prettyVar(val, var=None, key_fmt='{:15s}', digits=2, xmin=1e-16, **kwargs):
    s=''
    if var is not None:
        s+=key_fmt.format(var)+': '
    # Corner cases, being nice to user..
    if isinstance(val, str):
        s+=val
        return s
    if not hasattr(val,'__len__'):
        s+=pretty_num(val, digits=digits, **kwargs)
        return s

def printVar(val, var=None, key_fmt='{:15s}', digits=2, xmin=1e-16, **kwargs):
    var, val = _swapArgs(var, val)
    s = prettyVar(val, var=var, key_fmt=key_fmt, digits=digits, xmin=xmin, **kwargs)
    print(s)

def printVarTex(val, var=None, key_fmt='{:15s}', digits=2, xmin=1e-16, **kwargs):
    var, val = _swapArgs(var, val)
    s = ''
    if var is not None:
        var='$'+var+'$'
        s += key_fmt.format(var) + '&'   
    if isinstance(val, str):
        s+=val
    elif not hasattr(val,'__len__'):
        s+=pretty_num(val, digits=digits, **kwargs)
    s+='\\\\'
    return print(s)


def _swapArgs(var, val):
    if var is not None:
        if not isinstance(var, str):
            if isinstance(val, str):
                val, var = var, val # we swap 
    return var, val


# -------------------------------------------------------------------------
# Example usage
# -------------------------------------------------------------------------
if __name__ == "__main__":
    f= 10.**np.arange(-8,8,1)
    f1=10.**np.arange(-8,8,1)
    f2=-f1
    f3=f1*0
    M = np.stack((f1,f2,f3,f1))
    d=3
    nc=None
    d=None
    nc=12
    for x in f:
        print(pretty_num(x, digits=d, nchar=nc))
    for x in f:
        s=pretty_num(-x, digits=d, nchar=nc)
        print(s, len(s), -x)
    print(pretty_num(0, digits=d, nchar=nc))

    printMat(M, 'M', digits=1, align='right')


    FAIL("This is a failure message")
    WARN("This is a warning message")
    OK("This is a success message")
