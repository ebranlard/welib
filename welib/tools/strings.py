import numpy as np

def FAIL(msg):
    try: 
        from termcolor import cprint
        cprint('[FAIL] ' + msg , 'red', attrs=['bold'], file=sys.stderr)
    except:
        HEADER = '\033[95m'
        RED = '\033[91m'
        ENDC = '\033[0m'
        BOLD = '\033[1m'
        UNDERLINE = '\033[4m'
        print(RED+'[FAIL] ' + msg + ENDC)

def WARN(msg):
    try: 
        from termcolor import cprint
        cprint('[WARN] ' + msg , color='yellow', attrs=['bold'])
    except:
        ORAN = '\033[93m'
        ENDC = '\033[0m'
        print(ORAN+'[WARN] ' + msg + ENDC)

def OK(msg):
    try: 
        from termcolor import cprint
        cprint('[ OK ] ' + msg , 'green', attrs=['bold'])
    except:
        GREEN = '\033[92m'
        ENDC = '\033[0m'
        print(GREEN+'[ OK ] ' + msg + ENDC)

def pretty_num(x, digits=None, nchar=None, align='right', xmin=1e-16, center0=True):
    """ 
    Printing number with "pretty" formatting, either:
      - fixed number of decimals by setting digits
    OR
      - fixed number of characters by setting nchar

    """
    if np.abs(x)<xmin:
        x=0

    if nchar is None:
        nchar=7+digits
        method='fixed_number_of_digits'
    else:
        digits=int(nchar/2)
        method='fixed_number_of_char'
        if nchar<8:
            raise Exception('nchar needs to be at least 7 to accomodate exp notation')


    if x==0 and center0:
        s= ''.join([' ']*(nchar-digits-2))+ '0'+''.join([' ']*(digits+1))
    elif method=='fixed_number_of_digits':
        # --- Fixed number of digits
        if type(x)==int:
            raise NotImplementedError()
        if digits==6:
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
#         print(xlow, xhigh, sfmt, len(s))
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
    if not isinstance(var, str):
        if isinstance(M, str):
            # we swap 
            M, var = var, M
    M=np.asarray(M)
    print(prettyMat(M, var=var, **kwargs))

def printVec(M, var=None, newline=False, **kwargs):
    # Being nice if the user calls it by swapping the two arguments
    if var is not None:
        if not isinstance(var, str):
            if isinstance(M, str):
                # we swap 
                M, var = var, M
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


if __name__ == '__main__':
    f= 10.**np.arange(-8,8,1)
    f1=10.**np.arange(-8,8,1)
    f2=-f1
    f3=f1*0
    M = np.stack((f1,f2,f3,f1))
    d=3
    nc=None
#     d=None
#     nc=12
#     for x in f:
#         print(pretty_num(x, digits=d, nchar=nc))
#     for x in f:
#         s=pretty_num(-x, digits=d, nchar=nc)
#         print(s, len(s), -x)
#     print(pretty_num(0, digits=d, nchar=nc))

    printMat(M, 'M', digits=1, align='right')
