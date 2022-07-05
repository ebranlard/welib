import numpy as np

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
        if digits==4:
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

    if align=='right':
        return s.rjust(nchar)
    else:
        return s.ljust(nchar)

def prettyMat(M, var=None, digits=2, nchar=None, sindent='   ', align='right', center0=True, newline=True):
    s=''
    if var is not None:
        s=var+':'
        if newline:
            s+='\n'
    s+=sindent
    for iline,line in enumerate(M):
        s+='['+' '.join([pretty_num(v, digits=digits, nchar=nchar, align=align, center0=center0)  for v in line ])+']'
        if iline<M.shape[0]-1:
            s+='\n'+sindent
    return s

def printMat(*args, **kwargs):
    print(prettyMat(*args, **kwargs))

def printVec(M, *args, newline=False, **kwargs):
    M=np.atleast_2d(M)
    print(prettyMat(M, *args, newline=newline, **kwargs))

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

    print(printMat(M, 'M', digits=1, align='right'))
