import numpy as np

hline='\\hline\n'

def tolatex(M,linenames=None,fmt='{:5.2f}'):
    s=''

    if linenames is not None:
        slen = np.max([len(s) for s in linenames])
        sfmt = '{:'+str(slen)+'.'+str(slen)+'s}'

    for iline,line in enumerate(M):
        if linenames is not None:
            s+=sfmt.format(linenames[iline])+' & '
        s+=' & '.join([fmt.format(v) for v in line ])
        s+='\\\\\n'
    #s='\\\\\n'.join(['&'.join(["{:5.1f}".format(v) for v in line ]) for line in M]) 
    return s

if __name__=='__main__':
    import numpy as np
    M=np.zeros((2,3))
    M[0,:]=1
    M[1,:]=2
    print(tolatex(M, linenames=['Sim1','Sim2']))
