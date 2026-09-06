import numpy as np
from welib.tools.strings import prettyMat

def pretty_PrintMat(M,fmt='{:11.3e}',fmt_int='    {:4d}   ',sindent='   '):
    s = prettyMat(M, var=None, digits=2, nchar=None, sindent='   ', align='right', center0=True, newline=True, openChar='[',closeChar=']', sepChar=' ', xmin=1e-16)
    return s

class AugmentedLinModel():
    """
    Prepare the state matrices to be given to a Kalman Filter algorithm

    Open up a linear physical model, convert it to an augmented system

    Main Data:
        self.sQ  list of states            
        self.sQa list of augmented states  
        self.sY  list of measurements      
        self.sU  list of inputs            
        self.sS  list of additional storage
        self.sQd = None
        self.A   = None
        self.B   = None
        self.C   = None
        self.D   = None
    """
    def __init__(self):
        self.sQ  = None
        self.sQa = None
        self.sY  = None
        self.sU  = None
        self.sS  = None
        self.sQd = None
        self.A   = None
        self.B   = None
        self.C   = None
        self.D   = None

    def __repr__(self):
        s='<{} object> with attributes:\n'.format(type(self).__name__)
        s+='  sQ  : {} \n'.format(self.sQ)
        s+='  sQa : {} \n'.format(self.sQa)
        s+='  sY  : {} \n'.format(self.sY)
        s+='  sU  : {} \n'.format(self.sU)
        s+='  sS  : {} \n'.format(self.sS)
        if self.A is not None:
            s+=' A: State-State Matrix\n'
            s+=pretty_PrintMat(self.A)+'\n'
        if self.B is not None:
            s+=' B: State-Input Matrix\n'
            s+=pretty_PrintMat(self.B)+'\n'
        if self.C is not None:
            s+=' C: Output-State Matrix\n'
            s+=pretty_PrintMat(self.C)+'\n'
        if self.D is not None:
            s+=' D: Output-Input Matrix \n'
            s+=pretty_PrintMat(self.D)+'\n'
        return s

    @property
    def sX(self): return np.concatenate((self.sQ, self.sQa))
    @property
    def iX(self): return {lab: i   for i,lab in enumerate(self.sX)}
    @property
    def iY(self): return {lab: i   for i,lab in enumerate(self.sY)}
    @property
    def iU(self): return {lab: i   for i,lab in enumerate(self.sU)}
    @property
    def iS(self): return {lab: i   for i,lab in enumerate(self.sS)}
