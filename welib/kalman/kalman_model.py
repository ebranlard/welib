import numpy as np

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
