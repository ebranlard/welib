import numpy as np

class AugmentedLinModel():
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
