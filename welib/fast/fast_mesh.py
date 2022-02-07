import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local 
class PointMesh:
    def __init__(self, nPoints):
        # NOTE: for now, using Fortran convention for memory order
        self.Position         = np.zeros((3, nPoints))
        self.TranslationDisp  = np.zeros((3, nPoints))
        self.TranslationVel   = np.zeros((3, nPoints))
        self.TranslationAcc   = np.zeros((3, nPoints))
        self.RefOrientation = np.zeros((3,3, nPoints))
        self.Orientation    = np.zeros((3,3, nPoints))
        self.RotationVel = np.zeros((3, nPoints))
        self.RotationAcc = np.zeros((3, nPoints))
        self.Force  = np.zeros((3, nPoints))
        self.Moment = np.zeros((3, nPoints))
        self.nNodes = nPoints

    def mapLoadsToPoint(self,P):
        P = np.asarray(P)
        F = np.zeros(3)
        M = np.zeros(3)
        for j in range(self.nNodes):
            F0 = self.Force[:,j]
            M0 = self.Moment[:,j]
            P0 = self.Position[:,j] + self.TranslationDisp[:,j]
            r = P0-P
            dM = np.cross(r, F0)
            F+=F0
            M+=M0 + dM
            #M+=M0 
        return F, M

