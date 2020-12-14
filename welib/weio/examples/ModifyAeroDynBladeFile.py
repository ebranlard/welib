import weio
import numpy as np
Bld=weio.read('../../../data/NREL5MW/data/NREL5MW_AD15_blade.dat')
nSpan = 10

Spn   = np.linspace(0, 15, nSpan)       # BlSpn, radial stations [m]
CrvAC = np.zeros((nSpan,))              # BlCrvAC, prebend (usually <0) [m]
SwpAC = np.zeros((nSpan,))              # BlSwpC,  sweep                [m]
CrvAng = np.concatenate(([0], np.arctan2((CrvAC[1:]-CrvAC[:-1]),(Spn[1:]-Spn[:-1]))*180/np.pi))
Twist  = np.zeros((nSpan,)) + 1         # BlTwist [deg]
Chord  = np.zeros((nSpan,)) + 5         # BlChord [m]
AFID   = np.zeros((nSpan,)).astype(int) # BlAFID [-]
ADProp = np.column_stack((Spn,CrvAC,SwpAC,CrvAng,Twist,Chord,AFID))
Bld['NumBlNds']     = ADProp.shape[0]
Bld['BldAeroNodes'] = ADProp


Bld.write('_AeroDyn_Blade_Modified.dat')
