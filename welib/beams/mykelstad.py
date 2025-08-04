""" Mykelstad method to determine natural frequency of a cantilever
beam

TODO TODO TODO NEED TO BE REVISITED

"""

import numpy as np
import matplotlib.pyplot as plt
pi = np.pi
R_top = 3                                       # radius at top, m
t_top=.02                                       # wall thickness at top, m
R_base = 4                                      # radius at base, m
t_base =.03                                     # wall thickness at base, m
R_top_in = R_top-t_top                          # inner radius at top, m
R_base_in = R_base-t_base                       # inner radius at base, m
E = 210e9                                       # modulus of elasticity, Pa
L_overall = 80                                  # length, m
rho = 8500                                      # density of material, kg/m^3
Nsection = 100                                 # number of stations to use

# Initialize variables below
EI = np.zeros(Nsection)                        # section stiffness, N m^2
m = np.zeros(Nsection)                         # section mass, kg 
delta_x = L_overall/(Nsection)                 # dimensionless distance between stations
aLength = [delta_x] * Nsection
x = np.linspace(0, delta_x * (Nsection - 1) / L_overall, Nsection)

for i in range(Nsection):

    R_out = R_top  + x[i] * (R_base - R_top) 
    R_in = R_out - t_top - x[i] * (t_base - t_top)
    
    # below are used in Myklestad function
    m[i] = rho*delta_x*pi*(R_out**2 - R_in**2)
    EI[i] = E*pi*(R_out**4 - R_in**4)/4
aMass = m
 
def Influence(NSection):
# This is used to obtain the influence coefficients from the length and stiffness of each section

# Inputs 
# NSection = number of sections

# Inputs from main program
# aLength[] = length of sections
# EI[]  =  stiffness of sections

# Outputs returned
# SlopeFromMoment[]
# SlopeFromShear[]
# DeflecFromShear[]

    for i in range(NSection):
        SlopeFromMoment = aLength / EI
        SlopeFromShear = SlopeFromMoment * aLength / 2
        DeflecFromShear = SlopeFromShear * aLength * 2 / 3

    return SlopeFromMoment, SlopeFromShear, DeflecFromShear
    
def Myklestad_(NSections,f_start, f_final, delta_f):
# This function performs the calculations to find the natural
# frequencies of a non-uniform cantileverd beam

# Inputs
# NSections = number of sections
# aLength = length of section
# EI = stiffness of section
# aMass = mass of section
# f_start = starting frequency in calculations, Hz
# f_final = ending frequency in calculations, Hz
# delta_f = frequency step, Hz

# Outputs
# mode = number of vibration modes between omegaStart and omegaFinal
# freq = natural frequency of all modes
    modes = [] #added
    modesY=[]
    freq = []
    omegaStart, omegaFinal, deltaOmega  = f_start*2*pi, f_final*2*pi, delta_f*2*pi
    YInitial_1, YInitial_2 = 1 , 0 # intial displacement at free end and fixed end respectivelly
    ThetaInitial_1, ThetaInitial_2  = 0, 1
    mode = 0  
    SlopeFromMoment = np.zeros(NSections)
    SlopeFromShear = np.zeros(NSections)
    DeflecFromShear = np.zeros(NSections)
    """
    for rotating beam (not used in this version)
    # xFromAxisOfRotation = np.zeros(NSections)
    # FCent_1 = np.zeros(NSections)
    # FCent_2 = np.zeros(NSections)
    """
    ThetaSlope_1 = np.zeros(NSections)
    ThetaSlope_2 = np.zeros(NSections)
    Ydeflec_1 = np.zeros(NSections)
    Ydeflec_2 = np.zeros(NSections)
    VShear_1 = np.zeros(NSections)
    VShear_2 = np.zeros(NSections)
    MBend_1 = np.zeros(NSections)
    MBend_2 = np.zeros(NSections)
    deflection_full = np.zeros(NSections)
    # Get influence coefficients
    temp = Influence(NSections)

    SlopeFromMoment = temp[0]
    SlopeFromShear = temp[1]
    DeflecFromShear = temp[2]
    
    """
    for rotating beam (not used in this version)
    # xFromAxisOfRotation[NSections - 1] = aLength[1]
    # xFromAxisOfRotation[NSections - 1] = aLength[1] / 2
    # for i in range(NSections):
        # xFromAxisOfRotation[i] = xFromAxisOfRotation[i + 1] + aLength[i]
    """
    
    #omega1 = omegaStart-deltaOmega
    deflec1 = 0
    nsteps = int((omegaFinal-omegaStart)/deltaOmega)
    omega_range = np.linspace(omegaStart, omegaStart + deltaOmega * (nsteps -1), nsteps) #causes slight change after 13-14 decimal places, but much cleaner

    for omega1 in omega_range: #  = omegaStart To omegaFinal Step deltaOmega
    # while deflec1==0 or (abs(deflec) >0 and np.sign(deflec1)==np.sign(deflec)):
    # Note! this starts at the free end of the beam 
            """
            for rotating beam (not used in this version)
            # FCent_1[1] = 0
            # FCent_2[1] = 0
            """
            ThetaSlope_1[0], ThetaSlope_2[0] = ThetaInitial_1, ThetaInitial_2
            Ydeflec_1[0], Ydeflec_2[0] = YInitial_1, YInitial_2
            VShear_1[0], VShear_2[0], MBend_1[0], MBend_2[0]  = 0, 0, 0, 0  #adjust first parameter to simulate tip load

            for i in range(1,NSections):
                
                """
    for rotating beam (not used in this version)
    ------------
                    FCent[i, j] = FCent[i - 1, j] + omegaRotate ** 2 * aMass[i - 1] * xFromAxisOfRotation[i - 1]
                    cRotate1[j] = 1 - FCent[i, j] * aLength[i - 1] ** 2 / [2 * EI[i - 1]]
                    cRotate2[j] = FCent[i, j] * aLength[i - 1] ** 3 / [3 * EI[i - 1]]
                    # cRotate1[j] = 1 # test
   # ------------
                """

                VShear_1[i] = VShear_1[i - 1] - aMass[i - 1] * omega1 ** 2 * Ydeflec_1[i - 1] # - FCent_1[i] * ThetaSlope_1[i - 1]
                VShear_2[i] = VShear_2[i - 1] - aMass[i - 1] * omega1 ** 2 * Ydeflec_2[i - 1] # - FCent_2[i] * ThetaSlope_2[i - 1]
                    
                MBend_1[i] = MBend_1[i - 1] - VShear_1[i] * aLength[i - 1]    # - cRotate2_1[j]) + FCent_1[i, j] * ThetaSlope_1[i - 1, j] * aLength[i - 1]) / cRotate1_1[j]
                MBend_2[i] = MBend_2[i - 1] - VShear_2[i] * aLength[i - 1]    # - cRotate2_2[j]) + FCent_2[i, j] * ThetaSlope_2[i - 1, j] * aLength[i - 1]) / cRotate1_2[j]                                                       
                    
                ThetaSlope_1[i] = ThetaSlope_1[i - 1] + MBend_1[i] * SlopeFromMoment[i - 1] + VShear_1[i] * SlopeFromShear[i - 1]
                ThetaSlope_2[i] = ThetaSlope_2[i - 1] + MBend_2[i] * SlopeFromMoment[i - 1] + VShear_2[i] * SlopeFromShear[i - 1]
                    
                Ydeflec_1[i] = Ydeflec_1[i - 1] + ThetaSlope_1[i - 1] * aLength[i - 1] + MBend_1[i] * SlopeFromShear[i - 1] + VShear_1[i] * DeflecFromShear[i - 1]
                Ydeflec_2[i] = Ydeflec_2[i - 1] + ThetaSlope_2[i - 1] * aLength[i - 1] + MBend_2[i] * SlopeFromShear[i - 1] + VShear_2[i] * DeflecFromShear[i - 1]
                deflection_full[i] = Ydeflec_1[i] + Ydeflec_2[i] * -ThetaSlope_1[i]/ThetaSlope_2[i]
            slope = -ThetaSlope_1[NSections-1] / ThetaSlope_2[NSections-1]       # slope
            deflec = Ydeflec_1[NSections-1] + Ydeflec_2[NSections-1] * slope     # deflection (should ~ 0)
            
            if deflec1!= 0 and np.sign(deflec1)!= np.sign(deflec):
                freq.append(omega1/(2*pi))
                mode = mode+1
                modes.append(deflec)    #added
                modesY.append(deflection_full.copy())
            deflec1 = deflec     
            
    return freq,mode, modes, modesY

if __name__ == '__main__':
    freq, mode, modes, modesY= Myklestad_(100,.1, 25, .001)
# fn = np.zeros(10)
# n = freq[1]
# fn = freq[0]
    print("")
    print("Myklestad method")
    print(modes)
    print(modesY)
    print(freq)
    print(x)
    modesY = np.array(modesY)

    # # Tip loading condition ...
    # deflec_analytical = []
    # for i in range(Nsection-1):
    #      deflec_analytical.append(-100*(x[i]*L_overall)**2/ (6 * EI[i]) * (3*L_overall - x[i]*L_overall))

    plt.figure(figsize=(8, 5))
    plt.plot(x[1:], -modesY[0][1:], label="Deflection Mode Shape 1")
    plt.plot(x[1:], -modesY[1][1:], label="Deflection Shape Mode Shape 2")
    plt.plot(x[1:], -modesY[2][1:], label="Deflection Shape Mode Shape 3")
    #plt.plot(x[1:], deflec_analytical)
    plt.legend()
    plt.show()
