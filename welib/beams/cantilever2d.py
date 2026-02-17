import numpy as np

def deflection(py, pz, r, beta, EI1, EI2):
    n = len(r)
    
    # Data sanitization
    r    = np.asarray(r)
    EI1  = np.asarray(EI1)
    EI2  = np.asarray(EI2)
    beta = np.asarray(beta)

    # Initialization includes boundary condition at the tip
    My = np.zeros(n)
    Mz = np.zeros(n)
    Ty = np.zeros(n)
    Tz = np.zeros(n)
    M1 = np.zeros(n)
    M2 = np.zeros(n)
    ky = np.zeros(n)
    kz = np.zeros(n)
    
    # Loop to get T and M from p
    for i in range(n-1, 0, -1):
        Ty[i-1] = Ty[i] + 0.5 * (py[i-1] + py[i]) * (r[i] - r[i-1])
        Tz[i-1] = Tz[i] + 0.5 * (pz[i-1] + pz[i]) * (r[i] - r[i-1])
        
        My[i-1] = My[i] - Tz[i] * (r[i] - r[i-1]) - (1/6 * pz[i-1] + 1/3 * pz[i]) * (r[i] - r[i-1])**2
        Mz[i-1] = Mz[i] + Ty[i] * (r[i] - r[i-1]) + (1/6 * py[i-1] + 1/3 * py[i]) * (r[i] - r[i-1])**2
    #print('Ty', Ty)
    #print('Tz', Tz)
    #print('My', My)
    #print('Mz', Mz)
    
    # Going from principal directions to y and z
    M1 = My * np.cos(beta) - Mz * np.sin(beta)
    M2 = My * np.sin(beta) + Mz * np.cos(beta)
    k1 = M1 / EI1
    k2 = M2 / EI2
    kz = -k1 * np.sin(beta) + k2 * np.cos(beta)
    ky =  k1 * np.cos(beta) + k2 * np.sin(beta)
    
    # Initialization includes boundary condition at the hub
    thetay = np.zeros(n)
    thetaz = np.zeros(n)
    uy = np.zeros(n)
    uz = np.zeros(n)
    
    for i in range(n-1):
        thetay[i+1] = thetay[i] + 0.5 * (ky[i+1] + ky[i]) * (r[i+1] - r[i])
        thetaz[i+1] = thetaz[i] + 0.5 * (kz[i+1] + kz[i]) * (r[i+1] - r[i])
        
        uy[i+1] = uy[i] + thetaz[i] * (r[i+1] - r[i]) + (1/6 * kz[i+1] + 1/3 * kz[i]) * (r[i+1] - r[i])**2
        uz[i+1] = uz[i] - thetay[i] * (r[i+1] - r[i]) - (1/6 * ky[i+1] + 1/3 * ky[i]) * (r[i+1] - r[i])**2
    
    return uy, uz
# , thetay, thetaz, ky, kz, Ty, Tz, My, Mz

def compute_modes(n_modes, r, EI1, EI2, m, beta, maxiter=500, abs_tol=0.01):
    """ """
    modes = []
    freqs = []
    n = len(r)

    # Data sanitization
    r    = np.asarray(r)
    EI1  = np.asarray(EI1)
    EI2  = np.asarray(EI2)
    m    = np.asarray(m)
    beta = np.asarray(beta)

    
    for mode in range(n_modes):
        py = np.ones(n)
        pz = np.ones(n)
        omega2_prev = 0
        
        # --- Iteratin loop
        for step in range(maxiter):
            uy, uz = deflection(py, pz, r, beta, EI1, EI2)
            
            # Orthogonalization against previous modes
            for prev_mode in modes:
                uy_p = prev_mode[0]
                uz_p = prev_mode[1]
                const = (np.trapezoid(uy_p * m * uy, r) + np.trapezoid(uz_p * m * uz, r)) / (np.trapezoid(uy_p * m * uy_p, r) + np.trapezoid(uz_p * m * uz_p, r))
                uy -=  const * uy_p
                uz -=  const * uz_p
            omega2 = pz[-1] / (uz[-1] * m[-1])
            norm_factor = np.sqrt(uz[-1]**2 + uy[-1]**2)
            pz = omega2 * m * uz / norm_factor
            py = omega2 * m * uy / norm_factor
            if abs(omega2 - omega2_prev) < abs_tol:
                break
            # Prepare for next iteration
            omega2_prev = omega2

        if step == maxiter - 1:
            raise Exception(f'Mode {mode+1} did not converge')
        
        modes.append((uy, uz))
        freqs.append(np.sqrt(omega2))
        
    return freqs, modes

def compute_modes_seq(r, EI1, EI2, m, beta, maxiter=500):
    """ 
    NOTE: this is a debugging funciton that performs the same as compute_modes but sequentially
    """
    modes = []
    freqs = []
    n = len(r)
    R = np.max(r)

    # --- 1st Flapwise Mode
    py = np.ones(n)
    pz = np.ones(n)
    uy, uz = deflection(py, pz, r, beta, EI1, EI2)
    omega2 = 0
    for step in range(maxiter):
        uy, uz = deflection(py, pz, r, beta, EI1, EI2)
        if abs(pz[-1] / (uz[-1] * m[-1]) - omega2) < 0.001:
            break
        omega2 = pz[-1] / (uz[-1] * m[-1])
        pz = omega2 * m * uz / np.sqrt(uz[-1]**2 + uy[-1]**2)
        py = omega2 * m * uy / np.sqrt(uz[-1]**2 + uy[-1]**2)
    if step==maxiter-1:
        raise Exception('Not converged')
    uy1f, uz1f, pz1f, py1f = uy, uz, pz, py
    modes.append((uy, uz))
    freqs.append(np.sqrt(omega2))

    # --- 1st Edgewise mode
    py = np.ones(n)
    pz = np.ones(n)
    omega2 = 0
    for step in range(maxiter):
        # Getting uz and uy
        uy, uz = deflection(py, pz, r, beta, EI1, EI2)
        const = (np.trapezoid(uz1f * m * uz  , r) + np.trapezoid(uy1f * m * uy  , r )) / \
                (np.trapezoid(uz1f * m * uz1f, r) + np.trapezoid(uy1f * m * uy1f, r ))
        uz1e = uz - const * uz1f
        uy1e = uy - const * uy1f
        verif = np.trapezoid(r, uz1f * m * uz1e) + np.trapezoid(r, uy1f * m * uy1e)
        if abs(pz[-1] / (uz1e[-1] * m[-1]) - omega2) < 0.01:
            break
        omega2 = pz[-1] / (uz1e[-1] * m[-1])
        pz = omega2 * m * uz1e / np.sqrt(uz1e[-1]**2 + uy1e[-1]**2)
        py = omega2 * m * uy1e / np.sqrt(uz1e[-1]**2 + uy1e[-1]**2)
    if step==maxiter-1:
        raise Exception('Not converged')
    modes.append((uy1e, uz1e))
    freqs.append(np.sqrt(omega2))


    # --- Second flap 
    py = np.ones(n)
    pz = np.ones(n)
    omega2 = 0
    for step in range(maxiter):
        uy, uz = deflection(py, pz, r, beta, EI1, EI2)
        const2 = (np.trapezoid(uz1e * m * uz, r)   + np.trapezoid(uy1e * m * uy, r)) / \
                 (np.trapezoid(uz1e * m * uz1e, r) + np.trapezoid(uy1e * m * uy1e, r))
        const1 = (np.trapezoid(uz1f * m * uz, r)   + np.trapezoid(uy1f * m * uy, r)) / \
                 (np.trapezoid(uz1f * m * uz1f, r) + np.trapezoid(uy1f * m * uy1f, r))
        uz2f = uz - const2 * uz1e - const1 * uz1f
        uy2f = uy - const2 * uy1e - const1 * uy1f
        verif = np.trapezoid(uz1e * m * uz2f, r) + np.trapezoid(uy1e * m * uy2f, r)
        if abs(pz[-1] / (uz2f[-1] * m[-1]) - omega2) < 0.01:
            break
        omega2 = pz[-1] / (uz2f[-1] * m[-1])
        norm_factor = np.sqrt(uz2f[-1] ** 2 + uy2f[-1] ** 2)
        pz = omega2 * m * uz2f / norm_factor
        py = omega2 * m * uy2f / norm_factor
    if step==maxiter-1:
        raise Exception('Not converged')
    modes.append((uy2f, uz2f))
    freqs.append(np.sqrt(omega2))


    # --- Second edge 
    py = np.ones(n)
    pz = np.ones(n)
    omega2 = 0
    for step in range(maxiter):
        uy, uz = deflection(py, pz, r, beta, EI1, EI2)
        const1 = (np.trapezoid(uz1f * m * uz, r)   + np.trapezoid(uy1f * m * uy, r)) / \
                 (np.trapezoid(uz1f * m * uz1f, r) + np.trapezoid(uy1f * m * uy1f, r))
        const2 = (np.trapezoid(uz1e * m * uz, r)   + np.trapezoid(uy1e * m * uy, r)) / \
                 (np.trapezoid(uz1e * m * uz1e, r) + np.trapezoid(uy1e * m * uy1e, r))
        const3 = (np.trapezoid(uz2f * m * uz, r)   + np.trapezoid(uy2f * m * uy, r)) / \
                 (np.trapezoid(uz2f * m * uz2f, r) + np.trapezoid(uy2f * m * uy2f, r))
        uz2e = uz - const3 * uz2f - const2 * uz1e - const1 * uz1f
        uy2e = uy - const3 * uy2f - const2 * uy1e - const1 * uy1f
        if abs(pz[-1] / (uz2e[-1] * m[-1]) - omega2) < 0.01:
            break
        omega2 = pz[-1] / (uz2e[-1] * m[-1])
        norm_factor = np.sqrt(uz2e[-1] ** 2 + uy2e[-1] ** 2)
        pz = omega2 * m * uz2e / norm_factor
        py = omega2 * m * uy2e / norm_factor
    if step==maxiter-1:
        raise Exception('Not converged')
    modes.append((uy2e, uz2e))
    freqs.append(np.sqrt(omega2))

    # --- Third flap 
    py = np.ones(n)
    pz = np.ones(n)
    omega2 = 0
    for step in range(1, maxiter + 1):
        uy, uz = deflection(py, pz, r, beta, EI1, EI2)
        
        const1 = (np.trapezoid(uz1f * m * uz, r)   + np.trapezoid(uy1f * m * uy, r)) / \
                 (np.trapezoid(uz1f * m * uz1f, r) + np.trapezoid(uy1f * m * uy1f, r))
        const2 = (np.trapezoid(uz1e * m * uz, r)   + np.trapezoid(uy1e * m * uy, r)) / \
                 (np.trapezoid(uz1e * m * uz1e, r) + np.trapezoid(uy1e * m * uy1e, r))
        const3 = (np.trapezoid(uz2f * m * uz, r)   + np.trapezoid(uy2f * m * uy, r)) / \
                 (np.trapezoid(uz2f * m * uz2f, r) + np.trapezoid(uy2f * m * uy2f, r))
        const4 = (np.trapezoid(uz2e * m * uz, r)   + np.trapezoid(uy2e * m * uy, r)) / \
                 (np.trapezoid(uz2e * m * uz2e, r) + np.trapezoid(uy2e * m * uy2e, r))
        uz3f = uz - const4 * uz2e - const3 * uz2f - const2 * uz1e - const1 * uz1f
        uy3f = uy - const4 * uy2e - const3 * uy2f - const2 * uy1e - const1 * uy1f
        if abs(pz[-1] / (uz3f[-1] * m[-1]) - omega2) < 0.01 or step == 100:
            break
        omega2 = pz[-1] / (uz3f[-1] * m[-1])
        norm_factor = np.sqrt(uz3f[-1] ** 2 + uy3f[-1] ** 2)
        pz = omega2 * m * uz3f / norm_factor
        py = omega2 * m * uy3f / norm_factor
    if step==maxiter-1:
        raise Exception('Not converged')
    modes.append((uy3f, uz3f))
    freqs.append(np.sqrt(omega2))

    return freqs, modes



if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from scipy.interpolate import CubicSpline
    from welib.tools.colors import python_colors

    # Data Initialization
    r0= np.array([0.30, 0.50, 1.00, 1.50, 2.00, 2.50, 2.60, 2.70, 2.80, 2.90, 2.95])
    EI1 = np.array([500.00, 468.75, 390.62, 312.50, 234.38, 167.19, 151.56, 137.50, 121.88, 107.81, 100.00]) * 100
    EI2 = np.array([12000.00, 11250.00, 9375.00, 7500.00, 5625.00, 4012.50, 3637.50, 3300.00, 2925.00, 2587.50, 2400.00]) * 100
    m = np.array([1.70, 1.59, 1.33, 1.06, 0.80, 0.57, 0.52, 0.47, 0.41, 0.37, 0.34])
    beta = np.array([21.32, 14.66, 7.69, 5.17, 3.31, 1.38, 0.90, 0.42, -0.07, 0.03, 0.74])*np.pi/180 # [rad]
    maxiter=500;

    # --- Reinterpolate
    r    = np.linspace(0.3, 3, 200)
    m    = CubicSpline(r0, m)(r)
    beta = CubicSpline(r0, beta)(r)
    EI1  = CubicSpline(r0, EI1)(r)
    EI2  = CubicSpline(r0, EI2)(r)

    # --- Constant property beam
#     n=200;
#     r    = np.linspace(0, 3, n)
#     m    = np.ones(n)*np.mean(m);
#     beta = np.zeros(n)
#     EI1  = np.ones(n)*np.mean(EI1)
#     EI2  = np.ones(n)*np.mean(EI2)

    n_modes = 5
    freqs, modes = compute_modes(n_modes, r, EI1, EI2, m, beta)
    freqs_ref, m_seq = compute_modes_seq(r, EI1, EI2, m, beta)

    # Display results
    for i, (f,m) in enumerate(zip(freqs,modes)):
        fr = freqs_ref[i]
        print(f'Mode {i+1}: Frequency = {f:.3f} - {fr:.3f}')
        plt.figure()
        plt.plot(r/3, m_seq[i][0], '-' , label='Uy', c=python_colors(0))
        plt.plot(r/3, m_seq[i][1], '-' , label='Uz', c=python_colors(1))
#         plt.plot(r, modes[i][0], '+', label=None, c=python_colors(0))
#         plt.plot(r, modes[i][1], '+', label=None, c=python_colors(1))
        plt.xlabel('x/L [-]')
        plt.ylabel('Deflection [-]')
        plt.legend()

    plt.show()

