""" 
Lumped LL vortex, with free wake vortex points

The circulation of the lifting line and near wake is computed
algebraically using the equations given by Katz Plotkin, section 13.7


"""

import numpy as np
import pandas as pd

def vp_u11(x, z, x1, z1, Gamma):
    """Influence of a point vortex at (x1, z1) with strength gamma."""
    rx = x - x1
    rz = z - z1
    r  = np.sqrt(rx**2 + rz**2)
    if r < 1e-6: return 0.0, 0.0
    v = 0.5/np.pi * Gamma / r
    u = v * (rz / r)
    w = v * (-rx / r)
    return u, w

def vps_u1(x, z, xv, zv, Gammas):
    """Induced velocity at (x, z) from arrays of vortex positions and strengths."""
    rx = x - xv
    rz = z - zv
    r2 = rx**2 + rz**2
    mask = r2 > 1e-12
    u = np.zeros_like(rx, dtype=float)
    w = np.zeros_like(rx, dtype=float)
    if np.any(mask):
        r = np.sqrt(r2[mask])
        v = 0.5/np.pi * Gammas[mask] / r
        u[mask] = v * (rz[mask] / r)
        w[mask] = v * (-rx[mask] / r)
    return np.sum(u), np.sum(w)

def ui_all(CPs, r_FW, Gamma_FW, r_LL, Gamma_LL):
    """Induced velocity at control points from free and lifting line vortices."""
    ncp = CPs.shape[0]
    u_CP = np.zeros((ncp, 2))
    for i in range(ncp):
        u, w   = vp_u11(CPs[i,0], CPs[i,1], r_LL[0], r_LL[1], Gamma_LL)
        u1, w1 = vps_u1(CPs[i,0], CPs[i,1], r_FW[:,0], r_FW[:,1], Gamma_FW[:])
        u_CP[i,0] = u + u1
        u_CP[i,1] = w + w1
    return u_CP

def wake_rollup(r_FW, Gamma_FW, dt, r_LL, Gamma_LL, it):
    """Update wake vortex positions using induced velocities and rollup."""
    CPs   = r_FW[:it+1]
    u_CP2 = ui_all(CPs, r_FW[:it+1], Gamma_FW[:it+1], r_LL, Gamma_LL)
    r_FW[:it+1] += u_CP2 * dt

def flat_plate_lin_acc(
    nstep=200, alpha_deg=5.0, dt_u0_over_c=0.25, U0=50.0, c=1.0, rho=1.0, dxw_factor=0.3
):
    ITMAX    = max(300, nstep+10)
    r_FW     = np.zeros((ITMAX, 2)) # Wake vortex positions
    Gamma_FW = np.zeros(ITMAX)      # Wake vortex strengths
    alpha    = np.deg2rad(alpha_deg)
    dt       = dt_u0_over_c * c / U0
    dxw      = dxw_factor * U0 * dt

    out  = []
    out.append([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

    Gamma_prev = 0.0
    t = -dt
    for it in range(nstep):
        if it > 0:
            # Propagate wake downstream
            r_FW[1:it+1]     = r_FW[0:it]
            Gamma_FW[1:it+1] = Gamma_FW[0:it]

        Gamma_LL = np.nan # To be determined at this time step

        t += dt
        sn = np.sin(alpha) # TODO replace with arbitrary motion
        cs = np.cos(alpha)
        r_LE = np.array([-U0 * t, 0.0]) # TODO replace with arbitrary motion
        r_LL = np.array([0.25*c*cs, -0.25*c*sn]) + r_LE # Lifting line
        r_CP = np.array([0.75*c*cs, -0.75*c*sn]) + r_LE # Control point (3/4 chord)

        # Shed new wake vortex at index 0 (latest)
        r_FW[0,:] = np.array([(c+dxw)*cs,  -(c+dxw)*sn]) + r_LE

        # --- Solve for Gamma_FW[0] and Gamma_LL
        # See Katz Plotkin, section 13.7
        # Coefficients for algebraic solution
        a = -1/(np.pi*c)
        b = 1/(2.0*np.pi*(c/4.0+dxw))
        rhs2 = 0.0
        wwake = 0.0
        # Known wake
        if it > 0:
            u, w = vps_u1(r_CP[0], r_CP[1], r_FW[1:it+1,0], r_FW[1:it+1,1], Gamma_FW[1:it+1])
            wwake = u*sn + w*cs
            rhs2 = -np.sum(Gamma_FW[1:it+1])

        rhs1 = -U0*sn - wwake
        Gamma_FW[0] = 1/(b/a - 1.0)*(rhs1/a - rhs2)
        Gamma_LL    = rhs2 - Gamma_FW[0]
        #print(Gamma_FW[0], Gamma_LL, rhs2, Gamma_FW[0]+ Gamma_LL- rhs2)

        # --- Wake rollup (update position of all vortices 0..it)
        wake_rollup(r_FW, Gamma_FW, dt, r_LL, Gamma_LL, it)

        # Aerodynamic loads
        que    = 0.5*rho*U0**2
        dGamdt = (Gamma_LL - Gamma_prev)/dt
        u, w   = vps_u1(r_CP[0], r_CP[1], r_FW[:it+1,0], r_FW[:it+1,1], Gamma_FW[:it+1])
        ww     = u*sn + w*cs
        L      = rho*(U0*Gamma_LL + dGamdt*c)
        D      = rho*(-ww*Gamma_LL + dGamdt*c*sn)
        Cl     = L/que/c
        Cd     = D/que/c
        Cl_rel = Cl/(2.0*np.pi*sn)
        Gamma_rel   = Gamma_LL/(np.pi*U0*c*sn)
        sx1    = r_LE[0] - U0*dt

        out.append([-sx1, Cd, Cl, Gamma_rel, Cl_rel, L])

        # Prepare for next iteration
        Gamma_prev = Gamma_LL

    out = np.array(out)
    df  = pd.DataFrame(out, columns=['SX', 'Cd', 'Cl', 'Gamma_rel', 'Cl_rel', 'L'])
    df.to_csv("_uLumpedLL.csv", index=False, sep=',')

    # Output wake in reversed order (latest vortex at index 0)
    wake_out = np.column_stack((r_FW[:nstep,0], r_FW[:nstep,1], Gamma_FW[:nstep]))
    header = "#X,Z,GAMMA"
    np.savetxt("_uLumpedLL_Wake.csv", wake_out, delimiter=",", fmt="%g", header=header, comments="")

    return df, wake_out

# -------------------------------------------------------------------------
import unittest
class TestFlatPlateAcc(unittest.TestCase):
    def test_out_values(self):
        df, wake_out = flat_plate_lin_acc(nstep=200, alpha_deg=5.0, dt_u0_over_c=0.25, U0=50.0, c=1.0, rho=1.0, dxw_factor=0.3)
        ref_clt   = [1.9696, 0.91240, 0.98932]
        ref_Gammat = [0.393939, 0.497633, 0.98909]
        idxs      = [1, 2, -1]
        np.testing.assert_almost_equal(df['Cl_rel'].to_numpy()[idxs],   ref_clt,   decimal=4)
        np.testing.assert_almost_equal(df['Gamma_rel'].to_numpy()[idxs], ref_Gammat, decimal=4)

if __name__ == "__main__":
    unittest.main()
    #df, wake_out = flat_plate_acc()

