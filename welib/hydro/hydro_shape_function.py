import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from scipy.signal import find_peaks
import dill

from welib.yams.section_loads import beamSectionLoads1D
from welib.essentials import *

def hydro_shape_function_from_TS(zDepth, vTime, p_hydro, eta=None, zRef = 0, plot=False, plotExtra=False):
    """ 
    Computes phi such that p_hydro(z,t) = phi(z) * eta_dot(t).
    """
    # Finding index of sea level, used as a reference
    bPos = zDepth<=0
    zPos = zDepth[bPos]
    nPos = len(zPos)
    p_wet = p_hydro[bPos, :]

    out=dict()


    # --- Method based on eta_dot
    # Solve p = phi_unscaled * eta_dot via Least Squares
    # phi = dot(p, eta_dot) / dot(eta_dot, eta_dot)
    if eta is not None:
        dt = vTime[1]-vTime[0]
        eta_dot = np.gradient(eta, dt) # Velocity state
        numerator   = np.sum(p_wet * eta_dot, axis=1)
        denominator = np.sum(eta_dot**2)
        phi_eta   = numerator / denominator
        p_recon     = phi_eta[:, np.newaxis] * eta_dot[np.newaxis, :]
        phit_eta  = phi_eta / np.max(phi_eta)
        out['phi_eta']  = phi_eta
        out['phit_eta'] = phit_eta

    # --- Methods based on a reference loading
    izRef = np.abs(zDepth-(zRef)).argmin()
    p_ref = p_hydro[izRef, :]

    # -- Instantaneous ratio compared to ref
    # Making sure the reference loading at sea level is more or less positive
    p_ref2 = np.abs(p_ref)
    p_mean = np.mean(p_ref2)
    p_ref2[np.abs(p_ref2)< 0.2*p_mean]=np.nan
    phit_rat_t    = np.abs(p_hydro[bPos, :]) / p_ref2  # Shape (nPos, nTime)
    phit_rat      = np.nanmean(phit_rat_t, axis=1)
    phit_rat_std  = np.nanstd (phit_rat_t, axis=1)
    out['phit_rat_t']   = phit_rat_t
    out['phit_rat_std'] = phit_rat_std
    out['phit_rat']     = phit_rat

    # -- RMS
    rms_at_depth = np.std(p_hydro[bPos, :], axis=1)
    phit_rms = rms_at_depth / rms_at_depth[izRef]
    out['phit_rms'] = phit_rms

    # -- Least square regression
    p_wet = p_hydro[bPos, :]
    numerator = np.nansum(p_wet * p_ref, axis=1)
    denominator = np.nansum(p_ref * p_ref)
    phit_ls = numerator / denominator
    out['phit_ls'] = phit_ls

    # -- Peaks
    peaks, _ = find_peaks(p_ref, height=0.1*np.max(p_ref), distance=10)
    p_at_peaks = p_hydro[bPos, :][:, peaks]
    p_ref_peaks = p_ref[peaks]
    phit_pk_t = p_at_peaks / p_ref_peaks  
    phit_pk      = np.mean(phit_pk_t, axis=1)
    phit_pk_std  = np.std (phit_pk_t, axis=1)
    out['phit_pk'] = phit_pk

    # TODO TODO TODO
    #phi_re-scaled = phit * phi_unscaled[iz0]

    fig = None
    if plotExtra:
        fig, ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8))
        fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        ax.plot(vTime, p_hydro[izRef, :], label='Reference Pressure (Sea Level)', color='gray', alpha=0.6)
        ax.plot(vTime[peaks], p_hydro[izRef, peaks], "x", color='red', label='Detected Peaks')

        ax.set_title(f"Peak Detection Debug (Found {len(peaks)} peaks)")
        ax.set_xlabel("Time [s]")
        ax.set_ylabel("Pressure [Pa]")

        if eta is not None:
            iz0 = np.abs(zDepth-(0)).argmin()
            izsb = np.abs(zDepth-np.min(zDepth)).argmin()
            fig, axes = plt.subplots(2, 1, sharey=False, figsize=(6.4,4.8))
            fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
            ax = axes[0]
            ax.plot(p_wet  [iz0, :], label='True p_hydro (z=0)')
            ax.plot(p_recon[iz0, :], '--', label='Reconstructed (phi * eta_dot)')
            ax = axes[1]
            ax.plot(p_wet  [izsb, :], label='True p_hydro (z=sea bed)')
            ax.plot(p_recon[izsb, :], '--', label='Reconstructed (phi * eta_dot)')
            ax.legend()
            fig.suptitle(f'Verification of Unscaled Shape Function')

    if plot:
        fig, ax = plt.subplots(1, 1, sharey=True, figsize=(6.4,4.8)) # (6.4,4.8)
        fig.subplots_adjust(left=0.14, right=0.95, top=0.95, bottom=0.12, hspace=0.20, wspace=0.20)
        cScat = (0.5, 0.5, 0.5)
#         for iz,z in enumerate(zPos):
#             ax.plot(phi_t[iz,:], [z]*len(vTime), '.', c=cScat, ms=4, alpha=0.5)
#             ax.plot(phi_mean[iz]+np.array([-phi_std[iz],+phi_std[iz]]), [z]*2, '-', c='k')
        ax.plot(phit_rat, zPos,  '-'   , label='Comp Ref - Instantaneous fit', lw=2)
        ax.plot(phit_pk,  zPos,   '.'  , label='Comp Ref - Peaks ratio', lw=2)
        ax.plot(phit_rms, zPos,  ':'   , label='Comp Ref - RMS ratio', lw=2)
        ax.plot(phit_ls,  zPos,   '-.' , label='Comp Ref - Least Square', lw=2)
        if eta is not None:
            ax.plot(phit_eta, zPos,   '-'  , label='Eta Least Square', lw=2)
        ax.tick_params(direction='in')
        ax.set_ylabel('Vertical position [m]')
        ax.set_xlabel('Hydro shape function [-]')

    if eta is None:
        print('[WARN] Eta not provided, unscaled shape function is wrong')
        phit = phit_rms
        phi  = phit*np.nan
    else:
        phit = phit_eta
        phi  = phi_eta

        return zPos, phit, phi, out, fig


def hydro_shape_function(rho, h, D, CM, z, kp, fp, ap, sum_method='RMS', returnPos=True):
    """
    Analytical unscaled shape function such that:

    p_hydro(z, t) = phi(z) * eta_dot(t)

    Also returns:
      - phit: shape function normalized to unity at the sea level
    
    """
    # Sanitization
    z = np.atleast_1d(z)
    ap = np.asarray(ap).flatten()
    kp = np.asarray(kp).flatten()
    fp = np.asarray(fp).flatten()

    omegap = 2 * np.pi * fp
    # Under water only, newaxis to allow for broadcasting with freq components 
    b = z <= 0
    zb = z[b][:, np.newaxis]
    iz0 = np.abs(z-(0)).argmin()

    # Robust handling of D and CM (scalar or array)
    def mask_geom(val):
        arr = np.atleast_1d(val)
        if len(arr) == 1: return arr[0] # Let broadcasting handle it
        return arr.flatten()[b][:, np.newaxis]
    D_b  = mask_geom(D)
    CM_b = mask_geom(CM)
    #D_b = np.atleast_1d(D).flatten()[b][:, np.newaxis]
    #CM_b = np.atleast_1d(CM).flatten()[b][:, np.newaxis]

    # Kinematic Kernel (Inertia term: force is proportional to omega)
    # Note: We divide by omegap once because eta_dot = omega * eta
    with np.errstate(over='ignore'):
        kernel = np.where(kp * h > 50, np.exp(kp * zb), np.cosh(kp * (zb + h)) / np.sinh(kp * h))

    # Single frequency 'alpha' (Force / eta_dot)
    # Units: [N/(m * (m/s))] or [kg/m^2] # <<<< WHATCH OUT FOR UNITS
    alpha_matrix = rho * (np.pi / 4) * (D_b**2) * CM_b * omegap * kernel
    # Local Force Matrix for RMS/Linear suming
    # Units: [N/m] # <<<< WHAT OUT FOR UNITS
    #F_matrix     = rho * (np.pi / 4) * (D_b**2) * CM_b * (omegap**2) * ap * kernel
    F_matrix = alpha_matrix * omegap * ap

    phi = np.zeros_like(z, dtype=float)

    # --- Always compute spectral as it gives good value at the surface
    phi_spec = np.zeros_like(z, dtype=float)
    weight = (omegap * ap)**2
    phi_spec[b] = np.sum(alpha_matrix * weight, axis=1) / np.sum(weight)
    phi_surface_ref = phi_spec[iz0] 

    if sum_method=='spectral':
        # Spectral weighting (LSE weighting)
        # To match the LSE 'phi' from a TS, we weight by the variance of eta_dot 
        # at each frequency: Var(eta_dot_i) ~ (omega_i * a_i)^2

        # High-frequency components of the JONSWAP spectrum have large force amplitudes at the sea level 
        # but decay to almost zero at sea bed
        # When performing the "Spectral Weighting," these high-frequency components pull the "average" shape at the seabed down. 
        # However, in reality, the sea-level pressure is dominated by the peak period ($T_p$), and those long waves penetrate much deeper.
        #weight = (omegap * ap)**2
        #phi[b] = np.sum(alpha_matrix * weight, axis=1) / np.sum(weight)
        phi = phi_spec

    elif sum_method == 'force_weighted':
        # Weighting by the square of the Force Amplitude (Power) 
        # This prioritizes frequencies that actually contribute to the total force
        # Units [N/(m * m/s)] because we divide by the same weight

        # In a JONSWAP spectrum, the majority of the "Load" energy is near the peak.
        # Mathematically: Weighting by $F^2$ effectively ignores the high-frequency "tail" of the spectrum that only exists at the surface, 
        # allowing the shape function to look more like the long-period waves that dominate the seabed.

        weight = (alpha_matrix * omegap * ap)**2 
        phi[b] = np.sum(alpha_matrix * weight, axis=1) / np.sum(weight)
        phit = phi/np.max(phi)
        phi = phit * phi_surface_ref

    elif sum_method == 'narrow_band':
        fp_peak = np.max(fp)
        # Only use frequencies within 10% of the peak frequency
        mask = (fp > 0.9*fp_peak) & (fp < 1.1*fp_peak)
        weight = (omegap[mask] * ap[mask])**2
        phi[b] = np.sum(alpha_matrix[:, mask] * weight, axis=1) / np.sum(weight)

        phit = phi/np.max(phi)
        phi = phit * phi_surface_ref


    elif sum_method=='RMS':
        # Result: Force Amplitude Profile
        phi[b] = np.sqrt(np.sum(F_matrix**2, axis=1))
        phit = phi/np.max(phi)
        phi = phit * phi_surface_ref

    elif sum_method=='linear':
        phi[b] = np.sum(F_matrix, axis=1)
        phit = phi/np.max(phi)
        phi = phit * phi_surface_ref


    else:
        raise NotImplementedError(sum_method)

    # Normalization (to Sea Level if possible, otherwise max)
    iz0 = np.argmin(np.abs(z - 0))
    ref_val = np.abs(phi[iz0]) if b[iz0] else np.max(np.abs(phi))
    phit = phi / ref_val if ref_val != 0 else phi

    if returnPos:
        return z[b], phit[b], phi[b]

    return z, phit, phi



def verify_hydro_reconstruction(zDepth, p_hydro, eta, phi_unscaled, vTime, tRange=None, plot=False, moment=False):
    """
    Verifies the accuracy of the reconstruction p(z,t) = phi(z) * eta_dot(t).
    """
    
    # Compute the reference state (eta_dot)
    dt = vTime[1] - vTime[0]
    eta_dot = np.gradient(eta, dt)
    
    # Identify underwater nodes
    bPos = zDepth <= 0
    zPos = zDepth[bPos]
    p_true = p_hydro[bPos, :]
    p_recon = phi_unscaled[:, np.newaxis] * eta_dot[np.newaxis, :]
    if moment :
        M_true  = np.zeros_like(p_true)
        M_recon = np.zeros_like(p_true)
        for it in range(len(vTime)):
            F_true,  M_true[:,it]  =  beamSectionLoads1D(zPos, p_true [:, it].flatten())
            F_recon, M_recon[:,it] =  beamSectionLoads1D(zPos, p_recon[:, it].flatten())
        p_true  = M_true
        p_recon = M_recon

    
    # Comprehensive Error Analysis
    # RMS and Standard Deviation
    std_true = np.std(p_true, axis=1)
    std_recon = np.std(p_recon, axis=1)
    mab_true  = np.mean(np.abs(p_true), axis=1)
    mab_recon = np.mean(np.abs(p_recon), axis=1)

    rms_true = np.sqrt(np.mean(p_true**2, axis=1))
    rms_error = np.sqrt(np.mean((p_true - p_recon)**2, axis=1))
    
    # Metrics
    rel_rms_err = (rms_error / rms_true) * 100
    mab_err     = (mab_recon - mab_true) / mab_true * 100
    mab_rat     = mab_true / mab_recon
    std_err     = (std_recon - std_true) / std_true * 100
    
    # Correlation Coefficient (Pearson R) for each depth
    corr = np.array([np.corrcoef(p_true[i, :], p_recon[i, :])[0, 1] for i in range(len(zPos))])
    # Since this is a linear reconstruction, R2 is simply the square of the correlation
    r2_score = corr**2

    # Select indices
    #indices = np.linspace(0, len(zPos)-1, 5, dtype=int)
    indices = np.linspace(0, len(zPos)-1, 5, dtype=int)
    iz0 = np.abs(zPos-(0)).argmin()
    izb = np.abs(zPos-(np.min(zPos))).argmin()
    for i, idx in enumerate(indices[::-1]):
        z_label = f'z={zPos[idx]:5.0f}m'
        print(f'{z_label} std={std_recon[idx]/1e6:7.2f} {std_true[idx]/1e6:7.2f} {std_err[idx]:6.2f}% - mean={mab_recon[idx]/1e6:7.2f} {mab_true[idx]/1e6:6.1f} {mab_err[idx]:6.2f}%', )

    fig=None
    if plot:
        #  Plotting
        fig = plt.figure(figsize=(14, 12))
        gs = fig.add_gridspec(3, 3) # Increased to 3 columns for better profile spacing
        
        # --- Subplot A: Time Series (Stacked Waterfall) ---
        ax0 = fig.add_subplot(gs[0:2, :])
        colors = plt.cm.viridis(np.linspace(0, 0.8, 5))
        if tRange is not None:
            t_mask = (vTime >= tRange[0]) & (vTime <= tRange[1])
        else:
            t_mask = vTime >= 0

        max_val = np.max(np.abs(p_true))
        #p0 = 10**np.floor(np.log10(max_val*2))
        p0 = max_val*1.5
        #p0 *=2
        #if max_val / p0 < 2: p0 /= 2 
        for i, idx in enumerate(indices):
            offset = i * p0
            z_label = f'z={zPos[idx]:5.0f}m'
            ax0.axhline(offset, color='gray', linestyle=':', lw=0.8, alpha=0.5)
            ax0.plot(vTime[t_mask], p_true[idx, t_mask] + offset, color=colors[i], alpha=0.5, lw=2.5, label=f'True {z_label}')
            ax0.plot(vTime[t_mask], p_recon[idx, t_mask] + offset, '--', color='black', lw=1, label='Recon' if i==0 else "")
            ax0.text(vTime[t_mask][0], offset + p0*0.05, z_label, color=colors[i], fontweight='bold', fontsize=9)

        ax0.set_title(r"Time Series Reconstruction: $p(z,t) \approx \phi(z) \cdot \dot{\eta}(t)$")
        ax0.set_yticks([])
        ax0.set_xlabel("Time [s]")
        ax0.legend(loc='upper right', ncol=3, fontsize='small', frameon=True)

        # --- Subplot B: Error Profiles (Vertical) ---
        ax1 = fig.add_subplot(gs[2, 0])
        ax1.plot(rel_rms_err, zPos, 'ro-', label='Rel. RMS Err %', ms=4)
        ax1.plot(std_err, zPos, 'bo--', label='Rel. Std Err %', ms=4)
        ax1.plot(mab_err, zPos, 'ko--', label='Rel. Mean Abs Err %', ms=4)
        ax1.axvline(0, color='k', lw=0.8)
        ax1.set_xlabel("Error [%]")
        ax1.set_ylabel("Depth [m]")
        ax1.set_title("Moment & RMS Errors")
        ax1.legend(fontsize='x-small')
        ax1.grid(True, alpha=0.3)

        # --- Subplot C: Correlation Profile ---
        ax2 = fig.add_subplot(gs[2, 1])
        ax2.plot(corr, zPos, 'go-', ms=4, label='Corr')
        ax2.plot(r2_score, zPos, 'o-', label='$R^2$', ms=4)
        #ax2.set_xlim([min(0.9, np.min(corr)), 1.01])
        ax2.set_xlabel("Correlation $R$ [-]")
        ax2.set_title("Signal Correlation Profile")
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        # --- Subplot D: Scatter at Sea Bed ---
        ax3 = fig.add_subplot(gs[2, 2])
        sc = 1e6
        ax3.scatter(p_true[izb, ::2]/sc, p_recon[izb, ::2]/sc, alpha=0.3, s=2, color='teal')
        lims = [np.min(p_true[izb,:]/sc), np.max(p_true[izb,:]/sc)]
        ax3.plot(lims, lims, 'k--', alpha=0.8)
        ax3.set_xlabel("F True [MN/m]")
        ax3.set_ylabel("F Recon [MN/m]")
        ax3.set_title(f"Seabed Correlation (z={zPos[izb]:.1f}m)")

        plt.tight_layout()
    return mab_rat, mab_err, fig



# --------------------------------------------------------------------------------}
# --- debug plot
# --------------------------------------------------------------------------------{
def verify_phase_lag(vTime, eta, p_hydro, zDepth, fs, plot=True):
    """
    Verifies phase lag between eta and p_hydro using Hilbert Transform 
    and Cross-Spectral Density.
    """
    from scipy.signal import hilbert, csd
    bPos = zDepth <= 0
    zPos = zDepth[bPos]
    
    # 1. Method: Analytic Signal (Hilbert) - Instantaneous Phase
    phase_eta = np.angle(hilbert(eta))
    lags_hilbert = []
    
    # 2. Method: CSD (Frequency Domain) - Phase at Peak Frequency
    lags_csd = []
    
    for iz in range(len(zPos)):
        p_z = p_hydro[iz, :]
        # Hilbert method average phase diff
        phase_p = np.angle(hilbert(p_z))
        diff = np.mean(np.unwrap(phase_p - phase_eta))
        lags_hilbert.append(diff)
        
        # CSD method
        f, Pxy = csd(eta, p_z, fs=fs, nperseg=1024)
        peak_idx = np.argmax(np.abs(Pxy))
        lags_csd.append(np.angle(Pxy[peak_idx]))

    lags_hilbert = np.array(lags_hilbert)
    lags_csd = np.array(lags_csd)

    if plot:
        fig, ax = plt.subplots(figsize=(6, 5))
        ax.plot(lags_hilbert*180/np.pi, zPos, 'o-', label='Hilbert (Time Domain)')
        ax.plot(lags_csd*180/np.pi, zPos, 'x--', label='CSD (Freq Domain)')
        ax.axvline(np.pi/2*180/np.pi, color='r', linestyle=':', label='Theoretical (pi/2)')
        ax.set_xlabel('Phase Lag [deg]')
        ax.set_ylabel('Depth [m]')
        ax.set_title('Phase Lag: Wave Elevation to Pressure')
        ax.legend()
        
    return zPos, lags_hilbert, lags_csd

# --------------------------------------------------------------------------------}
# --- Old Wrapper
# --------------------------------------------------------------------------------{
def hydro_shape_from_pickle(fstFile, IPlot=None, export=False, shapeFile=None):
    if IPlot is None:
        IPlot =[]

    matchSS = re.search(r"Hs=[\d\.]+_Tp=[\d\.]+", fstFile)
    SSCase = matchSS.group(0) if matchSS else None
    if shapeFile is None:
        shapeFile = f'HydroShapeFunction_{SSCase}.csv'

    with open(fstFile.replace('.fst','_prescribed.dpkl'), 'rb') as f:
        D = dill.load(f)
    vTime   = D['vTime']
    zDepth  = D['zDepth']
    # p_hydro = D['p_hydro']
    p_hydro = D['p_hnoac']
    F_sec   = D['F_sec']
    eta   = D['eta_sim']
    pSS   = D['pSS']
    pST   = D['pST']
    pHD   = D['pHD']
    dt=vTime[1]-vTime[0]
    fs=1/dt

    print('ST keys',pST.keys())

    # Numerical shape functions
    zPos, phit_eta, phi_eta, out, fig = hydro_shape_function_from_TS(zDepth, vTime, p_hydro, eta=eta, zRef=0, plot=False, plotExtra=False)
    # Analytical shape function
    zPos, phit_rms, phi_rms = hydro_shape_function(pSS['rho'], h=pSS['WaterDepth'], D=pHD['D'], CM=pHD['CM'], z=pHD['zDepth'], kp=pSS['kp'], fp=pSS['fp'], ap=pSS['ap'], sum_method='RMS')

    # Calibrating numerical shape function
    mab_rat, _, _ = verify_hydro_reconstruction(zDepth, p_hydro, eta, phi_eta, vTime, plot=False, moment=False)
    phi_cal = phi_eta * mab_rat

    if export:
        phi2 = np.zeros_like(zDepth)
        phi2[:len(zPos)] = phi_cal
        phit2 = phi2/np.max(phi2)
        # Units [N/(m * m/s)] because we divide by the same weight
        df = pd.DataFrame(data=np.column_stack((zDepth, phi2, phit2)), columns=['z_[m]', 'phi_[Ns/m^2]', 'phit_[-]'])
        df.to_csv(shapeFile, index=False)
        print('>>> Written: ', shapeFile)


    # --------------------------------------------------------------------------------}
    # ---  
    # --------------------------------------------------------------------------------{
    if 1 in IPlot:
        fig, ax = plt.subplots(1, 1, sharey=True, figsize=(6.4,3.8)) # (6.4,4.8)
        fig.subplots_adjust(left=0.12, right=0.92, top=0.97, bottom=0.13, hspace=0.20, wspace=0.20)
        cScat = (0.5, 0.5, 0.5)
        zPos=(zPos-np.min(zPos))/(np.max(zPos)-np.min(zPos))
        ax.plot([np.nan], [np.nan], 'o', c=cScat, ms=4, alpha=0.5, label='Instantaneous ratio')
        scatter_x=[]
        scatter_y=[]
        for iz,z in enumerate(zPos):
            scatter_x .append(out['phit_rat_t'][iz,:])
            scatter_y .append([z]*len(vTime))
            #ax.plot(out['phit_rat_t'][iz,:], [z]*len(vTime), 'o', c=cScat, ms=3, alpha=0.05)
            #ax.plot(out['phit_rat'][iz]+np.array([-out['phit_rat_std'][iz],+out['phit_rat_std'][iz]]), [z]*2, '|', c='k', label='Std of fit' if iz==0 else None)
    #     ax.plot(out['phit_rat'], zPos,  '-'   , label='Comp Ref - Instantaneous ratio', lw=2)
        scatter_x = np.asarray(scatter_x).ravel()
        scatter_y = np.asarray(scatter_y).ravel()
        ax.scatter(scatter_x, scatter_y, c=cScat, s=9, alpha=0.05, rasterized=True, edgecolors='none')
        for iz,z in enumerate(zPos):
            ax.plot(out['phit_rat'][iz]+np.array([-out['phit_rat_std'][iz],+out['phit_rat_std'][iz]]), [z]*2, '|', c='k', label='Std. of inst. ratio' if iz==0 else None)

        ax.plot(phit_eta, zPos, 'k',  label='Least square fit')
        ax.plot(phit_rms, zPos, '--', label='Analytical (RMS)', lw=2, c=fColrs(1))
        ax.tick_params(direction='in')
        ax.set_ylabel('Vertical position [-]')
        ax.set_xlabel(r'Hydrodynamic shape function, $\phi_h$ [-]')
        ax.set_xlim([0,2])
        ax.set_xticks([0,0.5, 1, 1.5, 2])
        ax.set_ylim([0,1])
        ax.tick_params(direction='in', top=True, right=True)
        ax.legend(fontsize=12)
        fig._title='HydroShapeFunction_'+SSCase


    # --------------------------------------------------------------------------------}
    # --- DEBUG 
    # --------------------------------------------------------------------------------{
    if 0 in IPlot:
        # Analytical shape function (other)
        zPos, phit_lin, phi_lin = hydro_shape_function(pSS['rho'], h=pSS['WaterDepth'], D=pHD['D'], CM=pHD['CM'], z=pHD['zDepth'], kp=pSS['kp'], fp=pSS['fp'], ap=pSS['ap'], sum_method='linear')
        zPos, phit_spc, phi_spc = hydro_shape_function(pSS['rho'], h=pSS['WaterDepth'], D=pHD['D'], CM=pHD['CM'], z=pHD['zDepth'], kp=pSS['kp'], fp=pSS['fp'], ap=pSS['ap'], sum_method='spectral')
        zPos, phit_fwt, phi_fwt = hydro_shape_function(pSS['rho'], h=pSS['WaterDepth'], D=pHD['D'], CM=pHD['CM'], z=pHD['zDepth'], kp=pSS['kp'], fp=pSS['fp'], ap=pSS['ap'], sum_method='force_weighted')
        zPos, phit_nbd, phi_nbd = hydro_shape_function(pSS['rho'], h=pSS['WaterDepth'], D=pHD['D'], CM=pHD['CM'], z=pHD['zDepth'], kp=pSS['kp'], fp=pSS['fp'], ap=pSS['ap'], sum_method='narrow_band')

        fig, ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8))
        fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        ax.plot(phi_eta, zPos, label='Time series, from Eta')
        ax.plot(phi_spc, zPos, label='Analytical (spec)')
        ax.plot(phi_fwt, zPos, label='Analytical (fwt)')
        ax.plot(phi_nbd, zPos, label='Analytical (nbd)')
        ax.plot(phi_lin, zPos, label='Analytical (lin)')
        ax.plot(phi_rms, zPos, label='Analytical (rms)')
        ax.plot(phi_cal, zPos, label='Calibrated')
        ax.set_xlabel('')
        ax.set_ylabel('')
        fig.suptitle('Unscaled shape function')
        ax.legend()


        print('>>>>> cal')
        verify_hydro_reconstruction(zDepth, p_hydro, eta, phi_cal, vTime, plot=True)
        # print('>>>>> spc')
        # verify_hydro_reconstruction(zDepth, p_hydro, eta, phi_spc, vTime, plot=False)
        # print('>>>>> RMS')
        # verify_hydro_reconstruction(zDepth, p_hydro, eta, phi_rms, vTime, plot=False)


        print('-------------------------------- MOMENT')
        print('>>>>> cal')
        verify_hydro_reconstruction(zDepth, p_hydro, eta, phi_cal, vTime, plot=True, moment=True)
        print('>>>>> spc')
        verify_hydro_reconstruction(zDepth, p_hydro, eta, phi_spc, vTime, plot=False, moment=True)
        print('>>>>> RMS')
        verify_hydro_reconstruction(zDepth, p_hydro, eta, phi_rms, vTime, plot=False, moment=True)


    # --------------------------------------------------------------------------------}
    # ---  Misc debug plots
    # --------------------------------------------------------------------------------{
    # verify_phase_lag(vTime, eta, p_hydro, zDepth, fs, plot=True)
    # plot_Depth_force_ratio(vTime, zDepth, p_hydro, eta, DEPTH=DEPTH)
    return zPos, phit_eta, phi_eta, phit_rms, phi_rms, phi_cal
