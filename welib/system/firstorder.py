r"""
Set of tools/functions for a first order system expressed in scalar form. 
For multiple dimensions, see statespace, linearstatespace

 - With constant coefficients, the characteristic equation is:
     
       x_dot(t) + 1/tau x = f(t)                                                    (1)

   Generic solution:

      x(t) = e^(-1/tau(t-t0))*x(t0) + \int_t0^t e^{-1/tau(t-t')} 1/tau f(t) dt'     (2)


 - With time varying coefficient, characteristic equation:

     x_dot(t) = a(t) x(t) + b(t) u(t)                                               (3)

   Generic solution:

   
"""
import numpy as np

# --------------------------------------------------------------------------------}
# --- Unforced response
# --------------------------------------------------------------------------------{
def zero_response(time, tau, x0=1, t0=0):
    r"""
    Response to zero input, simple decay,  f(t>=t0)=0 in eq.(1)
    """
    t = np.asarray(time)-t0
    y = np.zeros(time.shape)+x0
    y += np.exp(-(time-time[0])/tau)*x0  
    return

# --------------------------------------------------------------------------------}
# --- Characteristics transient responses
# --------------------------------------------------------------------------------{
def impulse_response(time, tau, t0=0, A=1, x0=0):
    r"""
    Response to a delta dirac (impulse),  f(t)=delta(t-t0) in eq.(1)
    """
    t = np.asarray(time)-t0
    y = np.zeros(time.shape)
    y[t>=0] = 1/tau*np.exp(-t[t>=0]/tau)*A 
    y += np.exp(-(time-time[0])/tau)*x0  
    return y

def step_response(time, tau, t0=0, A=1, offset=0, x0=0):
    r"""
    Response to a unit step f(t)=A*Step(t-t0) in eq.(1)
    """
    t = np.asarray(time)-t0
    y = np.zeros(time.shape)+offset
    y[t>=0] = A*(1-np.exp(-t[t>=0]/tau)) + offset
    y += np.exp(-(time-time[0])/tau)*x0  
    return y

def ramp_response(time, tau, t0=0, A=1, offset=0, x0=0):
    r"""
    Response to a ramp,  f(t)=A*(t-t0) in eq.(1)
    """
    t = np.asarray(time)-t0
    y = np.zeros(time.shape)+offset
    y[t>=0] = A*(t[t>=0]-tau*(1-np.exp(-t[t>=0]/tau))) + offset
    y += np.exp(-(time-time[0])/tau)*x0  
    return y


def hat_response(time, tau, T=1, t0=0, A=1, x0=0, method='sum'):
    r""" 
    Response to a hat function:
      T : full time length of hat 
      A : Amplitude of hat
    """
    t = np.asarray(time)-t0
    # Decay of initial condition
    y = np.exp(-(time-time[0])/tau)*x0  
    # Solution should be summation of ramp responses
    y1= ramp_response(time, tau, t0=t0-T/2, A=A/T*2)
    y2= ramp_response(time, tau, t0=t0    , A=-A/T*4)
    y3= ramp_response(time, tau, t0=t0+T/2, A=A/T*2)
    y+= y1+y2+y3
    return y

def sine_response(time, tau, t0=0, A=1, omega=1, phi=0, x0=0):
    r"""
    Response to a sine,  f(t)=A*sin(omega*(t-t0)+phi) in eq.(1) only for t>=t0 
    input is assumed to be 0 before t0
    """
    # Decay of initial condition
    y = np.exp(-(time-time[0])/tau)*x0  

    t = np.asarray(time)-t0
    b =t>=0
    t=t[b]
    # --- Formulation 1 (obtained by integrating sin replaced by exponentials)
    # Decay of steady state equation value at "t=t0"
    #y[b] += 1/(1+tau**2 *omega**2) * np.exp(-t/tau) * ( - np.sin(phi)+omega*tau*np.cos(phi))
    # Steady state response
    #y[b] += 1/(1+tau**2 *omega**2) * ( np.sin(omega*t+phi)  -omega*tau *np.cos(omega*t + phi) )

    # --- Formulation 2 (simplifying the sum of sin and cos in Formulation 1
    ot      = omega*tau
    phi_out = np.arctan2((-ot *np.sin(phi) - np.cos(phi)),(-ot*np.cos(phi) + np.sin(phi))  )
    A_out   = 1/np.sqrt(1+ot**2)
    # Decay of steady state equation value at "t=t0"
    y[b] += - A_out *np.exp(-t/tau) * np.cos(phi_out)
    # Steady state response
    y[b] += A_out * np.cos(omega*t + phi_out)
    return y


def integrate(t_eval, tau, u, x0=0, method='LSODA', **options):

    if method in ['LSODA','RK45','RK23','DOP853','Radau','BDF']:
        # NOTE:
        #  - 'DOP853', 'Radau' work good for varying input
        #  - 'LSODA' good for adaptative time step
        from scipy.interpolate import interp1d
        from scipy.integrate import  solve_ivp #odeint
        # Create an interpolant for u
        fU = interp1d(t_eval, u)

        odefun = lambda t, x : 1/tau*(fU(t) - x)

        res = solve_ivp(fun=odefun, t_span=[t_eval[0], t_eval[-1]], y0=[x0], t_eval=t_eval, method=method, vectorized=False, **options)   
        return res.y.ravel()

    elif method.lower() in ['duhamel','convolution']:

        return duhamel(t_eval, tau, u, x0=x0)

    else:
        raise NotImplementedError('Method {}'.format(method))


def duhamel(time, tau, u, x0=0):
    r""" 
    Compute time response of a single DOF system using Duhamel's integral

        x(t)    = \int_0^t Phi(t-t') *B* u(t')  dt' + Phi(t-t0)x0
        x(t)    = \int_0^t H(t-t')    * u(t')   dt' + Phi(t-t0)x0
    with
     - F: force
     - B: 1/tau
     - Phi: state transition, Phi(t) = e^(-t/tau)
     - H  : impulse response, H(t)   = e^(-t/tau)/tau

    INPUTS:
      - time : 1d-array of time stamps
      - tau  : time constant
      - u    : 1d-array of inputs at each time step
    OUTPUTS:
      - x 1d-array, dof value 

    """
    H = impulse_response(time, tau, t0=time[0], A=1)
    #dt = t_eval[1]-t_eval[0]
    #x  = np.convolve(u.ravel(), H.ravel()  )[:len(t_eval)]*dt
    from welib.tools.signal_analysis import convolution_integral
    x  = convolution_integral(time, u.ravel(), H )
    x += tau*H*x0
    return x
# --------------------------------------------------------------------------------}
# ---  
# --------------------------------------------------------------------------------{
def frequency_response(Omega, tau, num=None, method='analytical' ):
    """ 
    Frequency respone of a first order system defined by :
       H(s)  =  b2 s + b1  / (1 + a),   num=[b2, b1]
       tau = 1/a
    """
    # --- Sanity checks
    if num is None:
        num = [0, 1]
    if len(num)==1:
        num = [0] + list(num)
    if len(num)!=2:
        raise Exception('Function only works for a numerator of dimension 2 max')

    # --- Usefule parameters
    params ={}
    b_1 = num[1]
    b_2 = num[0]
    a   = 1/tau
    params['omega_c'] = 1/tau
    params['fc'] = params['omega_c']/(2*np.pi)
    if b_2 ==0:
        params['hmag_0']   = np.abs(b_1/a)
        params['hmag_inf'] = np.nan
        params['phi_0']   = 0
        params['phi_inf'] = -np.pi/2 *180/np.pi
    else:
        params['hmag_0']   = np.abs(b_1/a)
        params['hmag_inf'] = np.abs(b_2)
        params['phi_0']   = 0
        params['phi_inf'] = 0

    # --- Frequency response
    if method=='analytical':
        # Method 1 - analytical expression
        Hmag = np.sqrt(b_1**2 + b_2**2*Omega**2)/np.sqrt(a**2+Omega**2)
        phi  = np.arctan( ( (a*b_2-b_1)*Omega) /(a*b_1+b_2*Omega**2))

    elif method=='control':
        # Method 2 - using Control toolbox
        import control as ct
        den = [1, a]
        sys = ct.tf(num, den)
        Hdat = sys.frequency_response(Omega)
        H = Hdat.fresp[0,0,:]
        Hmag, phi = np.abs(H), np.angle(H)
    else:
        raise Exception('Wrong method', method)
    return Hmag, phi, params


def bode_plot(Omega, tau, num=None, fig=None, Hscale='0', fscale='fc'):
    import matplotlib.pyplot as plt
    freq = Omega / (2*np.pi)
    omega_c = 1/tau
    ic = np.argmin(np.abs(Omega-omega_c))

    Hmag_th, phi_th, params = frequency_response(Omega, tau, num=num, method='analytical')
    #Hmag, phi, params = frequency_response(Omega, tau, num=None, method='control')
    fc = params['fc']

    # --- Scales
    if Hscale is None:
        Hscale='none'
    if fscale is None or fscale=='none':
        fscale='hertz'
    Hscale=Hscale.lower()
    H0 = {'none':1, '0':params['hmag_0'], 'inf':params['hmag_inf']}[Hscale]
    fscale=fscale.lower()
    fc = {'hertz':1, 'fc':fc}[fscale]


    fig,axes = plt.subplots(2, 1, sharex=True, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    # Amplitudes

    axes[0].plot(freq     /fc, Hmag_th/H0 , ls='-', color='k')
#     axes[0].plot(freq[:ic]/fc, freq[:ic]*0+params['hmag_0']  /H0, ls=':', color='k')
#     axes[0].plot(freq[ic:]/fc, freq[ic:]*0+params['hmag_inf']/H0, ls=':', color='k')
    # Phases
    axes[1].plot(freq     /fc, phi_th*180/np.pi, ls='-', color='k')
#     axes[1].plot(freq[:ic]/fc, freq[:ic]*0+params['phi_0']  , ls=':', color='k')
#     axes[1].plot(freq[ic:]/fc, freq[ic:]*0+params['phi_inf'], ls=':', color='k')
#     axes[0].axvline(fc /fc, c='k', ls=':')
#     axes[1].axvline(fc /fc, c='k', ls=':')
    ylabel = {'none':r' $H$', '0':'ratio $H/H_0$ [-]', 'inf':r'ratio $H/H_\infty$ [-]'}[Hscale]
    xlabel = {'hertz':r' $f$ [Hz]', 'fc':'ratio $f/f_c$ [-]', 'rad':r' $\omega$ [rad/s]'}[fscale]
    axes[0].set_ylabel(r'Amplitude'+ylabel)
    axes[1].set_xlabel(r'Frequency '+xlabel)    
    axes[1].set_ylabel('Phase [deg]')    
    axes[0].set_xscale('log')
    axes[1].set_xscale('log')
    axes[1].set_xlim([np.min(freq/fc), np.max(freq/fc)])
    axes[0].set_xlim([np.min(freq/fc), np.max(freq/fc)])

    xl = axes[1].get_xticklabels()
    for ix, x in enumerate(xl):
        if str(x) == r"Text(1.0, 0, '$\\mathdefault{10^{0}}$')":
            xl[ix]=plt.Text(1,0, r'$1$')
    axes[1].set_xticklabels(xl)
    return fig
