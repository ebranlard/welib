""" 
Simple tools for transfer functions

Tries to respect some interfaces from matalb and python control toolbox
See python control toolbox, xferfcn for more
"""
import numpy as np


def sine_perturb_one_input(component, nInputs, A, omega, uDC=None, tmin=50, nPerPeriod=50, nPeriods=3, interpolant=False):
    """ 
    """
    # Figure out a time long enough 
    T    = 2*np.pi/omega
    dt   = T/nPerPeriod
    tmax = nPeriods*T + tmin
    t    = np.arange(0,tmax+dt,dt)
    perturb = A * np.sin(omega*t)

    # Inputs
    Mu = np.zeros((nInputs, len(t)))
    if uDC is not None:
        for iu in range(nInputs):
            Mu[iu,:] = uDC[iu]
    Mu[component,:] += perturb

    if interpolant:
        return t, interp1d(t, Mu), perturb
    else:
        return t, Mu, perturb

def numerical_frequency_response(omega, calcOutput, nInputs, nOutputs, A_inputs=1, uDC=None, tmin=20, nPeriods=4, nPerPeriod=51, deg=True):
    """ 
    Compute the frequency response (amplitude ratio and phase shift) for a given function
    using sinusoidal signals as inputs, and comparing the amplitude and phase of the output
    wrt the input.

    # TODO detect lack of convergence/instability..

    INPUTS:
     - omega: array of frequencies where the response is to be computed
     - calcOutput: function with interface Y = calcOutput(t, U)
                   where t is an array of time of length nt, starting at 0
                         U is nu x nt array of inputs
                         Y is ny x nt array of outputs
            If the function contains states, it's the responsability of the function to integrate these states
     - nInputs: number of inputs nu
     - nOuputs: number of outputs ny
     - A_inputs: amplitude to use for the inputs
            if scalar provided: the same is used for all inputs
            otherwise `A_inputs` should be a list of length nu.
            the perturbation of input iu, will be A_inputs[iu] * np.sin(2 pi omega[j])
     - uDC: references values of all the inputs, array of size nu. By default, the inputs are 0.
            if a reference value is provided, the sinusoidal pertubation is added to the DC component
     - tmin: minimum time before comparing input and output (to elliminate transients)
     - nPeriods: how many periods are used for comparison for each frequency
     - nPerPeriod: number of time steps per frequencies

    """
    from welib.tools.signal_analysis import input_output_amplitude_phase

    if not hasattr(A_inputs,'__len__'):
        A_inputs=[A_inputs]*nInputs


    nFreq   = len(omega)
    G   = np.zeros((nOutputs, nInputs, nFreq))
    phi = np.zeros((nOutputs, nInputs, nFreq))

    # Loop on frequencies
    for ifreq, om in enumerate(omega):
        # loop on inputs
        for iu, A_u in enumerate(A_inputs):
            # Perturb input
            t, MU, perturb = sine_perturb_one_input(iu, nInputs, A_u, om, uDC=uDC, tmin=tmin, nPeriods=nPeriods, nPerPeriod=nPerPeriod)
            # Compute outputs using user define function
            MY = calcOutput(t, MU)
            # Phase and amplitude for each outputs
            for iy in range(nOutputs): 
                G_, phi_ = input_output_amplitude_phase(t, perturb, MY[iy,:], A_u=A_u, omega_u=om, mask=t>=tmin, deg=deg)
                G[iy, iu, ifreq], phi[iy, iu, ifreq] =  G_, phi_
    return G, phi


class TransferFunction():
    pass


