import numpy as np



def dqdt_lorenz(t, q, sigma, beta, rho):
    """ 
    Lorenz system
    """
    x, y, z = q
     
    dx = sigma * (y - x)
    dy = x * (rho - z) - y
    dz = x * y - beta * z
     
    return [dx, dy, dz]

 
def generatePNG(filename, dpi=400):
    from scipy.integrate import odeint
    import matplotlib.pyplot as plt
    # --- Time integration
    p = (10, 8./3, 28)  # Parameters of the system
    res = odeint(dqdt_lorenz, [1.0,1.0,1.0], np.arange(0.0, 40, 0.01), p, tfirst=True) # NOTE: uses LSODA method
    # --- Plot
    fig,ax = plt.subplots(1, 1, sharey=False, figsize=(2.6,5.8)) # (6.4,4.8)
    ax.plot( res[:, 0],   res[:, 1],'k', lw=0.3)
    ax.set_xlim([-20, 20])
    ax.set_ylim([-27.5, 27.5])
    ax.axis('off')
    fig.tight_layout()
    fig.savefig(filename, dpi=dpi)

