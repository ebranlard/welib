""" 
Plot C_P,opt (lambda) according to actuator disc momentum theory with wake rotation

Reproduces:
 - figure 9.7 from [1] or 3.5 from [2]
 - Figure 9.8 from [1] or 3.6 from [2]

See also welib.BEM.idealrotors

References:
   [1]  Branlard (2017) Wind turbine aerodynamics and vorticity-based methods, Springer
   [2]  Manwell, et al. (2023) Wind Energy Explained, 3rd Edition, Wiley

"""
from welib.wt_theory.idealrotors import *

def plot_CPopt(lambda_):
    CP, a    = ADMTO_CP(lambda_, method='analytical')

    fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,3.6)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.97, top=0.96, bottom=0.16, hspace=0.10, wspace=0.13)
    ax.plot(lambda_, 16/27+lambda_*0, 'k-'   , label='Betz limit (no wake rotation)')
    ax.plot(lambda_, CP             , 'k--'  , label='Including wake rotation')
    ax.set_xlabel(r'Tip speed ratio, $\lambda$ [-]')
    ax.set_ylabel(r'Power coefficient, $C_P$ [-]')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_ylim([0,0.6])
    ax.set_xlim([0,10.0])
    ax.set_yticks([0,0.1,0.2,0.3,0.4,0.5,0.6])

    bbox2 = dict(boxstyle='square,pad=-0.1', fc='w', ec='none')
    ax.annotate( r'$C_P=\frac{16}{27}$' , xy=(0.78, 16/27+0.005), xycoords='data', ha="center", va="center", xytext=(0.78, 0.52), textcoords='data', arrowprops=dict(arrowstyle="->"),bbox=bbox2, fontsize=14)

    ax.legend(frameon=False, loc='center right')
    lw=0.9
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(lw)
    ax.tick_params(length=6, width=lw)
    ax.set_title('WT Theory - CP optimal')
    return fig


def plot_a_ap_opt(lambda_):
    TSR=np.max(lambda_)
    a, ap    = ADMTO_inductions(lambda_, method='analytical')
    fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,3.6)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.97, top=0.96, bottom=0.16, hspace=0.10, wspace=0.13)
    ax.plot(lambda_/TSR, 1/3+lambda_*0, 'k-'   , label=r'$a$ - Betz limit (no wake rotation)')
    ax.plot(lambda_/TSR, a            , 'k--'  , label=r'$a$ - Including wake rotation')
    ax.plot(lambda_/TSR, ap           , 'k:'   , label=r"$a'$ - Including wake rotation")
    ax.set_xlabel(r'Radius $r/R$ [-]')
    ax.set_ylabel(r'Induction factors [-]')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_ylim([0,0.35])
    ax.set_xlim([0,1.0])

    ax.legend(frameon=False, loc='center right')
    lw=0.9
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(lw)
    ax.tick_params(length=6, width=lw)
    ax.set_title('WT Theory - CP optimal')

    return fig


if __name__ == '__main__':
    lambda_ = np.linspace(0.01,10,1000)
    fig = plot_CPopt(lambda_)
    fig = plot_a_ap_opt(lambda_)
    plt.show()

if __name__=="__test__":
    pass # See tests/test_idealrotors.py

if __name__=="__export__":
    lambda_ = np.linspace(0.01,10,1000)
    fig = plot_CPopt(lambda_)
    fig = plot_a_ap_opt(lambda_)
    from welib.tools.repo import export_figs_callback
    export_figs_callback(__file__)

