""" 
Suite of tools to handle the welib repository

"""

import numpy as np
import os
import re
import matplotlib.pyplot as plt
import builtins
from welib.tools.figure import *
from welib.tools.strings import OK, WARN, print_red, print_green

FIG_MD=[]
TIT_MD=[]
FIG_NM=[]

MyDir=os.path.dirname(__file__)

DESCRIPTIONS={}
DESCRIPTIONS['Airfoils'] = ('airfoils', 'airfoil and polar applications', r"""
Examples of applications:
- Manipulation of airfoil curves, find slopes, interpolate (see [airfoils](welib/airfoils/examples/))
- Run different dynamic stall models (e.g Oye or MHH/HGM model) (see [airfoils/DS](welib/airfoils/examples/))

Sample figures from examples:
""")
DESCRIPTIONS['BEM']          = ('BEM'      , 'Blade Element Momentum Theory', r"""
Examples of Blade Element Momentum (BEM) Theory applications.
- Run steady state BEM simulations (see [BEM/steady 1-2](welib/BEM/examples)
- Run unsteady BEM simulations (see [BEM/unsteady 1-2](welib/BEM/examples/)

Sample figures from examples:
""")

DESCRIPTIONS['Beam']          = ('beam'      , 'beam theory and numerical analyses', r"""
Sample figures from examples:
""")


# - Controls applications (packages `ctrl`, `kalman`):
#     - Run a kalman filter to estimate states of a system (see [kalman](welib/kalman/))

# Category: DynamicInflow   (1 entries)
# Category: FAST            (2 entries)
# Category: FEM             (2 entries)

DESCRIPTIONS['FAST']          = ('fast'      , 'OpenFAST tools', r"""
This package led to the `openfast\_toolbox`
Sample figures from examples:
""")

DESCRIPTIONS['FEM']          = ('FEM'      , 'finite element method', r"""
Examples of applications:
- Perform 2d/3d FEM analyses using beam/frame elements (see [FEM](welib/FEM/examples))
- Craig-Bampton / Guyan reduction of a structure (see [FEM](welib/FEM/examples))

Sample figures from examples:
""")

DESCRIPTIONS['Hydro'] = ('hydro', 'hydrodynamics applications', r"""
Examples of applications:
- Wave kinematics for linear waves (see [hydro/Ex1](welib/hydro/examples/Ex1_WaveKinematics.py))
- Generation of wave time series from a given spectrum (see [hydro/Ex3](welib/hydro/examples/Ex3_WaveTimeSeries.py))
- Computation of wave loads on a monopile (see [hydro/Ex4](welib/hydro/examples/Ex4_WaveLoads.py))

Sample figures from examples:
""")

DESCRIPTIONS['IECStandards'] = ('standards', '', r"""
Sample figures from examples:
""")


DESCRIPTIONS['PartDyn'] = ('partdyn', 'particle dynamics', r"""
Sample figures from examples:
""")

DESCRIPTIONS['Plot'] = ('plot', '', r"""
Sample figures from examples:
""")

DESCRIPTIONS['Stochastic'] = ('stoch', '', r"""
Manipulate stochastic variables.
Sample figures from examples:
""")

DESCRIPTIONS['System'] = ('system', 'system dynamics applications', r"""
Examples of applications:
- Linearize a non-linear system defined by a state and output equation (implicit or explicit) (see [system](welib/system/tests))
- Perform time integration of mechanical systems (see [system](welib/system/examples))

Sample figures from examples:
""")

DESCRIPTIONS['Tools'] = ('tools', 'misc tools', r"""
-  Spectral analyses, signal processing, time integration, vector analyses
Sample figures from examples:
""")

DESCRIPTIONS['Vortilib'] = ('vortilib', 'vortex dynamics, theory and methods', r"""
Sample figures from examples:
""")

DESCRIPTIONS['Wind'] = ('wind and welib/ws\_estimator', 'wind generation and estimation', r"""
- Generate stochastic [wind](welib/hydro/examples/Ex3_WaveTimeSeries.py) times series
- Estimate wind speed (see 'welib\ws\_estimator`))

Sample figures from examples:
""")

DESCRIPTIONS['WTTheory'] = ('wt_theory', 'Analytical wind turbine solutions', r""" 
Examples of applications:
- Theory of optimal circulation

Sample figures from examples:
""")

DESCRIPTIONS['weio'] = ('weio', 'Wind energy IO', r""" 
- Read and write common wind energy file formats (see [weio](welib/weio), a clone of [weio](http://github.com/ebranlard/weio/))
Sample figures from examples:
""")

DESCRIPTIONS['yams'] = ('yams', 'Yet Another Multibody Solver', r""" 
Set of tools to work with structural dyanmics.
Examples of applications:
- Setup the equation of motions for a multibody system with flexible members analytically or numerically (see [yams](welib/yams/tests))

""")






def export_figs_callback(filename):
    from welib.tools.repo import FIG_MD, TIT_MD, FIG_NM
    script_dir = os.path.dirname(filename)
    setFigurePath('_figs/')
    figNames, filenames, titles = export2png(print_latex=False, verbose=False)
    print('filename:',filename)
    print('figNames:',figNames)
    print('titles:  ',titles)
    for fign, fn, t in zip(figNames,filenames,titles):
        TIT_MD+=['[{}](/{})'.format(t, filename.replace('\\','/'))]
        FIG_MD+=['![{}](/../figs/{})'.format(t, fn)]
        FIG_NM+=[fign]


def myprint(*args, **kwargs):
    # lazy initialization on first call
    if not hasattr(myprint, "_file"):
        myprint._file = open("README_OUTPUT.md", "w", encoding="utf-8")

    # print to screen (normal print)
    builtins.print(*args, **kwargs)

    # also print to file (same formatting)
    file_kwargs = dict(kwargs)
    file_kwargs["file"] = myprint._file
    builtins.print(*args, **file_kwargs)

def export_figs_rec(maindir, dry=False):
    """ 
    Recursively loop in directory structure, look for example files
    Call them, with the "__export__" name, so that export_figs_callback is called.
    """
    from welib.tools.repo import FIG_MD, TIT_MD, FIG_NM
    try:
        os.mkdir('_figs')
    except:
        pass
    FIG_MD.clear()
    TIT_MD.clear()
    FIG_NM.clear()
    HAS_FIG=[]
    HAS_NOFIG=[]
    reobj = re.compile('[a-zA-Z0-9][a-zA-Z0-9_]*.py')
    for root,dirnames,filenames in os.walk(maindir):
        sp = re.split(r'/|\\', root)
        if any([s.startswith('_') for s in sp]):
            #print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> SKIPPING',root)
            continue

        if os.path.basename(root)=='examples':
            for f in filenames:
                if reobj.match(f):
                    print('--------------------------------------------------------------')
                    fullpath=os.path.join(root,f)
                    print('Running: {}'.format(os.path.relpath(f, MyDir)))
                    plt.close('all')
                    n1=len(TIT_MD)
                    execfile(fullpath, {'__name__': '__export__', 'print': lambda *_:None})
                    n2=len(TIT_MD)
                    if n2>n1:
                        HAS_FIG.append(fullpath)
                        OK('{} figure(s)'.format(n2-n1))
                    else:
                        HAS_NOFIG.append(fullpath)
                        WARN('No figure: {}'.format(fullpath))
                    #if len(TIT_MD)>7:
                    #    break
    print('--------------------------------------------------------------')
    # --- print a summary
    print('Scripts with figures:')
    for f in HAS_FIG:
        print_green(f)
    print('Scripts without figures:')
    for f in HAS_NOFIG:
        print_red(f)
    print('')

    idx = sorted(range(len(FIG_NM)), key=FIG_NM.__getitem__)
    FIG_NM = [FIG_NM[i] for i in idx]
    FIG_MD = [FIG_MD[i] for i in idx]
    TIT_MD = [TIT_MD[i] for i in idx]

    # Extract category (string before first '-')
    categories = [s.split('-', 1)[0] for s in FIG_NM]
    uniq_cats = sorted(set(categories))
    print(f"Categories: {len(uniq_cats)} {uniq_cats}")
    for cat in uniq_cats:
        ind = [i for i, c in enumerate(categories) if c == cat]
        fig_nm_cat = [FIG_NM[i] for i in ind]
        print(f"Category: {cat:15s} ({len(ind)} entries)")

    # --- Print MarkDown to Screen



















    for cat in uniq_cats:
        ind = [i for i, c in enumerate(categories) if c == cat]
        fig_nm_cat = [FIG_NM[i] for i in ind]
        fig_md_cat = [FIG_MD[i] for i in ind]
        tit_md_cat = [TIT_MD[i] for i in ind]
        # --- Generate markdown for README.md FOR ALL
        if cat in DESCRIPTIONS:
            package = DESCRIPTIONS[cat][0]
            sdescr  = DESCRIPTIONS[cat][1]
            descr   = DESCRIPTIONS[cat][2]
            if len(sdescr)>0:
                myprint('## welib/{}: {}'.format(package, sdescr))
            else:
                myprint('## welib/{}'.format(package))
            if len(descr)>0:
                myprint(descr)
            del DESCRIPTIONS[cat]
        else:
            print('## welib/{}'.format(cat))
        print_MD_Figs(tit_md_cat, fig_md_cat, nCols=5)
        myprint('')
        myprint('')

    for k,v in DESCRIPTIONS.items():
        package = v[0]
        sdescr  = v[1]
        descr   = v[2]
        if len(sdescr)>0:
            myprint('## welib/{}: {}'.format(package, sdescr))
        else:
            myprint('## welib/{}'.format(package))
        if len(descr)>0:
            myprint(descr)
        myprint('')
        myprint('')
        



    # --- Generate markdown for README.md FOR ALL
    #print_MD_Figs(TIT_MD, FIG_MD, nCols=5)


def print_MD_Figs(Titles, Figures_MD, nCols=5):
    if len(Titles)<nCols:
        nCols = len(Titles)
    nRow= int(np.ceil(len(Titles)/nCols))
    k=0
    kk=0
    myprint(''.join(['| ']*nCols) +  ' |')
    myprint(''.join(['| :-------------------------: ']*nCols) +' |')
    for i in np.arange(nRow):
        kk=k
        myprint('| ',end='')
        for j in np.arange(nCols):
            if k<len(Titles):
                tit=Titles[k]
                myprint(tit, end='')
            if j<nCols-1:
                myprint(' | ', end='')
            else:
                myprint(' |')
            k=k+1
        k=kk
        myprint('| ',end='')
        for j in np.arange(nCols):
            if k<len(Titles):
                fig=Figures_MD[k]
                myprint(fig, end='')
            if j<nCols-1:
                myprint(' | ', end='')
            else:
                myprint(' |')
            k=k+1


def execfile(filepath, globals=None, locals=None):
    """ Execute a given python file """
    if globals is None:
        globals = {"__name__": "__main__"}
    globals.update({
        "__file__": filepath,
    })
    with open(filepath, 'rb') as file:
        exec(compile(file.read(), filepath, 'exec'), globals, locals)

