""" 
Suite of tools to handle the welib repository

"""

import numpy as np
import os
import re
from welib.tools.figure import *
import matplotlib.pyplot as plt
from termcolor import colored, cprint

FIG_MD=[]
TIT_MD=[]

def export_figs_callback(filename):
    from welib.tools.repo import FIG_MD, TIT_MD
    script_dir = os.path.dirname(filename)
    setFigurePath('_figs/')
    figNames, filenames, titles = export2png(print_latex=False, verbose=False)
    print('filename:',filename)
    print('figNames:',figNames)
    print('titles:  ',titles)
    for fign, fn, t in zip(figNames,filenames,titles):
        TIT_MD+=['[{}](/{})'.format(t, filename.replace('\\','/'))]
        FIG_MD+=['![{}](/../figs/{})'.format(t, fn)]

def export_figs_rec(maindir):
    """ 
    Recursively loop in directory structure, look for example files
    Call them, with the "__export__" name, so that export_figs_callback is called.
    """
    try:
        os.mkdir('_figs')
    except:
        pass
    FIG_MD.clear()
    TIT_MD.clear()
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
                        cprint('[ OK ] {} figure(s)'.format(n2-n1), 'green')
                    else:
                        HAS_NOFIG.append(fullpath)
                        cprint('[INFO] No figure: {}'.format(fullpath), 'red')
    print('--------------------------------------------------------------')
    nCols= 5
    nRow= int(np.ceil(len(TIT_MD)/nCols))
    # --- print a summary
    print('Scripts with figures:')
    for f in HAS_FIG:
        cprint(f, 'green')
    print('Scripts without figures:')
    for f in HAS_NOFIG:
        cprint(f, 'red')
    print('')

    # --- Generate markdown for README.md
    k=0
    kk=0
    print(''.join(['| ']*nCols) +  ' |')
    print(''.join(['| :-------------------------: ']*nCols) +' |')
    for i in np.arange(nRow):
        kk=k
        print('| ',end='')
        for j in np.arange(nCols):
            if k<len(TIT_MD):
                tit=TIT_MD[k]
                print(tit, end='')
            if j<nCols-1:
                print(' | ', end='')
            else:
                print(' |')
            k=k+1
        k=kk
        print('| ',end='')
        for j in np.arange(nCols):
            if k<len(TIT_MD):
                fig=FIG_MD[k]
                print(fig, end='')
            if j<nCols-1:
                print(' | ', end='')
            else:
                print(' |')
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

