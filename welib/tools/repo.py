""" 
Suite of tools to handle the welib repository

"""

import numpy as np
import os
import re
from welib.tools.figure import *
import matplotlib.pyplot as plt

FIG_MD=[]
TIT_MD=[]

def export_figs_callback(filename):
    from welib.tools.repo import FIG_MD, TIT_MD
    script_dir = os.path.dirname(filename)
    setFigurePath('_figs/')
    figNames, filenames, titles = export2png(print_latex=False)
    print('filename:',filename)
    print('figNames:',figNames)
    print('titles:',titles)
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
    print(FIG_MD)
    reobj = re.compile('[a-zA-Z0-9][a-zA-Z0-9_]*.py')
    for root,dirnames,filenames in os.walk(maindir):
        if os.path.basename(root)=='examples':
            for f in filenames:
                if reobj.match(f):
                    print('--------------------------------------------------------------')
                    fullpath=os.path.join(root,f)
                    print('Running example script: {}'.format(fullpath))
                    plt.close('all')
                    execfile(fullpath, {'__name__': '__export__', 'print': lambda *_:None})
    print('--------------------------------------------------------------')
    nCols= 5
    nRow= np.int(np.ceil(len(TIT_MD)/nCols))

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

