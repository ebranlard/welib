"""
Tools to manipulate and export figures

"""
import numpy as np
from numpy import mod
import matplotlib as mpl
import matplotlib.pyplot as plt
import re

# --- On load, set default rcParams
def defaultRC():
    # --- Ticks
    # ax.tick_params(direction='in', top=True, right=True, labelright=False, labeltop=False, which='both')
    mpl.rcParams['xtick.direction'] = 'in'
    mpl.rcParams['ytick.direction'] = 'in'
    mpl.rcParams['xtick.top']  = True
    mpl.rcParams['ytick.right'] = True

    # --- Axes limits
    # ax.autoscale(enable=True, axis='both', tight=True)
    mpl.rcParams['axes.autolimit_mode'] ='round_numbers'
    mpl.rcParams['axes.xmargin'] = 0
    mpl.rcParams['axes.ymargin'] = 0

    # --- Grid
    # ax.grid(True, linestyle=':')
    #mpl.rcParams['axes.grid']     = True
    #mpl.rcParams['grid.alpha']     = 1.0,
    mpl.rcParams['grid.color']     = '#b0b0b0'
    mpl.rcParams['grid.linestyle'] = ':'
    mpl.rcParams['grid.linewidth'] = 0.5 # default 0.8

    # --- Linewidth
    mpl.rcParams['lines.linewidth'] = 1.2  # default 1.5

    # --- Fontsize
    #mpl.rcParams['font.size'] = 15

    # --- Colors
    # from cycler import cycler
    # from .colors import rgb2hex, fColrs
    # Should be placed in some kind of InitClear
    # # mpl.rcParams['axes.color_cycle']=[rgb2hex((c*255).astype(int)) for c in fColrs()]
    # mpl.rcParams['axes.prop_cycle']=cycler(color=[rgb2hex((c*255).astype(int)) for c in fColrs()])
    # mpl.rcParams['font.family'] = 'helvetica'
    # line.set_dashes([8, 4, 2, 4, 2, 4]) 
    #lg = a.legend()
    # fr = lg.get_frame()
    # fr.set_lw(0.2)

# --------------------------------------------------------------------------------
# ---  
# --------------------------------------------------------------------------------
def dash(p):
    seq = [13, 5]
    p[0].set_dashes(seq)

def dotdash(p):
#     seq = [2, 4, 7, 4]
    seq = [3, 4, 13, 4]
    p[0].set_dashes(seq)


# --------------------------------------------------------------------------------}
# --- Color bar 
# --------------------------------------------------------------------------------{

def axesBB(axes):
    """ Bounding boxes of axes
    p: x1, y1, x2, y2
    """
    pmin=[1,1,1,1]
    pmax=[0,0,0,0]
    for ax in np.asarray(axes).flatten():
        chartBox = ax.get_position()
        #x, y, w, h = chartBox.x0, chartBox.y0, chartBox.width, chartBox.height
        #print('xywh',x,y,w,h)
        p    = np.around(ax.get_position().get_points().flatten(),3) # BBox: x1, y1,  x2 y2
        pmin = np.around([min(v1, v2) for v1,v2 in zip(pmin, p)] ,3)
        pmax = np.around([max(v1, v2) for v1,v2 in zip(pmax, p)] ,3)
    return pmin, pmax

def add_colorbar(fig, commonBar=True, commonLevels=True, 
        position='right',
        size=2, pad=5,
        title=None, title_kwargs=None,
        levels=None, vmin=None, vmax=None, 
        verbose=False,
        mappables=None, cax=None):
    """ 
    Add colorbar(s) to figure axes

    INPUTS:
     - commonBar:  one colorbar for all the mappable
     - commonLevels: use the same levels for all. 
          NOTE: it's important to pass the levels to the mappable BEFORE.
              e.g.:  contourf(x, y, M, levels=levels)
     - size: 
         when position in ['left']          colorbar width  in % of plot area width 
         when positin in ['bottom', 'top']  colorbar height in % of plot area height 
     - pad: 
     - levels: levels passed to contourf, typically: linspace(0,10,100) or int
    ADVANCED:
      - cax : axis used for the colorbar itself)
      - mappables: list of mappables, e.g. output of mappable=plt.contourf()

    """
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    orientation = {'top':'horizontal','bottom':'horizontal','right':'vertical','left':'vertical'}[position]
    # --- Misc init
    axes = fig.get_axes()
    axes_flat = np.asarray(axes).flatten()
    fig = axes_flat[0].figure
    # Get all contourf mappable objects
    if mappables is None:
        mappables=[]
        for ax in axes_flat:
            mappable = ax.collections
            if len(mappable)!=1:
                raise NotImplementedError('Number of mappable is not one for axis')
            mappables.append(mappable[0])
    assert(len(mappables)==len(axes_flat))

    # --- Check consistencies between all plots
    is_polar = any([isinstance(ax, plt.PolarAxes) for ax in axes_flat])
    if verbose:
        print('is_polar', is_polar, [isinstance(ax, plt.PolarAxes) for ax in axes_flat])
    cmaps = [mappable.get_cmap().name for mappable in mappables]
    if verbose:
        print('cmaps', cmaps)

    min_values = min(mappable.get_array().min() for mappable in mappables)
    max_values = max(mappable.get_array().max() for mappable in mappables)
    if verbose:
        print('Super Min:',min_values)
        print('Super Max:',max_values)
    if levels is not None and not isinstance(levels, int):
        vmax = np.max(levels)
        vmin = np.min(levels)
    else:
        if vmin is None:
            vmin=min_values
        if vmax is None:
            vmax=max_values

    # --- Color range
    if commonLevels:
#        if levels is None:
#            levels = np.linspace(vmin, vmax, 10)
#        elif isinstance(levels, int):
#            levels = np.linspace(vmin, vmax, levels)
        if verbose:
            print('>>> setting levels', [vmin, vmax])
        for mappable in mappables:
            mappable.set_clim([vmin, vmax])

    if commonBar:
        # --- Sanity checks
        if len(np.unique(cmaps)) != 1:
            raise Exception('Not all colormaps are the same', cmaps)
        # --- Setup cax : axes where the colorbar will be placed
        if cax is None:
            # Compute Bounding box of plot
            #p: x1, y1, x2, y2
            pmin, pmax = axesBB(axes)
            if verbose:
                print('pmin:',pmin)
                print('pmax:',pmax)
            # Use bounding box to figure out where to put the colorbar
            pad =pad/100 
            w = size/100
            plots_height = pmax[3]-pmin[1]
            plots_width  = pmax[2]-pmin[0]
            # BoundingBox: Left (x), Bottom (y), Width, Height
            if position == 'right' :
                BB = [pmax[2]+pad, pmin[1], w, plots_height]
            elif position == 'left' :
                BB = [pmin[0]-pad-w, pmin[1], w, plots_height]
            elif position == 'bottom' :
                BB = [pmin[0], max(pmin[1]-1.5*pad-w,0), plots_width, w]
            elif position == 'top' :
                BB = [pmin[0], pmax[3]+pad, plots_width, w]
            cax = fig.add_axes(BB)
            pcax = np.around(cax.get_position().get_points().flatten(),3)
            if any(pcax>1): 
                print('[WARN] Colorbar not fitting in plot area, use subplot_adjust to make room for it')
            #print('pmin',pmin)
            #print('pmax',pmax)
            if verbose:
                print('BB  ',BB)
                print('pcax', pcax)

        cbars = fig.colorbar(mappable, cax=cax, orientation=orientation)
        if isinstance(title, list):
            title = title[0]
        if title is not None:
            cbars.ax.set_title(title, **title_kwargs)
    else:
        # --- One CBar per plot
        cbars = []
        if not isinstance(title, list):
            titles = [title] * len(mappables)
        else:
            titles = title
        for ax,mappable,title in zip(axes_flat, mappables, titles):
            is_polar = isinstance(ax, plt.PolarAxes)
            if not is_polar:
                # Don't know how to create a cax when axes is polar
                divider = make_axes_locatable(ax)
                if orientation=='vertical':
                    cax = divider.append_axes("right", size= "{}%".format(size), pad="{}%".format(pad))
                else:
                    cax = divider.append_axes("bottom", size="{}%".format(size), pad="{}%".format(pad))
            cbar = fig.colorbar(mappable, cax=cax, orientation=orientation)
            if title is not None:
                cbar.ax.set_title(title, **title_kwargs)
            cbars.append(cbar)

    return cbars


# --------------------------------------------------------------------------------
# --- Export library to help export figures
# --------------------------------------------------------------------------------
""" 
 Export library to help export figures

Main entry point is export or export2pdf

"""
class FigureExportParams:
    def __init__(self):
        self.path=['./']
        self.width=6.4    # default matplotlib
        self.height=4.8   # default matplotlib
        self.title=None
        self.font=15
        self.btitle=False
        self.blatex=False
        self.bconstrained=False

_global_params = FigureExportParams()  

def default_fig_size(fig):
    # default is (6.4,4.8)
    fig.subplots_adjust(top=0.96,bottom=0.12,left=0.14,right=0.96,wspace=0.05,hspace=0.05)
    fig.set_size_inches(_global_params.width,_global_params.height,forward=True) 

def setFigureFont(font):
    font=int(font)
    _global_params.font=font
    mpl.rcParams.update({'font.size': font})
    mpl.rcParams['font.size'] = font

def setFigureWidth(width):
    _global_params.width=width

def setFigureHeigth(height):
    _global_params.height=height

def setFigurePath(path):
    if type(path)==list:
        _global_paramGeompath=path
    else:
        _global_params.path=[path]
    for i,p in enumerate(_global_params.path):
        if p[-1]=='/' or p[-1]=='\\':
            pass
        else:
            _global_params.path[i] = p+'/'

def setFigureTitle(btitle):
    _global_params.btitle=btitle

def title2filename(title):
    title=''.join([s[0].capitalize()+s[1:] for s in title.split()])
    return re.sub(r'[%|:;.\[ \]\\=^*_/]','',title)

def findtitle(fig):
    axTitle=None
    title=''
    # storing the title, figure name
    try:
        title=fig._title
    except:
        title=fig._suptitle
        if title is not None  and len(title.get_text())>0:
            title=title.get_text()
        else:
            for ax in fig.get_axes():
                title=ax.get_title()
                if title is not None and len(title)>0:
                    axTitle=ax
                    break
    return title,axTitle

class FigureExporter:
    @staticmethod
    def print1figure(figName,titleLatexSafe,script_name,script_run_dir,script_run_date):
            print('in \\autoref{fig:%s}'%figName);
            print('% ---------------------------------- FIGURE --------------------------------------')
            print('%% From script: %s, folder: %s, %s'%(script_name,script_run_dir,script_run_date));
            print('\\noindent\\begin{figure}[!htb]\\centering%');
            print('  \\includegraphics[width=0.49\\textwidth]{%s}'%figName);
            print('  \\caption{%s}\\label{fig:%s}%%'%(titleLatexSafe,figName));
            print('\\end{figure}');
            print('% --------------------------------------------------------------------------------')
            print(' ')

    @staticmethod
    def print2figures(figName,figNameLast,titleLatexSafe,script_name,script_run_dir,script_run_date):
            print('in \\autoref{fig:%s}'%figNameLast);
            print('% ---------------------------------- FIGURES -------------------------------------')
            print('%% From script: %s, folder: %s, %s'%(script_name,script_run_dir,script_run_date));
            print('\\noindent\\begin{figure}[!htb]\\centering%%');
#             print('  \\begin{subfigure}[b]{0.49\\textwidth}\\centering \\includegraphics[width=\\textwidth]{%s}\\caption{}\\label{fig:%s}\\end{subfigure}%%'%(figNameLast,figNameLast));
#             print('  \\begin{subfigure}[b]{0.49\\textwidth}\\centering \\includegraphics[width=\\textwidth]{%s}\\caption{}\\label{fig:%s}\\end{subfigure}%%'%(figName,figName));
            print('  \\hfill\\includegraphics[width=0.49\\textwidth]{%s}%%'%(figNameLast));
            print('  \\hfill\\includegraphics[width=0.49\\textwidth]{%s}\\hfill'%(figName));
            print('  \\caption{%s}\\label{fig:%s}%%'%(titleLatexSafe,figNameLast));
            print('\\end{figure}');
            print('% --------------------------------------------------------------------------------')
            print(' ')

    @staticmethod
    def export(fig,figformat,i=1,n=1,width=None,height=None,figNameLast='',script_name='',script_run_dir='',script_run_date='', print_latex=True, verbose=False, twoByTwo=True):
        if i is None:
            i=1
        # params (for now, using global params)
        params=_global_params
        title,axTitle=findtitle(fig)
        titleLatexSafe = re.sub(r"[_%^]", "", title)
        #print('>>>>> TITLE',title)
        # figure name from title or figure number
        if title=='' or (title is None):
            figName='%d'%i
        else:
            figName=title2filename(title);

        # remove figure title if needed
        if not params.btitle:
            if axTitle is not None:
                axTitle.set_title('')
            else:
                fig.suptitle('')
        
        if params.blatex :
            pass
#             for iax =1:length(axList)
#                 ax=axList(iax);
#                 format_ticks(ax,0.02,0.009,'FontSize',11); % font size useless, hgexport takes care of it

        if params.bconstrained:
            pass
        # Actually just using the loose figure is enough to keep what the user has inputted it seems
            #xlims=get(gca,'XLim');
            #ylims=get(gca,'YLim');
            #pause(1)
            #set(gca,'XLim',xlims);
            #set(gca,'YLim',ylims);

        # --- Exporting in figure pathc
        if verbose:
            print('Export path: ',params.path)
        for ifp in range(len(params.path)):
            filename='%s%s.%s'%(params.path[ifp],figName,figformat);
            fig.savefig(filename)
            if verbose:
                print('Figure file: ',filename)
            # restoring the title
            if axTitle is not None:
                axTitle.set_title(title)
            else:
                fig.suptitle(title)

        # --- Generating latex code 
        if print_latex:
            if twoByTwo:
                if mod(n,2)==0:
                    if mod(i,2)==0:
                        FigureExporter.print2figures(figName,figNameLast,titleLatexSafe,script_name,script_run_dir,script_run_date)
                else:
                    if mod(i,2)==0:
                        FigureExporter.print2figures(figName,figNameLast,titleLatexSafe,script_name,script_run_dir,script_run_date)
                    else:
                        if i==n:
                            FigureExporter.print1figure(figName,titleLatexSafe,script_name,script_run_dir,script_run_date)
            else:
                FigureExporter.print1figure(figName,titleLatexSafe,script_name,script_run_dir,script_run_date)

        figNameLast=figName;
        return figNameLast, filename, title

# --- Export call wrapper 
def export(figformat, fig=None, i=None, width=None, height=None, print_latex=True, verbose=True, twoByTwo=True, path=None):
    import pylab
    import inspect
    import os
    import os.path
    import datetime
    # --- Default argument:
    if path is not None:
        setFigurePath(path)

    frame=inspect.stack()[2]
    script_name=os.path.basename(frame[0].f_code.co_filename)
    script_run_dir=os.getcwd()
    script_run_date=datetime.datetime.now().strftime('%Y/%m/%d')
    __exporter = FigureExporter()  

    figNames  = []
    fileNames = []
    titles = []
    if fig is None:
        # We'll loop over all figures
        figures=[manager.canvas.figure for manager in pylab.matplotlib._pylab_helpers.Gcf.get_all_fig_managers()]
        figNameLast=''
        for i, figure in enumerate(figures):
            figNameLast, filename, title=__exporter.export(fig=figure,figformat=figformat,i=(i+1),n=len(figures),width=width,height=height,figNameLast=figNameLast,script_name=script_name,script_run_dir=script_run_dir,script_run_date=script_run_date,print_latex=print_latex, verbose=verbose, twoByTwo=twoByTwo)
            figNames.append(figNameLast)
            fileNames.append(filename)
            titles.append(title)

    else:
        figNameLast, filename, title = __exporter.export(fig=fig,figformat=figformat,i=i,width=width,height=height,script_name=script_name,script_run_dir=script_run_dir,script_run_date=script_run_date, print_latex=print_latex, verbose=verbose, twoByTwo=twoByTwo)
        figNames.append(figNameLast)
        fileNames.append(filename)
        titles.append(title)
    
    if verbose:
        for ifp in range(len(_global_params.path)):
            print('Figure saved in: %s'%_global_params.path[ifp]);
        print(' ');

    return figNames, fileNames, titles

def export2pdf(**kwargs):
    return export('pdf',**kwargs)
def export2png(**kwargs):
    return export('png',**kwargs)
def export2eps(**kwargs):
    return export('png',**kwargs)


# --------------------------------------------------------------------------------}
# --- Tools to move figures, set them on a grid on different screens
# --------------------------------------------------------------------------------{
""" 
Set of tools to move figures

"""
def getMonitors():
    from screeninfo import get_monitors
    return get_monitors()

def getRightScreenArea(LeftPanel=0,BottomPanel=0,TopPanel=0,RightPanel=0):
    Monitors=getMonitors()
    if len(Monitors)>1:
        for m in Monitors:
            if m.__repr__().find('+0+0')<=0:
                SA = m
    else:
        SA=Monitors[0]
    #try:
    #except:
    #    SA         = get_monitors()[0]
    SA.y      += LeftPanel
    SA.x      += TopPanel
    SA.width  = SA.width-LeftPanel-RightPanel
    SA.height = SA.height-BottomPanel-TopPanel
    return SA

def getLeftScreenArea(LeftPanel=105,BottomPanel=0,TopPanel=0,RightPanel=0):
    SA         = getMonitors()[0]
    SA.y      += LeftPanel
    SA.x      += TopPanel
    SA.width  = SA.width-LeftPanel-RightPanel
    SA.height = SA.height-BottomPanel-TopPanel
    return SA

def _getArea(AreaName=None,ScreenName=None):
    # Determining Area from a acombination of area name and screen name
    if ScreenName.lower()=='rightscreen':
        ScreenArea=getRightScreenArea()
    elif ScreenName.lower()=='leftscreen':
        ScreenArea=getLeftScreenArea()
    else :
        raise Exception('Unknown ScreenName `{}`'.format(ScreenName))

    #print('Screen Position  : X={} Y={}'.format(ScreenArea.x,ScreenArea.y))
    #print('Screen Dimensions: Width={} Height={}'.format(ScreenArea.width,ScreenArea.height))
    if AreaName is not None:
        if AreaName.lower()=='leftscreen':
            Area=getLeftScreenArea()
        elif AreaName.lower()=='rightscreen':
            Area=getRightScreenArea()
        elif AreaName.lower()=='top':
            Area=ScreenArea
            Area.height /=2
        elif AreaName.lower()=='bottom':
            Area=ScreenArea
            Area.height/=2
            Area.y +=Area.height
        elif AreaName.lower()=='left':
            Area=ScreenArea
            Area.width/=2
        elif AreaName.lower()=='right':
            Area=ScreenArea
            Area.width/=2
            Area.x +=Area.width
        elif AreaName.lower()=='center':
            Area=ScreenArea
            Area.width/=2
            Area.height /=2
            Area.x +=Area.width/2
            Area.y +=Area.height/2
        elif AreaName.lower()=='topright':
            Area=ScreenArea
            Area.width/=2
            Area.height /=2
            Area.x +=Area.width
        elif AreaName.lower()=='topleft':
            Area=ScreenArea
            Area.width/=2
            Area.height /=2
        elif AreaName.lower()=='bottomright':
            Area=ScreenArea
            Area.width/=2
            Area.height /=2
            Area.x +=Area.width
            Area.y +=Area.height
        elif AreaName.lower()=='bottomleft':
            Area=ScreenArea
            Area.width/=2
            Area.height /=2
            Area.y +=Area.height
        else:
            raise Exception('Screen area name not supported `{}`'.format(Name))
    else: 
        Area = ScreenArea
    return Area

def fig_grid(nX=1,nY=1,AreaName=None,ScreenName='rightscreen',Area=None):
    # Creating figures
    figs=[]
    fig1=plt.figure()
    #bKeepFirstFig=fig1.number!=1
    bKeepFirstFig=False
    for i in range(nX):
        for j in range(nY):
            if i==0 and j==0 and bKeepFirstFig:
                figs.append(fig1)
            else:
                figs.append(plt.figure())
            
    Area=_getArea(AreaName,ScreenName)
    fW=int(Area.width/nY)
    fH=int(Area.height/nX)
    #print('Window Dimensions: Width={} Height={}'.format(fW,fH))
    c=0
    for i in range(nX):
        for j in range(nY):
            fx=int(i*fH+Area.y)
            fy=int(j*fW+Area.x)
            _move_fig(figs[c], fW, fH, fx, fy)
            c += 1
    if not bKeepFirstFig:
        plt.close(fig1)
    return figs

def fig_move(fig, AreaName=None,ScreenName='rightscreen',Area=None):
    """ Moves a figure to a given place on the screen """ 
    if Area is None:
        Area=_getArea(AreaName,ScreenName)
    _move_fig(fig, Area.width, Area.height, Area.x, Area.y) 
    return fig
    
def _move_fig(fig, w, h, x, y):
    if isinstance(fig,int) : fig = plt.figure(fig)
    manager = fig.canvas.manager
    backend = mpl.get_backend()
    if backend == 'TkAgg':
        s='{:d}x{:d}+{:d}+{:d}'.format(int(w),int(h),int(x),int(y))
        manager.window.wm_geometry("+%d+%d" % (x, y))
    elif backend == 'WXAgg':
        manager.window.SetPosition((x, y))
    else:
        # This works for QT and GTK
        manager.window.setGeometry(y,x,w,h)
        #manager.window.move(x, y)


# --------------------------------------------------------------------------------
# --- Example/tests
# --------------------------------------------------------------------------------
def test_fig_move():
    import numpy as np
    import matplotlib.pyplot as plt
    # --- Test fig_move
#     fig=plt.figure()
#     fig_move(fig,AreaName='Right')
#     plt.show()
    # --- Test fig_grid
    plt.ion()
    figs=fig_grid(2,2)
#     # figs=[]
#     # figs+=fig_grid(AreaName='topright')
#     # figs+=fig_grid(AreaName='topleft')
#     # figs+=fig_grid(AreaName='bottomleft')
#     # figs+=fig_grid(AreaName='bottomright')
#     #figs=create_fig_grid(2,2,Name='rightscreen')
    for f in figs:
        f.add_subplot(111)
    figs[0].axes[0].plot(np.linspace(0,1,10))
    figs[2].axes[0].plot(np.linspace(0,1,10))

    plt.ioff()
    plt.show()
# 
# 
#     n=1000
#     for i in range(1):
#         # Plotting 
#         x=np.random.rand(1,10)
#         y=np.random.rand(1,10)
#         for f in figs:
#             f.axes[0].clear()
#             f.axes[0].scatter(x,y,s=10)
#         plt.pause(0.001)
#         # Computing something heavy
#         a = np.random.rand(n,n); ainv = inv(a)
#         # Plotting additional data
#         x=np.random.rand(1,10)
#         y=np.random.rand(1,10)
#         for f in figs:
#             f.axes[0].scatter(x,y,color='k',marker='+',s=80)
#         plt.pause(0.001)
#         # Computing something heavy
#         a = np.random.rand(n,n); ainv = inv(a)
#     plt.pause(3)    

def test_export():
    #from ebra.export import *
    from numpy import linspace,sin,pi
    from matplotlib import pyplot as plt
    setFigurePath('./')

    x=linspace(0,2*pi,100);
    plt.figure()
    plt.title('First Example Figure')
    plt.grid()
    plt.plot(x,sin(x),'-')
    plt.xlabel('x coordinate [m]')
    plt.ylabel('Velocity  U_i [m/s]')
    plt.xlim([0,2*pi])

    #plt.figure()
    #plt.title('Second Example Figure')
    #plt.grid()
    #plt.plot(x,sin(x),'-')
    #plt.xlabel('x coordinate [m]')
    #plt.ylabel('Velocity  U_i [m/s]')
    #plt.xlim([0,2*pi])

    export2pdf()

if __name__ == "__main__":
    test_fig_move()
    #test_export()
