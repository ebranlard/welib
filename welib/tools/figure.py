"""
Tools to manipulate and export figures

"""
from numpy import mod
import matplotlib as mpl
import matplotlib.pyplot as plt
# from cycler import cycler
# from .colors import rgb2hex, fColrs
# Should be placed in some kind of InitClear
# mpl.rcParams['lines.linewidth'] = 1.55
# mpl.rcParams['font.size'] = 15
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
        _global_params.path=path
    else:
        _global_params.path=[path]

def setFigureTitle(btitle):
    _global_params.btitle=btitle

def title2filename(title):
    title=''.join([s[0].capitalize()+s[1:] for s in title.split()])
    return sub(r'[%|:;.\[ \]\\=^*_/]','',title)

def findtitle(fig):
    axTitle=None
    title=''
    # storing the title, figure name
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
    def export(fig,figformat,i=1,n=1,width=None,height=None,figNameLast='',script_name='',script_run_dir='',script_run_date=''):
        if i is None:
            i=1
        # params (for now, using global params)
        params=_global_params


        title,axTitle=findtitle(fig)
        titleLatexSafe = sub(r"[_%^]", "", title)
        print('>>>>> TITLE',title)
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
        print('Export path: ',params.path)
        for ifp in range(len(params.path)):
            filename='%s%s.%s'%(params.path[ifp],figName,figformat);
            fig.savefig(filename)
            print('Figure file: ',filename)
            # restoring the title
            if axTitle is not None:
                axTitle.set_title(title)
            else:
                fig.suptitle(title)

        # --- Generating latex code 
        if mod(n,2)==0:
            if mod(i,2)==0:
                FigureExporter.print2figures(figName,figNameLast,titleLatexSafe,script_name,script_run_dir,script_run_date)
        else:
            if mod(i,2)==0:
                FigureExporter.print2figures(figName,figNameLast,titleLatexSafe,script_name,script_run_dir,script_run_date)
            else:
                if i==n:
                    FigureExporter.print1figure(figName,titleLatexSafe,script_name,script_run_dir,script_run_date)

        figNameLast=figName;
        return figNameLast

# --- Export call wrapper 
def export(figformat,fig=None,i=None,width=None,height=None):
    import pylab
    import inspect
    import os
    import os.path
    import datetime
    frame=inspect.stack()[2]
    script_name=os.path.basename(frame[0].f_code.co_filename)
    script_run_dir=os.getcwd()
    script_run_date=datetime.datetime.now().strftime('%Y/%m/%d')
    __exporter = FigureExporter()  

    if fig is None:
        # We'll loop over all figures
        figures=[manager.canvas.figure for manager in pylab.matplotlib._pylab_helpers.Gcf.get_all_fig_managers()]
        figNameLast=''
        for i, figure in enumerate(figures):
            figNameLast=__exporter.export(fig=figure,figformat=figformat,i=(i+1),n=len(figures),width=width,height=height,figNameLast=figNameLast,script_name=script_name,script_run_dir=script_run_dir,script_run_date=script_run_date)

    else:
        __exporter.export(fig=fig,figformat=figformat,i=i,width=width,height=height,script_name=script_name,script_run_dir=script_run_dir,script_run_date=script_run_date)
        pass
    
    for ifp in range(len(_global_params.path)):
        print('Figure saved in: %s'%_global_params.path[ifp]);
    print(' ');

    pass

def export2pdf(fig=None,i=None,width=None,height=None):
    export('pdf',fig=fig,i=i,width=width,height=height)
def export2png(fig=None,i=None,width=None,height=None):
    export('png',fig=fig,i=i,width=width,height=height)
def export2eps(fig=None,i=None,width=None,height=None):
    export('png',fig=fig,i=i,width=width,height=height)


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
    Monitors=get_monitors()
    if len(Monitors)>1:
        for m in get_monitors():
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
    SA         = get_monitors()[0]
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
            
    fW=int(W/nY)
    fH=int(H/nX)
    #print('Window Dimensions: Width={} Height={}'.format(fW,fH))
    c=0
    for i in range(nX):
        for j in range(nY):
            fx=int(i*fH+Area.y)
            fy=int(j*fW+Area.x)
            s='{:d}x{:d}+{:d}+{:d}'.format(fW,fH,fy,fx)
            figs[c].canvas.manager.window.wm_geometry(s) # NOTE for TK, QT TODO
            c += 1
    if not bKeepFirstFig:
        plt.close(fig1)
    return figs

def fig_move(fig, AreaName=None,ScreenName='rightscreen',Area=None):
    """ Moves a figure to a given place on the screen """ 
    if Area is None:
        Area=_getArea(AreaName,ScreenName)
    s='{:d}x{:d}+{:d}+{:d}'.format(int(Area.width),int(Area.height),int(Area.x),int(Area.y))
    fig.canvas.manager.window.wm_geometry(s) # NOTE for TK, QT TODO
    return fig
    

# --------------------------------------------------------------------------------
# --- Example/tests
# --------------------------------------------------------------------------------
def test_fig_move()
    import matplotlib.pyplot as plt
    # --- Test fig_move
    fig=plt.figure()
    fig_move(fig,AreaName='Right')
    plt.show()
    # --- Test fig_grid
#     plt.ion()
#     figs=fig_grid(2,2)
#     # figs=[]
#     # figs+=fig_grid(AreaName='topright')
#     # figs+=fig_grid(AreaName='topleft')
#     # figs+=fig_grid(AreaName='bottomleft')
#     # figs+=fig_grid(AreaName='bottomright')
#     #figs=create_fig_grid(2,2,Name='rightscreen')
#     plt.show()
# 
#     for f in figs:
#         f.add_subplot(111)
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
    test_export()
