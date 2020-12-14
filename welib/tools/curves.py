""" 
Tools to manipulate 2D curves

In particular: 

 - tools to interpolate a curve at equidistant curvilinear spacing.
 - improved sream Quiver compare to native pyplot tool, for better handling of the arrows placements

"""

import numpy as np
import unittest

def curve_coord(line=None):
    """ return curvilinear coordinate """
    x=line[:,0]
    y=line[:,1]
    s     = np.zeros(x.shape)
    s[1:] = np.sqrt((x[1:]-x[0:-1])**2+ (y[1:]-y[0:-1])**2)
    s     = np.cumsum(s)                                  
    return s


def curve_extract(line,spacing,offset=None):
    """ Extract points at equidistant space along a curve"""
    x=line[:,0]
    y=line[:,1]
    if offset is None:
        offset=spacing/2
    # Computing curvilinear length
    s = curve_coord(line)
    offset=np.mod(offset,s[-1]) # making sure we always get one point
    # New (equidistant) curvilinear coordinate
    sExtract=np.arange(offset,s[-1],spacing)
    # Interpolating based on new curvilinear coordinate
    xx=np.interp(sExtract,s,x);
    yy=np.interp(sExtract,s,y);
    return np.array([xx,yy]).T


def curve_interp(x=None,y=None,n=3,line=None):
    """ Interpolate a curves to equidistant curvilinear space between points"""
    MatOut=False
    if line is None:
        line = np.array([x,y]).T
    else:
        x=line[:,0]
        y=line[:,1]
        MatOut=True
    if type(x) is not np.ndarray:
        x=np.array(x)
        y=np.array(y)
    if len(x)!=len(y):
        raise Exception('Inputs should have the same length')
    if len(x)<2 or len(y)<2:
        if MatOut:
            return np.array([x,y]).T
        else:
            return x,y
    # Computing curvilinear length
    s = curve_coord(line)
    # New (equidistant) curvilinear coordinate
    sNorm=np.linspace(0,s[-1],n);
    # Interpolating based on new curvilinear coordinate
    xx=np.interp(sNorm,s,x);
    yy=np.interp(sNorm,s,y);
    if MatOut:
        return np.array([xx,yy]).T
    else:
        return xx,yy


# --------------------------------------------------------------------------------}
# --- Streamlines and quiver 
# --------------------------------------------------------------------------------{
def lines_to_arrows(lines,n=5,spacing=None,normalize=True):
    """ Extract "streamlines" arrows from a set of lines 
    Either: `n` arrows per line
        or an arrow every `spacing` distance
    If `normalize` is true, the arrows have a unit length
    """
    if spacing is None:
        # if n is provided we estimate the spacing based on each curve lenght)
        if type(n) is int:
            n=[n]*len(lines)
        n=np.asarray(n)
        spacing = [ curve_coord(l)[-1]/nn for l,nn in zip(lines,n)]
    try:
        len(spacing)
    except:
        spacing=[spacing]*len(lines)

    lines_s=[curve_extract(l,spacing=sp,offset=sp/2)         for l,sp in zip(lines,spacing)]
    lines_e=[curve_extract(l,spacing=sp,offset=sp/2+0.01*sp) for l,sp in zip(lines,spacing)]
    arrow_x  = [l[i,0] for l in lines_s for i in range(len(l))]
    arrow_y  = [l[i,1] for l in lines_s for i in range(len(l))]
    arrow_dx = [le[i,0]-ls[i,0] for ls,le in zip(lines_s,lines_e) for i in range(len(ls))]
    arrow_dy = [le[i,1]-ls[i,1] for ls,le in zip(lines_s,lines_e) for i in range(len(ls))]

    if normalize:
        dn = [ np.sqrt(ddx**2 + ddy**2) for ddx,ddy in zip(arrow_dx,arrow_dy)]
        arrow_dx = [ddx/ddn for ddx,ddn in zip(arrow_dx,dn)] 
        arrow_dy = [ddy/ddn for ddy,ddn in zip(arrow_dy,dn)] 
    return  arrow_x,arrow_y,arrow_dx,arrow_dy 


def seg_to_lines(seg):
    """ Convert list of segments to list of continuous lines"""
    def extract_continuous(i):
        x=[]
        y=[]
        # Special case, we have only 1 segment remaining:
        if i==len(seg)-1:
            x.append(seg[i][0,0])
            y.append(seg[i][0,1])
            x.append(seg[i][1,0])
            y.append(seg[i][1,1])
            return i,x,y
        # Looping on continuous segment
        while i<len(seg)-1:
            # Adding our start point
            x.append(seg[i][0,0])
            y.append(seg[i][0,1])
            # Checking whether next segment continues our line
            #print('cur start',seg[i][0,:]  , ' end ',seg[i][1,:])
            Continuous= all(seg[i][1,:]==seg[i+1][0,:])
            #print('nxt start',seg[i+1][0,:], ' end ',seg[i+1][1,:],'Conti:',Continuous)
            if not Continuous:
                # We add our end point then
                x.append(seg[i][1,0])
                y.append(seg[i][1,1])
                break
            elif i==len(seg)-2:
                # we add the last segment
                x.append(seg[i+1][0,0])
                y.append(seg[i+1][0,1])
                x.append(seg[i+1][1,0])
                y.append(seg[i+1][1,1])
            i=i+1
        return i,x,y
    lines=[]
    i=0
    while i<len(seg):
        iEnd,x,y=extract_continuous(i)
        lines.append(np.array( [x,y] ).T)
        i=iEnd+1
    return lines

def streamQuiver(ax,sp,*args,**kwargs):
    """ Plot arrows from streamplot data  
    The number of arrows per streamline is controlled either by `spacing` or by `n`.
    See `lines_to_arrows`.
    """
    # TODO this was changed for python 2.7
    if 'spacing' in kwargs.keys():
        spacing=kwargs['spacing']
    else:
        spacing=None
    if 'n' in kwargs.keys():
        n=kwargs['n']
        del kwargs['n']
    else:
        n=5

    # --- Main body of streamQuiver
    # Extracting lines
    seg   = sp.lines.get_segments() # list of (2, 2) numpy arrays
    # lc = mc.LineCollection(seg,color='k',linewidth=0.7)
    lines = seg_to_lines(seg)       # list of (N,2) numpy arrays
    # Convert lines to arrows
    ar_x, ar_y, ar_dx, ar_dy = lines_to_arrows(lines,spacing=spacing,n=n,normalize=True)
    # Plot arrows
    qv=ax.quiver(ar_x, ar_y, ar_dx, ar_dy, *args, **kwargs)
    return qv


# --------------------------------------------------------------------------------}
# --- TEST 
# --------------------------------------------------------------------------------{
class TestCurves(unittest.TestCase):
    def test_seg_to_lines(self):
        # --- Useful variables
        Seg01=np.array([[0,0],[1,1]])
        Seg12=np.array([[1,1],[2,2]])
        Seg21=np.array([[2,2],[1,1]])
        Seg02=np.array([[0,0],[2,2]])

        # --- Empty segment > empty line
        lines= seg_to_lines([])
        self.assertEqual(seg_to_lines([]),[])
        # --- One segment >  one line
        lines= seg_to_lines([Seg01])
        np.testing.assert_equal(lines,[Seg01])
        # --- One continuous line
        lines= seg_to_lines([Seg01,Seg12,Seg21])
        np.testing.assert_equal(lines[0],np.array([ [0,1,2,1],[0,1,2,1]]).T)
        # --- One segment and one continuous lines
        lines= seg_to_lines([Seg02,Seg01,Seg12])
        np.testing.assert_equal(lines[0],Seg02)
        np.testing.assert_equal(lines[1],np.array([ [0,1,2],[0,1,2]]).T)
        # --- One continuous lines, one segment
        lines= seg_to_lines([Seg01,Seg12,Seg02])
        np.testing.assert_equal(lines[1],Seg02)
        np.testing.assert_equal(lines[0],np.array([ [0,1,2],[0,1,2]]).T)
        # --- Two continuous lines
        lines= seg_to_lines([Seg01,Seg12,Seg02,Seg21])
        np.testing.assert_equal(lines[0],np.array([ [0,1,2],[0,1,2]]).T)
        np.testing.assert_equal(lines[1],np.array([ [0,2,1],[0,2,1]]).T)

    def test_curv_interp(self):
        # --- Emtpy or size 1
#         np.testing.assert_equal(curve_interp([],[],5),([],[]))
#         np.testing.assert_equal(curve_interp([1],[0],5),([1],[0]))
#         np.testing.assert_equal(curve_interp([1],[0],5),([1],[0]))
        # --- Interp along x
        x=[0,1,2]
        y=[0,0,0]
        xx,yy=curve_interp(x,y,5)
        np.testing.assert_equal(xx,[0,0.5,1,1.5,2])
        np.testing.assert_equal(yy,xx*0)
        # --- Interp diag
        x=[0,1]
        y=[0,1]
        xx,yy=curve_interp(x,y,3)
        np.testing.assert_equal(xx,[0,0.5,1])
        np.testing.assert_equal(yy,[0,0.5,1])
        # --- Interp same size
        xx,yy=curve_interp(x,y,2)
        np.testing.assert_equal(xx,x)
        np.testing.assert_equal(yy,y)
        # --- Interp lower size
        xx,yy=curve_interp(x,y,1)
        np.testing.assert_equal(xx,x[0])
        np.testing.assert_equal(yy,y[0])


# --------------------------------------------------------------------------------}
# --- Examples 
# --------------------------------------------------------------------------------{
def example_streamquiver():
    import numpy as np
    import matplotlib.pyplot as plt

    x = np.linspace(-1,1,100)
    y = np.linspace(-1,1,100)
    X,Y=np.meshgrid(x,y)
    u = -np.sin(np.arctan2(Y,X))
    v =  np.cos(np.arctan2(Y,X))

    xseed=np.linspace(0.1,1,4)

    fig=plt.figure()
    ax=fig.add_subplot(111)
    sp = ax.streamplot(x,y,u,v,color='k',arrowstyle='-',start_points=np.array([xseed,xseed*0]).T,density=30)
    qv = streamQuiver(ax,sp,spacing=0.5, scale=50)
    plt.axis('equal')
    plt.show()


if __name__ == "__main__":
    TestCurves().test_curv_interp()
#     unittest.main()
