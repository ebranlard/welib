""" 
Tools to manipulate 2D curves

In particular: 

 - tools to interpolate a curve at equidistant curvilinear spacing.
 - improved sream Quiver compare to native pyplot tool, for better handling of the arrows placements

"""

import numpy as np
import unittest
_DEFAULT_REL_TOL=0.001

def get_line(x, y, close_it=False):
    if close_it:
        x, y = close_contour(x, y)

    line = np.zeros((len(x), 2))
    line[:,0] = x
    line[:,1] = y
    return line
    #return np.array([x,y]).T

def curve_coord(x=None, y=None, line=None):
    """ return curvilinear coordinate """
    if line is not None:
        x = line[:,0]
        y = line[:,1]
    x = np.asarray(x)
    y = np.asarray(y)
    s     = np.zeros(x.shape)
    s[1:] = np.sqrt((x[1:]-x[0:-1])**2+ (y[1:]-y[0:-1])**2)
    s     = np.cumsum(s)                                  
    return s


def curve_extract(line, spacing, offset=None):
    """ Extract points at equidistant space along a curve
    NOTE: curve_extract and curve_interp do basically the same, one uses "spacing", the other uses n
    """
    x=line[:,0]
    y=line[:,1]
    if offset is None:
        offset=spacing/2
    # Computing curvilinear length
    s = curve_coord(line=line)
    offset=np.mod(offset,s[-1]) # making sure we always get one point
    # New (equidistant) curvilinear coordinate
    sExtract=np.arange(offset,s[-1],spacing)
    # Interpolating based on new curvilinear coordinate
    xx=np.interp(sExtract,s,x);
    yy=np.interp(sExtract,s,y);
    return np.array([xx,yy]).T


def curve_interp(x=None, y=None, n=3, s=None, line=None):
    """ Interpolate a curves to equidistant curvilinear space between points

    INPUTS:
      either:
       -  x,y : 1d arrays 
        or
       -  line: (n_in x 2) array
      either
        - n : number of points to interpolate
         or
        - s : array of curvilinear coordinates where well interpolate
    """
    # --- Sanitization
    MatOut=False
    if line is None:
        line = get_line(x,y)
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
            return get_line(x,y)
        else:
            return x,y

    # Computing curvilinear length
    s_old = curve_coord(line=line)
    if s is None:
        # New (equidistant) curvilinear coordinate
        sNorm = np.linspace(0, s_old[-1], n);
    else:
        # New curvilinear coordinate
        if s[0]<s_old[0] or s[-1]>s_old[-1]:
            raise Exception('curve_interp: Curvilinear coordinate needs to be between 0 and {}, currently it is between {} and {}'.format(s_old[-1], s[0], s[-1]))
        sNorm = s
    # Interpolating based on new curvilinear coordinate
    xx = np.interp(sNorm, s_old, x)
    yy = np.interp(sNorm, s_old, y)
    if MatOut:
        return get_line(xx, yy) 
    else:
        return xx,yy



# --------------------------------------------------------------------------------
# --- Contour tools 
# --------------------------------------------------------------------------------
# NOTE: merge with airfoils/shape.py

def contour_length_scale(x, y):
    """ return representative length scale of contour """
    lx = np.max(x)-np.min(x)
    ly = np.max(y)-np.min(y)
    return  max(lx, ly)

def contour_isClosed(x, y, reltol=_DEFAULT_REL_TOL):
    """ Return true if contour is closed """
    l = contour_length_scale(x, y)
    return np.abs(x[0]-x[-1])<l*reltol or np.abs(y[0]-y[-1])<l*reltol

def contour_remove_duplicates(x, y, reltol=_DEFAULT_REL_TOL):
    l = contour_length_scale(x, y)
    unique_points = []
    duplicate_points = []
    for x,y in zip(x, y):
        if all(np.sqrt((x-p[0])**2 + (y-p[1])**2) > reltol*l for p in unique_points):
            unique_points.append((x,y))
        else:
            duplicate_points.append((x,y))
    x = np.array([p[0] for p in unique_points])
    y = np.array([p[1] for p in unique_points])
    return x, y, duplicate_points

def close_contour(x, y, reltol=_DEFAULT_REL_TOL, force=False):
    """ Close contour, unless it's already closed, alwasys do it if force is True"""
    x = np.asarray(x)
    y = np.asarray(y)
    isClosed = contour_isClosed(x, y, reltol=reltol)
    if isClosed or force:
        x = np.append(x, x[0])
        y = np.append(y, y[0])
    return x, y
    
def reloop_contour(x, y, i):
    """
    Reloop a contour array so that it starts at a specific index.
    NOTE: duplicates should preferably be removed
    INPUTS:
    - contour_array: Numpy array of shape (n, 2) representing the contour of point coordinates.
    - i: Index where the relooped contour should start.
    OUTPUTS:
    - Relooped contour array.
    """
    #relooped_contour = np.concatenate((contour_array[i:], contour_array[:i]))
    x2 = np.concatenate((x[i:], x[:i]))
    y2 = np.concatenate((y[i:], y[:i]))
    return x2, y2

def opposite_contour(x, y, reltol = _DEFAULT_REL_TOL):
    """
    Make a clockwise contour counterclockwise and vice versa
    INPUTS:
    - contour_array: Numpy array of shape (n, 2) representing the contour of point coordinates.
    OUTPUTS:
    - opposite contour
    """
    isClosed = contour_isClosed(x, y, reltol=reltol)
    if not isClosed:
        # we close the contour
        x, y = close_contour(x, y, force=True, reltol=reltol)
    xopp=x[-1::-1]
    yopp=y[-1::-1]
    # If it was not closed, we remove duplicates
    if not isClosed:
        xopp, yopp, dupli = contour_remove_duplicates(xopp, yopp, reltol=reltol)
    return xopp, yopp

def contour_orientation(x, y):
    """
    Determine if a contour is clockwise or counterclockwise.
    
    INPUTS:
    - x: 1D array containing the x-coordinates of the contour nodes.
    - y: 1D array containing the y-coordinates of the contour nodes.
    
    OUTPUTS:
    - 'clockwise' if the contour is clockwise, 'counterclockwise' if it's counterclockwise.
    """
    # Compute the signed area
    signed_area = np.sum(x[:-1] * y[1:] - x[1:] * y[:-1])
    
    # Determine orientation
    if signed_area > 0:
        return 'counterclockwise' # Positive about z
    elif signed_area < 0:
        return 'clockwise' # Negative about z
    else:
        return 'undetermined'




def closest_point(x0, y0, x, y):
    """ return closest point to curve and index"""
    i = np.argmin((x - x0)**2 + (y - y0)**2)
    return x[i], y[i], i

def find_closest(X, Y, point, xlim=None, ylim=None):
    """Return closest point(s), using norm2 distance 
    if xlim and ylim is provided, these are used to make the data non dimensional.
    """
    # NOTE: this will fail for datetime
    if xlim is not None:
        x_scale = (xlim[1]-xlim[0])**2
        y_scale = (ylim[1]-ylim[0])**2
    else:
        x_scale = 1
        y_scale = 1

    norm2 = ((X-point[0])**2)/x_scale + ((Y-point[1])**2)/y_scale
    ind = np.argmin(norm2, axis=0)
    return X[ind], Y[ind], ind


def point_in_contour(X, Y, contour, method='ray_casting'):
    """
    Checks if a point is inside a closed contour.
 
    INPUTS:
      - X: scalar or array of point coordinates
      - Y: scalar or array of point coordinates
      - contour: A numpy array shape (n x 2), of (x, y) coordinates representing the contour.
            [[x1, y1]
               ...
             [xn, yn]]
         or [(x1,y1), ... (xn,yn)]
 
    OUTPUTS:
      - logical or array of logical: True if the point is inside the contour, False otherwise.
    """
    def __p_in_c_ray(x, y, contour):
        # --- Check if a point is inside a polygon using Ray Casting algorithm.
        n = len(contour)
        inside = False
        p1x, p1y = contour[0]
        for i in range(n+1):
            p2x, p2y = contour[i % n]
            if y > min(p1y, p2y):
                if y <= max(p1y, p2y):
                    if x <= max(p1x, p2x):
                        if p1y != p2y:
                            xinters = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                        if p1x == p2x or x <= xinters:
                            inside = not inside
            p1x, p1y = p2x, p2y
        return inside

    def __p_in_c_cv2(x, y, contour):
       # --- CV2
       import cv2 # pip install opencv-python
       # NOTE: not working well...
       contour = contour.reshape((-1,1,2))
       dist = cv2.pointPolygonTest(contour.T, P, True)
       if dist > 0:
           return True
       else:
           return False
    # --- End subfunctions

    # --- Sanity
    contour = np.asarray(contour)
    assert(contour.shape[1]==2)

    if np.isscalar(X):
        # Check if the point is inside the bounding box of the contour.
        if X > np.max(contour[:,0]) or X < np.min(contour[:,0]) or Y > np.max(contour[:,1]) or Y < np.min(contour[:,1]):
            return False
        if method=='cv2':
            return __p_in_c_cv2(X, Y, contour)
        elif method=='ray_casting':
            return __p_in_c_ray(X, Y, contour)
        else:
            raise NotImplementedError()
    # --- For arrays
    #assert(X.shape==Y.shape)
    shape_in = X.shape
    Xf = np.array(X).flatten()
    Yf = np.array(Y).flatten()
    Bf = np.zeros(Xf.shape, dtype=bool)
    if len(Xf)!=len(Yf):
        raise Exception('point_in_contour: when array provided, X and Y must have the same length')
    # Quickly eliminate pints outside of bounding box (vectorial calculation)
    bbbox = (Xf <= np.max(contour[:,0])) & (Xf >= np.min(contour[:,0])) & (Yf <= np.max(contour[:,1])) & (Yf >= np.min(contour[:,1]))
    Bf[bbbox] = True
    for i, (x,y,b) in enumerate(zip(Xf, Yf, Bf)):
        if b: # If in Bounding box, go into more complicated calculation
            Bf[i] = __p_in_c_ray(x, y, contour)
    B =  Bf.reshape(shape_in)
    return B


# --------------------------------------------------------------------------------}
# --- Streamlines and quiver 
# --------------------------------------------------------------------------------{
def lines_to_arrows(lines, n=5, offset=None, spacing=None, normalize=True):
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
        spacing = [ curve_coord(line=l)[-1]/nn for l,nn in zip(lines,n)]
    try:
        len(spacing)
    except:
        spacing=[spacing]*len(lines)
    if offset is None:
        lines_s=[curve_extract(l, spacing=sp, offset=sp/2)         for l,sp in zip(lines,spacing)]
        lines_e=[curve_extract(l, spacing=sp, offset=sp/2+0.01*sp) for l,sp in zip(lines,spacing)]
    else:
        lines_s=[curve_extract(l, spacing=sp, offset=offset)         for l,sp in zip(lines,spacing)]
        lines_e=[curve_extract(l, spacing=sp, offset=offset+0.01*sp) for l,sp in zip(lines,spacing)]
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

def streamQuiver(ax, sp, n=5, spacing=None, offset=None, normalize=True, **kwargs):
    """ Plot arrows from streamplot data  
    The number of arrows per streamline is controlled either by `spacing` or by `n`.
    See `lines_to_arrows`.
    """
    # --- Main body of streamQuiver
    # Extracting lines
    seg   = sp.lines.get_segments() 
    if seg[0].shape==(2,2):
        #--- Legacy)
        # list of (2, 2) numpy arrays
        # lc = mc.LineCollection(seg,color='k',linewidth=0.7)
        lines = seg_to_lines(seg) # list of (N,2) numpy arrays
    else:
        lines = seg

    # Convert lines to arrows
    ar_x, ar_y, ar_dx, ar_dy = lines_to_arrows(lines, offset=offset, spacing=spacing, n=n, normalize=normalize)
    # Plot arrows
    qv=ax.quiver(ar_x, ar_y, ar_dx, ar_dy, **kwargs)
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


