import numpy as np
import matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import itertools
from colorsys import rgb_to_hls as rgb_to_hls_scalar
from colorsys import hls_to_rgb as hls_to_rgb_scalar
import unittest

# --------------------------------------------------------------------------------}
# --- COLOR TOOLS 
# --------------------------------------------------------------------------------{
def rgb(r,b,g):
    if (r+b+g)<=3:
        return np.array([r,g,b])
    else:
        return np.array([r/255.,g/255.,b/255.])

def standardize(col, rgb_int=False):
    """ given a color return a rgb tuple between 0 and 1 """
    if rgb_int:
        return np.asarray(standardize(col), rgb_int=False)*255

    if isinstance(col, str):
        if len(col) in [6,7]:
            return hex2rgb(col)
        else:
            raise NotImplementedError('Strings of length {}'.format(len(col)))
    # 
    col = np.asarray(col)
    if len(col.shape)==1:
        # NOTE: this is approximative...
        if np.any(col>1):
            col = col/255
        return col
    else:
        #  TODO for loop and call standardize
        raise NotImplementedError('Arrays of colors')


def rgb2hex(C, g=None, b=None):
    if len(C)==3 :
        r=C[0]
        g=C[1]
        b=C[2]
    if r<1.1 and g<1.1 and b<1.1:
        r=r*255
        g=g*255
        b=b*255
        
    return '#%02X%02X%02X' % (r,g,b)

def hex2rgb(hx):
    """ Hex string to RGB tuple between 0 and 1 """
    hx = hx.lstrip('#')
    rgb = tuple(int(hx[i:i+2], 16)/255. for i in (0, 2, 4))
    return rgb

def adjust_color_lightness_scalar(r, g, b, factor):
    """
    r,g,b between 0 and 1
    factor between 0 and +infty, but lightness bounded between 0 and 1
    """
    h, l, s = rgb_to_hls_scalar(r, g, b)
    l = max(min(l * factor, 1.0), 0.0)
    r, g, b = hls_to_rgb_scalar(h, l, s)
    return r,g,b

def adjust_color_lightness(rgb, factor):
    """
    r,g,b between 0 and 1
    factor between 0 and +infty, but lightness bounded between 0 and 1
    """
    hls = rgb_to_hls(rgb)
    hls[...,1] = np.multiply(hls[...,1] , factor)
    hls[ hls[...,1]<0, 1] = 0
    hls[ hls[...,1]>1, 1] = 1
    rgb = hls_to_rgb(hls)
    return rgb


def lighten_color(rgb, factor=0.1):
    rgb = standardize(rgb)
    if factor ==0:
        return rgb
    return adjust_color_lightness(rgb, 1 + factor)

def darken_color(rgb, factor=0.1):
    rgb = standardize(rgb)
    return adjust_color_lightness(rgb, 1 - factor)

# --------------------------------------------------------------------------------}
# --- COLOR MAPS 
# --------------------------------------------------------------------------------{
def python_colors(i=None):
    if i is None:
        return plt.rcParams['axes.prop_cycle'].by_key()['color']
    else:
        Colrs=plt.rcParams['axes.prop_cycle'].by_key()['color']
        return Colrs[ np.mod(i,len(Colrs)) ]

# ---- ColorMap
def make_colormap(seq,values=None,name='CustomMap'):
    """Return a LinearSegmentedColormap
    seq: RGB-tuples. 
    values: corresponding values (location betwen 0 and 1)

    cmap=make_colormap([MW_Blue, [1.0,1.0,1.0],MW_Red])
    """
    hasAlpha=len(seq[0])==4
    if hasAlpha:
        nComp=4
    else:
        nComp=3

    n=len(seq)
    if values is None:
        values=np.linspace(0,1,n)

    doubled     = list(itertools.chain.from_iterable(itertools.repeat(s, 2) for s in seq))
    doubled[0]  = (None,)* nComp
    doubled[-1] = (None,)* nComp
    cdict = {'red': [], 'green': [], 'blue': []}
    if hasAlpha:
        cdict['alpha']=[]
    for i,v in enumerate(values):
        if hasAlpha:
            r1, g1, b1, a1 = doubled[2*i]
            r2, g2, b2, a2 = doubled[2*i + 1]
        else:
            r1, g1, b1 = doubled[2*i]
            r2, g2, b2 = doubled[2*i + 1]
        cdict['red'].append((v, r1, r2))
        cdict['green'].append((v, g1, g2))
        cdict['blue'].append((v, b1, b2))
        if hasAlpha:
            cdict['alpha'].append((v, a1, a2))
    print(cdict)
    return mcolors.LinearSegmentedColormap(name, cdict)

def cmap_colors(n, name='viridis'):
    try:
        cmap = matplotlib.colormaps[name]
        COLRS = [(cmap(v)[0],cmap(v)[1],cmap(v)[2]) for v in np.linspace(0,1,n)]
    except:
#         try:
#             import matplotlib.cm
#             cmap = matplotlib.cm.get_cmap('viridis')
#             COLRS = [(cmap(v)[0],cmap(v)[1],cmap(v)[2]) for v in np.linspace(0,1,n)]
#     except:
        print('[WARN] colors.py: cmap_colors failed.')
        COLRS = color_scales(n, color='blue')
    return COLRS

def color_scales(n, color='blue'):
    maps={
    'blue':mpl.cm.Blues,
    'purple':mpl.cm.Purples,
    'orange':mpl.cm.Oranges,
    'red':mpl.cm.Reds,
    'green':mpl.cm.Greens,
    'gray':mpl.cm.Greys,
    'grey':mpl.cm.Greys,
    }
    norm = mpl.colors.Normalize(vmin=0, vmax=n)
    cmap = mpl.cm.ScalarMappable(norm=norm, cmap=maps[color])
    cmap.set_array([])
    colrs = [cmap.to_rgba(i) for i in np.arange(n)]
    colrs.reverse()
    return colrs



# --- Diverging Brown Green
DV_BG=[rgb(140,81,10), rgb(191,129,4), rgb(223,194,1), rgb(246,232,1), rgb(245,245,2), rgb(199,234,2), rgb(128,205,1), rgb(53,151,14), rgb(1,102,94)]

# --- Diverging Red Blue
DV_RB=[ rgb(215,48,39), rgb(244,109,67), rgb(253,174,97), rgb(254,224,144), rgb(255,255,191), rgb(224,243,248), rgb(171,217,233), rgb(116,173,209), rgb(69,117,180) ]
 
# --- Diverging Purple Green
DV_PG=[rgb(118,42,131), rgb(153,112,171), rgb(194,165,207), rgb(231,212,232), rgb(247,247,247), rgb(217,240,211), rgb(166,219,160), rgb(90,174,97), rgb(27,120,55) ]

# Maureen Stone, for line plots
MW_Light_Blue    = rgb(114,147,203)
MW_Light_Orange  = rgb(225,151,76)
MW_Light_Green   = rgb(132,186,91)
MW_Light_Red     = rgb(211,94,96)
MW_Light_Gray    = rgb(128,133,133)
MW_Light_Purple  = rgb(144,103,167)
MW_Light_DarkRed = rgb(171,104,87)
MW_Light_Kaki    = rgb(204,194,16)

MW_Blue     =     rgb(57,106,177)
MW_Orange   =     rgb(218,124,48)
MW_Green    =     rgb(62,150,81)
MW_Red      =     rgb(204,37,41)
MW_Gray     =     rgb(83,81,84)
MW_Purple   =     rgb(107,76,154)
MW_DarkRed  =     rgb(146,36,40)
MW_Kaki     =     rgb(148,139,61)

MathematicaBlue       = rgb(63 ,63 ,153 );
MathematicaRed        = rgb(153,61 ,113 );
MathematicaGreen      = rgb(61 ,153,86  );
MathematicaYellow     = rgb(152,140,61  );
MathematicaLightBlue  = rgb(159,159,204 );
MathematicaLightRed   = rgb(204,158,184 );
MathematicaLightGreen = rgb(158,204,170 );
# 
ManuDarkBlue    = np.array([0   ,0  ,0.7 ])     ;
ManuDarkRed     = rgb(138 ,42 ,93  );
# ManuDarkOrange  = rgb(245 ,131,1   );
ManuDarkOrange  = rgb(198 ,106,1   );
ManuLightOrange = rgb(255.,212,96  );
# 
Red    = np.array([1  ,0  ,0]);
Blue   = np.array([0  ,0  ,1]);
Green  = np.array([0  ,0.6,0]);
Yellow = np.array([0.8,0.8,0]);

MatlabGreen   = np.array([0         ,0.5      ,1        ]);
MatlabCyan    = np.array([0.0e+0    ,750.0e-03,750.0e-03]);
MatlabMagenta = np.array([ 750.0e-03,0.0e+0   ,750.0e-03]);
MatlabYellow  = np.array([750.0e-03 ,750.0e-03,0.0e+0   ]);
MatlabGrey    = np.array([250.0e-03 ,250.0e-03,250.0e-03]);


# --- Pastel

Pastel_Beige     = rgb(250, 240, 215) # #FAF0D7
Pastel_Peach     = rgb(255, 217, 192) # #FFD9C0
Pastel_LigtRed   = rgb(246, 214, 214) # #F6D6D6
Pastel_Red       = rgb(223, 110, 110) # #df6e6e
Pastel_Yellow    = rgb(246, 247, 196) # #F6F7C4
Pastel_Bleue1    = rgb(140, 192, 222) # #8CC0DE
Pastel_Green1    = rgb(204, 238, 188) # #CCEEBC
Pastel_Green2    = rgb(158, 210, 190) # #9ED2BE
Pastel_Green3    = rgb(126, 170, 146) # #7EAA92
Pastel_Green4    = rgb( 27, 156, 133) # #1B9C85
Pastel_Green5    = rgb( 62, 142, 126) # #3E8E7E


# cRed=plt.cm.Reds(np.linspace(0.9,1,2))
# cGreen=plt.cm.Greens(np.linspace(0.9,1,2))
# cPur=plt.cm.Purples(np.linspace(0.9,1,2))
# cGray=plt.cm.Greys(np.linspace(0.9,1,2))

# --- Mathematica darkrainbow colormap:
darkrainbow=mcolors.LinearSegmentedColormap('CustomMap',
        {'red': [[0.0, None, 0.23529411764705882], [0.1111111111111111, 0.25098039215686274, 0.25098039215686274], [0.2222222222222222, 0.2627450980392157, 0.2627450980392157], [0.3333333333333333, 0.2901960784313726, 0.2901960784313726], [0.4444444444444444, 0.41568627450980394, 0.41568627450980394], [0.5555555555555556, 0.6235294117647059, 0.6235294117647059], [0.6666666666666666, 0.8117647058823529, 0.8117647058823529], [0.7777777777777777, 0.8745098039215686, 0.8745098039215686], [0.8888888888888888, 0.807843137254902, 0.807843137254902], [1.0, 0.7294117647058823, None]],
        'green': [[0.0, None, 0.33725490196078434], [0.1111111111111111, 0.3411764705882353, 0.3411764705882353], [0.2222222222222222, 0.4196078431372549, 0.4196078431372549], [0.3333333333333333, 0.4745098039215686, 0.4745098039215686], [0.4444444444444444, 0.5529411764705883, 0.5529411764705883], [0.5555555555555556, 0.6705882352941176, 0.6705882352941176], [0.6666666666666666, 0.7647058823529411, 0.7647058823529411], [0.7777777777777777, 0.7294117647058823, 0.7294117647058823], [0.8888888888888888, 0.5019607843137255, 0.5019607843137255], [1.0, 0.23921568627450981, None]] ,
        'blue': [[0.0, None, 0.5725490196078431], [0.1111111111111111, 0.5568627450980392, 0.5568627450980392], [0.2222222222222222, 0.3843137254901961, 0.3843137254901961], [0.3333333333333333, 0.27058823529411763, 0.27058823529411763], [0.4444444444444444, 0.23921568627450981, 0.23921568627450981], [0.5555555555555556, 0.2627450980392157, 0.2627450980392157], [0.6666666666666666, 0.30196078431372547, 0.30196078431372547], [0.7777777777777777, 0.3254901960784314, 0.3254901960784314], [0.8888888888888888, 0.2980392156862745, 0.2980392156862745], [1.0, 0.22745098039215686, None]]})
# NOTE: generated with:
MathematicaDarkRainbow=[(60 /255,86 /255,146/255), (64 /255,87 /255,142/255), (67 /255,107/255,98 /255), (74 /255,121/255,69 /255), (106/255,141/255,61 /255), (159/255,171/255,67 /255), (207/255,195/255,77 /255), (223/255,186/255,83 /255), (206/255,128/255,76 /255), (186/255,61 /255,58 /255)]
# darkrainbow2= make_colormap(MathematicaDarkRainbow)
# --- Another rainbow:
# Colrs=[(152/255,0,0),(152/255,69 /255,0 /255), (167/255,127/255,3  /255), (12 /255,137/255,0  /255), (0  /255,75 /255,131/255)]
# Colrs.reverse()

def fColrs_hex(*args):
    return rgb2hex(fColrs(*args))

def fGray(x):
    return [x,x,x]

def fColrs(i=-1, n=-1, bBW=True, cmap=None):
    # Possible calls
    # M=fColrs()  : returns a nx3 matrix of RBG colors
    # C=fColrs(i) : cycle through colors, modulo the number of color
    # G=fColrs(i,n) : return a grey color (out of n), where i=1 is black
    # % Thrid argument add a switch possibility between black and white or colors:
    # % G=fColrs(i,n,1) : return a grey color (out of n), where i=1 is black
    # % G=fColrs(i,n,0) : cycle through colors

    # Table of Color used
    if cmap is None:
        mcolrs=np.array([
                MathematicaBlue,
                MathematicaGreen,
                ManuDarkRed,
                ManuDarkOrange,
                MathematicaLightBlue,
                MathematicaLightGreen,
                MathematicaLightRed,
                ManuLightOrange,
                Blue,
                Green,
                Red,
                Yellow,
                MatlabCyan,
                MatlabMagenta ]);
    elif cmap=='darker':
        mcolrs=np.array(
                ['firebrick','darkgreen', MathematicaBlue, ManuDarkOrange], dtype=object )
    else:
        raise NotImplementedError()
        
    # 
    if i==-1:
        return mcolrs
    elif (i!=-1 and n==-1):
        return mcolrs[np.mod(i-1,len(mcolrs))];
    elif (i!=-1 and n!=-1):
        if bBW:
            if n==1:
                return [0,0,0]
            else:
                return [0.55,0.55,0.55]*(v-1)/(n-1); #grayscale
        else:
            return mcolrs[mod(i-1,len(mcolrs,1))]
    else:
        return mcolrs

# --------------------------------------------------------------------------------}
# --- Colorbar 
# --------------------------------------------------------------------------------{
def manual_colorbar(fig, cmap, norm=None, **kwargs):
    """ Adds a colorbar to a plot without linking it to a specific axis """
    if norm is None:
        norm=mcolors.Normalize(vmin=0, vmax=1)
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    return fig.colorbar(sm, **kwargs)
        
# 
def test_colrs():
    from matplotlib import pyplot as plt

    x=np.linspace(0,2*np.pi,100);
    plt.figure()
    plt.title('fig_python')
    plt.grid()
    for i in range(30):
        plt.plot(x,np.sin(x)+i,'-',color=fColrs(i))
    plt.xlabel('x coordinate [m]')
    plt.ylabel('Velocity  U_i [m/s]')
    plt.xlim([0,2*pi])

    plt.show()


def rgb_to_hls(rgb):
    """
    convert float rgb values (in the range [0, 1]), in a numpy array to hls
    values.

    Parameters
    ----------
    rgb : (..., 3) array-like or tuple
       All values must be in the range [0, 1]

    Returns
    -------
    hls : (..., 3) ndarray
       Colors converted to hls values in range [0, 1]
    """
    # --- Handling of arguments
    rgb = np.asarray(rgb)
    # check length of the last dimension, should be _some_ sort of rgb
    if rgb.shape[-1] != 3:
        raise ValueError("Last dimension of input array must be 3; "
                         "shape {shp} was found.".format(shp=rgb.shape))
    # if we got passed a 1D array, try to treat as a single color and reshape as needed
    in_ndim = rgb.ndim
    if in_ndim == 1:
        rgb = np.array(rgb, ndmin=2)

    hls = np.zeros_like(rgb) # We rely on zeros values for hue and saturation

    maxc  = rgb.max(-1) #  maxc = max(r, g, b)
    minc  = rgb.min(-1) #  minc = min(r, g, b)
    try:
        delta = np.ptp(rgb, -1) # max-min
    except:
        delta = rbg.ptp(-1) # max-min

    # --- Lightness
    hls[...,1] = (minc+maxc)/2.0
    # --- Saturation (HLS)
    idx = (maxc > 0) & (delta > 0)
    hls[idx,2] = delta[idx]/(1-np.abs(2*hls[idx,1]-1)) # s=(max-min)/(1-|2l-1|)
    # --- Hue
    ipos = delta > 0
    # red is max
    idx = (rgb[..., 0] == maxc) & ipos
    hls[idx, 0] = (rgb[idx, 1] - rgb[idx, 2]) / delta[idx]
    # green is max
    idx = (rgb[..., 1] == maxc) & ipos
    hls[idx, 0] = 2. + (rgb[idx, 2] - rgb[idx, 0]) / delta[idx]
    # blue is max
    idx = (rgb[..., 2] == maxc) & ipos
    hls[idx, 0] = 4. + (rgb[idx, 0] - rgb[idx, 1]) / delta[idx]
    hls[..., 0] = (hls[..., 0] / 6.0) % 1.0 

    if in_ndim == 1:
        hls.shape = (3,)
    return hls


def hls_to_rgb(hls):
    """
    convert hls values in a numpy array to rgb values
    all values assumed to be in range [0, 1]

    Parameters
    ----------
    hls : (..., 3) array-like or tuple
       All values assumed to be in range [0, 1]

    Returns
    -------
    rgb : (..., 3) ndarray
       Colors converted to RGB values in range [0, 1]
    """
    hls = np.asarray(hls)

    # check length of the last dimension, should be _some_ sort of rgb
    if hls.shape[-1] != 3:
        raise ValueError("Last dimension of input array must be 3; "
                         "shape {shp} was found.".format(shp=hls.shape))

    # if we got passed a 1D array, try to treat as
    # a single color and reshape as needed
    in_ndim = hls.ndim
    if in_ndim == 1:
        hls = np.array(hls, ndmin=2)

    # make sure we don't have an int image
    hls = hls.astype(np.promote_types(hls.dtype, np.float32))

    h = hls[..., 0]
    l = hls[..., 1]
    s = hls[..., 2]

    r = np.zeros_like(h)
    g = np.zeros_like(h)
    b = np.zeros_like(h)

    c  = (1-np.abs(2*l-1))*s
    hp = h*6.0
    x  = c * (1-np.abs( (hp % 2) -1 ))
    i  = hp.astype(int)
    m  = l-c/2

    idx = i % 6 == 0
    r[idx] = c[idx]
    g[idx] = x[idx]

    idx = i == 1
    r[idx] = x[idx]
    g[idx] = c[idx]

    idx = i == 2
    g[idx] = c[idx]
    b[idx] = x[idx]

    idx = i == 3
    g[idx] = x[idx]
    b[idx] = c[idx]

    idx = i == 4
    r[idx] = x[idx]
    b[idx] = c[idx]

    idx = i == 5
    r[idx] = c[idx]
    b[idx] = x[idx]

    r=r+m
    g=g+m
    b=b+m

    rgb = np.stack([r, g, b], axis=-1)

    if in_ndim == 1:
        rgb.shape = (3,)

    return rgb


def PRINT_COLORS(msg):
    # See welib.tools.strings
    HEADER = '\033[95m'
    RED = '\033[91m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    ORAN = '\033[93m'
    GREEN = '\033[92m'
    ENDC = '\033[0m'



class TestColors(unittest.TestCase):
    def test_rgb_hsv(self):
        # --- Test Data 
        RGB=np.zeros((11,3))
        RGB[0,:]=np.array([0.0,0.0,0.0])
        RGB[1,:]=np.array([1.0,1.0,1.0])
        RGB[2,:]=np.array([0.3,0.5,0.4])
        RGB[3,:]=np.array([0.9,0.9,0.9])
        RGB[4,:]=np.array([0.1,0.9,0.9])
        RGB[5,:]=np.array([0.9,0.1,0.9])
        RGB[6,:]=np.array([0.9,0.9,0.1])
        RGB[7,:]=np.array([0.9,0.1,0.1])
        RGB[8,:]=np.array([0.1,0.9,0.1])
        RGB[9,:]=np.array([0.1,0.1,0.9])
        RGB[10,:]=np.array([0.1,0.1,0.1])
        # --- Converting back and forth, 
        HLS  = rgb_to_hls(RGB)
        RGB2 = hls_to_rgb(HLS)
        np.testing.assert_almost_equal(RGB,RGB2)

        # --- Comparing results with scalar version
        for i in np.arange(RGB.shape[0]):
            h,l,s=rgb_to_hls_scalar(RGB[i,0],RGB[i,1],RGB[i,2])
            np.testing.assert_array_equal(HLS[i,:],[h,l,s])
        # --- Calling with tuple
        hls_arr       = rgb_to_hls((0.3,0.5,0.4))
        hls_tuple_ref = rgb_to_hls_scalar(0.3,0.5,0.4)
        np.testing.assert_equal(hls_arr,np.asarray(hls_tuple_ref))

    def test_adjust_lightness(self):
        RGB=np.zeros((10,3))
        RGB[0,:]=np.array([0.0,0.0,0.0])
        RGB[1,:]=np.array([1.0,1.0,1.0])
        RGB[2,:]=np.array([0.3,0.5,0.4])
        RGB[3,:]=np.array([0.9,0.9,0.9])
        RGB[4,:]=np.array([0.1,0.9,0.9])
        RGB[5,:]=np.array([0.9,0.1,0.9])
        RGB[6,:]=np.array([0.9,0.9,0.1])
        RGB[7,:]=np.array([0.9,0.1,0.1])
        RGB[8,:]=np.array([0.1,0.9,0.1])
        RGB[9,:]=np.array([0.1,0.1,0.9])
        factor=1.5
        RGB_out=adjust_color_lightness(RGB,factor)
        for i in np.arange(RGB.shape[0]):
            r,g,b=adjust_color_lightness_scalar(RGB[i,0],RGB[i,1],RGB[i,2],factor)
            np.testing.assert_almost_equal(RGB_out[i,:],[r,g,b])
        #  --- Calling with tuple
        rgb_in = (0.1,0.3,0.5)
        rgb_arr       = adjust_color_lightness(rgb_in,factor)
        rgb_tuple_ref = adjust_color_lightness_scalar(rgb_in[0],rgb_in[1],rgb_in[2],factor)
        np.testing.assert_almost_equal(rgb_arr,np.asarray(rgb_tuple_ref))

    def test_standardize(self):
        # 
        rgb_in = (255,0,0)
        rgb_out = standardize(rgb_in)
        np.testing.assert_almost_equal(rgb_out, [1,0,0], 5)

        rgb_in = (1,0,0)
        rgb_out = standardize(rgb_in)
        np.testing.assert_almost_equal(rgb_out, [1,0,0], 5)

        rgb_in = '#FF0000'
        rgb_out = standardize(rgb_in)
        np.testing.assert_almost_equal(rgb_out, [1,0,0], 5)

        rgb_in = '#FF0000'
        rgb_out = standardize(rgb_in)
        np.testing.assert_almost_equal(rgb_out, [1,0,0], 5)



if __name__ == "__main__":
#     test_colrs()
    unittest.main()
