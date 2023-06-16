import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local 


def normalize(x,y):
    c  = np.max(x)-np.min(x)
    xc = (x-np.min(x))/c # between 0 and 1
    yc = (y)/c
    return xc, yc

class AirfoilShape():
    def __init__(self, x=None, y=None, name=None, filename=None):
        # 
        if filename is not None:
            import welib.weio as weio
            df = weio.read(filename).toDataFrame()
            x = df.iloc[:,0]
            y = df.iloc[:,1]
            if name is None:
                name = os.path.splitext(os.path.basename(filename))[0]
        if name is None:
            name =''
        x = np.asarray(x)
        y = np.asarray(y)
        self.name = name
        self.filename = filename
        self.x = x
        self.y = y
        self._x = x.copy()
        self._y = y.copy()

    def write(self, filename, format='csv', delim=' '):
        s='# x/c_[-] y/c_[-]\n'
        for xx,yy in zip(self.x,self.y):
            s += '{:12.7f}{:s}{:12.7f}\n'.format(xx, delim,yy)
        with open ( filename, 'w' ) as f:
              f.write ( s )

    def leading_edge(self):
        xc, yc = normalize(self.x, self.y)
        xcLE, ycLE, i = find_closest(xc, yc, (0,0))
        return self.x[i], self.y[i], i

    def trailing_edge(self):
        # TODO consider interpolating and returning upper and lower neighbors 
        xc, yc = normalize(self.x, self.y)
        xcTE, ycTE, i = find_closest(xc, yc, (1,0))
        return self.x[i], self.y[i], i

    def split_surfaces(self):
        # TODO TODO Lot more work do be done here
        xLE, yLE, iLE = self.leading_edge ()
        xTE, yTE, iTE = self.trailing_edge()
        # Naive implementation
        if iLE==0:
            I1 = np.arange(iLE+1, iTE)
            I2 = np.arange(iTE+1, len(self.x))
        elif iTE==0:
            I1 = np.arange(iTE+1, iLE)
            I2 = np.arange(iLE+1, len(self.x))
        else:
            raise NotImplementedError()
        y1 =np.mean(self.y[I1])
        y2 =np.mean(self.y[I2])
        if y1>y2:
            Iu=I1
            Il=I2
        else:
            Il=I1
            Iu=I2

        xu = self.x[Iu]
        xl = self.x[Il]
        # Making sure x is increasing
        IIu = np.argsort(xu)
        IIl = np.argsort(xl)
        Iu = Iu[IIu]
        Il = Il[IIl]
        xu = self.x[Iu]
        xl = self.x[Il]
        yu = self.y[Iu]
        yl = self.y[Il]

        return xu, yu, xl, yl, Iu, Il

    def plot_surfaces(self):
        xu, yu, xl, yl, Iu, Il = self.split_surfaces()
        x0, y0, y0u, y0l = self.camberline()
        fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
        fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        ax.plot(self.x, self.y, 'k.', label='All')
        ax.plot(xu, yu, 'o'         , label='Upper', markerfacecolor='none')
        ax.plot(xl, yl, 'd'         , label='Lower', markerfacecolor='none')
        ax.plot(x0, y0, '--'  )
        #ax.plot(x0, y0u, ':'  )
        #ax.plot(x0, y0l, '-.'  )
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_title(name)
        plt.axis ( 'equal' )

    def camberline(self):
        xu, yu, xl, yl, Iu, Il = self.split_surfaces()
        x1 = min(np.min(xu), np.min(xl))
        x2 = max(np.max(xu), np.max(xl))
        x0 = np.linspace(x1, x2, int(len(self.x)/2)+1)
        y0u = np.interp(x0, xu, yu)
        y0l = np.interp(x0, xl, yl)

        y0 = (y0u+y0l)/2
        return x0, y0, y0u, y0l

    def plot(self):
        fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
        fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
        ax.plot(self.x, self.y, label='')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        plt.axis ( 'equal' )

    




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



def debug_slope():

    Up=np.array([[0.9947532, 0.0019938],
                 [0.9976658, 0.0015870],
                 [0.9994161, 0.0013419],
                 [1.0000000, 0.0012600]])

    Do=np.array([ [0.9947532, -.0019938],
                  [0.9976658, -.0015870],
                  [0.9994161, -.0013419],
                  [1.0000000, -.0012600]])

    ##X/c_[-]   Y/c_[-]
    #1.000000   0.001260
    #0.992704   0.002274
    #0.979641   0.004079
    #0.964244   0.006169
    #
    #0.964244  -0.006169
    #0.979641  -0.004079
    #0.992704  -0.002274
    #1.000000  -0.001260

    DxU = Up[:,0]-Up[-1,0]
    DyU = Up[:,1]-Up[-1,1]

    DxD = Do[:,0]-Do[-1,0]
    DyD = Do[:,1]-Do[-1,1]

    aU = np.arctan(DyU/DxU)*180/np.pi
    aD = np.arctan(DyD/DxD)*180/np.pi
    print('aU',aU[:-1])
    print('aD',aD[:-1])
    alpha=np.mean(-aU[:-1]+aD[:-1])
    print(alpha)


    import matplotlib.pyplot as plt

    fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
    ax.plot(Up[:,0], Up[:,1]    , label='')
    ax.plot(Do[:,0],Do[:,1]    , label='')


    x_=np.linspace(-1,0,10)
    ax.plot(x_+1,  np.tan(alpha/2*np.pi/180)*x_)
    ax.plot(x_+1, -np.tan(alpha/2*np.pi/180)*x_)


    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.legend()
    plt.show()

if __name__ == '__main__':
    plt.show()

