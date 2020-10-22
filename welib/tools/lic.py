import numpy as np
from colors import adjust_color_lightness
from curves import streamQuiver
try:
    from welib.tools.external.lic_internal  import line_integral_convolution
except:
    print("""
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Error loading lic_internal, please compile it using cython.

Find the package `pybra/external/` and follow the compilation instructions.

A substitude dummy function is provided for the current script to run
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

""")
    def line_integral_convolution(v,t,k):
        return np.ones(v[:,:,0].shape)*len(k)/2




def licImage(x,y,u,v,nLICKernel=31,kernel=None,texture=None,minSpeed=None,maxSpeed=None,accentuation=1.0,offset=0.1,spread=1,axial=True,cmap=None):
    u=np.asarray(u)
    v=np.asarray(v)

    if kernel is None:
        kernel = np.sin(np.arange(nLICKernel)*np.pi/nLICKernel)

    if texture is None:
        texture = np.random.rand(u.shape[1],u.shape[0]).astype(np.float32)
    nLICKernel=len(kernel)

    kernel = kernel*accentuation
    kernel = kernel.astype(np.float32)

    image=lic(u,v,texture=texture,kernel=kernel)
    image=image-np.mean(image)+nLICKernel/2 # Making sure all images have the same mean
    image=image/nLICKernel # scaling between 0 and 1
    image=offset+spread*image

    if axial:
        Speed=u
    else:
        Speed=np.sqrt((u**2+v**2))
    if minSpeed is None:
        minSpeed = Speed.min()
        maxSpeed = Speed.max()
    Speed[Speed>maxSpeed]=maxSpeed
    Speed[Speed<minSpeed]=minSpeed

    COL=cmap((Speed-minSpeed)/(maxSpeed-minSpeed)) # cmap requires value between 0 and 1
    COL=COL[:,:,0:3]
    MyImage=adjust_color_lightness(COL, image)

    return MyImage

def licPlot(x,y,u,v,ax,nLICKernel=31,texture=None,kernel=None,minSpeed=None,maxSpeed=None,accentuation=1.0,offset=0,spread=1,axial=True,nStreamlines=0,cmap=None):

    MyImage=licImage(x,y,u,v,nLICKernel=nLICKernel,texture=texture,kernel=kernel,minSpeed=minSpeed,maxSpeed=maxSpeed,accentuation=accentuation,offset=offset,spread=spread,axial=axial,cmap=cmap)

    # Background
    im=ax.imshow(MyImage,extent=[min(x),max(x),max(y),min(y)])

    if nStreamlines>0:
        yseed=np.linspace(min(y)*0.9,max(y)*0.9,nStreamlines)
        start=np.array([yseed*0,yseed])
        sp=ax.streamplot(x,y,u,v,color='k',start_points=start.T,linewidth=0.7,density=30,arrowstyle='-')
        qv=streamQuiver(ax,sp,n=5,scale=40,angles='xy')

    return im


def lic(u,v,texture=None,kernel=31):
    nx,ny=u.shape

    if texture is None:
        texture = np.random.rand(nx,ny).astype(np.float32)

    # Kernel
    if type(kernel) is int:
        n=kernel
        kernel = np.sin(np.arange(n)*np.pi/n)
        #For movies:
        #   kernel = np.sin(np.arange(n)*np.pi/n)*(1+np.sin(2*np.pi*freq*(np.arange(n)/float(n)+t)))
        #   kernel = np.sin(np.arange(n)*np.pi/n)*(1+np.sin(2*np.pi*freq*(np.arange(n)/float(n)+t)))
    elif type(kernel) is np.ndarray:
        pass
    else:
        raise Exception('Unsupported kernel type. Kernel should be an integer or a numpy array')
    kernel = kernel.astype(np.float32)

    vectors = np.zeros((nx,ny,2),dtype=np.float32)
    vectors[:,:,0]=u
    vectors[:,:,1]=v

    # calling internal function
    image = line_integral_convolution(vectors, texture, kernel)

    return image



def lic_flow(vectors,len_pix=10):
    """ MANU: was side to side with lic_internal????? """
    vectors = np.asarray(vectors)
    m,n,two = vectors.shape
    if two!=2:
        raise ValueError
    result = np.zeros((2*len_pix+1,m,n,2),dtype=np.int32) # FIXME: int16?
    center = len_pix
    result[center,:,:,0] = np.arange(m)[:,np.newaxis]
    result[center,:,:,1] = np.arange(n)[np.newaxis,:]
    for i in range(m):
        for j in range(n):
            y = i
            x = j
            fx = 0.5
            fy = 0.5
            for k in range(len_pix):
                vx, vy = vectors[y,x]
                print(x, y, vx, vy)
                if vx>=0:
                    tx = (1-fx)/vx
                else:
                    tx = -fx/vx
                if vy>=0:
                    ty = (1-fy)/vy
                else:
                    ty = -fy/vy
                if tx<ty:
                    print("x step")
                    if vx>0:
                        x+=1
                        fy+=vy*tx
                        fx=0.
                    else:
                        x-=1
                        fy+=vy*tx
                        fx=1.
                else:
                    print("y step")
                    if vy>0:
                        y+=1
                        fx+=vx*ty
                        fy=0.
                    else:
                        y-=1
                        fx+=vx*ty
                        fy=1.
                if x<0: x=0
                if y<0: y=0
                if x>=n: x=n-1
                if y>=m: y=m-1
                result[center+k+1,i,j,:] = y, x
    return result
