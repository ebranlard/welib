""" 
Read/Write Mann Box

Part of weio library: https://github.com/ebranlard/weio

Mann box: 
  - z is the fast index, then y, then x
  - y from from  ly/2 to -ly/2  !<<<<<<< IMPORTANT, we will flip it 
  - z from from -lz/2 to  lz/2 
  - ix==1..nx corresponds to it = nt..1

The field stored in this class has the following properties:
  - shape: nx x ny x nz
  - y: goes from -ly/2 to ly/2  ! <<<<< IMPORTANT subtlety, it has been flipped
  - z: goes from -lz/2 to lz/2  
  - ix==1..nx corresponds to it = nt..1 (it has not been flipped)
"""
import pandas as pd
import numpy as np
import os
import re
import struct

try:
    from .file import File, EmptyFileError 
except:
    EmptyFileError = type('EmptyFileError', (Exception,),{})
    File=dict

class MannBoxFile(File):

    @staticmethod
    def defaultExtensions():
        return ['.u','.v','.w','.bin']

    @staticmethod
    def formatName():
        return 'Hawc2 turbulence box'

    def __init__(self,filename=None,**kwargs):
        self.filename = None
        if filename:
            self.read(filename=filename,**kwargs)

    def read(self, filename=None, N=None, dy=1, dz=1, y0=None, z0=0, zMid=None):
        """ read MannBox
             field (nx x ny x nz)
             NOTE: y-coord in Mann Box goes from Ly/2 -> -Ly/2 but we flip this to -Ly/2 -> Ly/2
        """
        if filename:
            self.filename = filename
        if not self.filename:
            raise Exception('No filename provided')
        if not os.path.isfile(self.filename):
            raise OSError(2,'File not found:',self.filename)
        if os.stat(self.filename).st_size == 0:
            raise EmptyFileError('File is empty:',self.filename)

        if N is None:
            # try to infer N's from filename with format 'stringN1xN2xN3'
            basename = os.path.splitext(os.path.basename(self.filename))[0]
            splits = basename.split('x')
            temp = re.findall(r'\d+', basename)
            res = list(map(int, temp))
            if len(res)>=3:
                N=res[-3:]
            else:
                raise BrokenFormatError('Reading a Mann box requires the knowledge of the dimensions. The dimensions can be inferred from the filename, for instance: `filebase_1024x32x32.u`. Try renaming your file such that the three last digits are the dimensions in x, y and z.')
        nx,ny,nz=N
        def _read_buffered():
            data=np.zeros((nx,ny,nz),dtype=np.float32)
            with open(self.filename, mode='rb') as f:            
                for ix in range(nx):
                    Buffer = np.frombuffer(f.read(4*ny*nz), dtype=np.float32) # 4-bytes
                    data[ix,:,:] =  np.flip(Buffer.reshape((ny,nz)),0)
            return data

        def _read_nonbuffered():
            data = np.fromfile(self.filename, np.dtype('<f'), -1)
            assert len(data) == nx*ny*nz, "Size of turbulence box (%d) does not match nx x ny x nz (%d)" % (len(data),nx*ny*nz)
            # Fortran order z the fastest, then y then x. We then flip that back to nx, ny, nz
            data = np.transpose(data.reshape(nz, ny, nx, order='F'), axes=(2,1,0))
            # The issue is the y-coordinate in Mann Boxes go from Ly/2 -> -Ly/2
            # So we flip the y-axis, so that the field is consistent with typical y values
            return np.flip(data, 1) # i.e. data=data[:,::-1,:]

#         self['field']= _read_nonbuffered()
        self['field']= _read_buffered()
        self['dy']=dy
        self['dz']=dz
        self['y0']=y0
        self['z0']=z0
        self['zMid']=zMid
#         print('1',self['field'][:,-1,0])
#         print('2',self['field'][0,-1::-1,0])
#         print('3',self['field'][0,-1,:])


    def write(self, filename=None):
        """ Write mann box """
        if filename:
            self.filename = filename
        if not self.filename:
            raise Exception('No filename provided')
        nx,ny,nz = self['field'].shape
        sfmt='<{:d}f'.format(ny*nz)
        with open(self.filename, mode='wb') as f:            
            for ix in np.arange(nx):
                data = np.flip(self['field'][ix,:,:],0).ravel() # We have to flip the y axis again
                f.write(struct.pack(sfmt, *data))

    
    def __repr__(self):
        s='<{} object> with keys:\n'.format(type(self).__name__)
        s+='| - filename: {}\n'.format(self.filename)
        s+='| - field:  shape {}x{}x{}\n'.format(self['field'].shape[0],self['field'].shape[1],self['field'].shape[2])
        s+='|   min: {}, max: {}, mean: {} \n'.format(np.min(self['field']), np.max(self['field']), np.mean(self['field']))
        s+='| - dy, dz:  {}, {}\n'.format(self['dy'], self['dz'])
        s+='| - y0, z0 zMid:  {}, {}, {}\n'.format(self['y0'], self['z0'], self['zMid'])
        s+='|useful getters: y, z, _iMid, fromTurbSim\n'
        z=self.z
        y=self.y
        s+='|   y: [{} ... {}],  dy: {}, n: {} \n'.format(y[0],y[-1],self['dy'],len(y))
        s+='|   z: [{} ... {}],  dz: {}, n: {} \n'.format(z[0],z[-1],self['dz'],len(z))
        return s

    def _iMid(self):
        _, ny, nz = self['field'].shape
        return int(ny/2), int(nz/2)

    @property
    def z(self):
        zmax = self['z0'] + (self['field'].shape[2]-1+0.1)*self['dz']
        z = np.arange(self['z0'], zmax, self['dz'])
        if self['zMid'] is not None:
            z+= self['zMid']-np.mean(z)
        return z

    @property
    def y(self):
        if self['y0'] is not None:
            ymax = self['y0'] + (self['field'].shape[1]-1+0.1)*self['dy']
            y = np.arange(self['y0'], ymax, self['dy'])
        else:
            ymax = (self['field'].shape[1]-1+0.1)*self['dy']
            y = np.arange(0, ymax, self['dy'])
            y -= np.mean(y)
        return y


    @property
    def vertProfile(self):
        iy, iz = self._iMid()
        m = np.mean(self['field'][:,iy,:], axis=0)
        s = np.std (self['field'][:,iy,:], axis=0)
        return self.z,m,s


    def toDataFrame(self):
        pass

    # Useful converters
    def fromTurbSim(self,u,icomp=0, removeConstant=None, removeAllMean=False):
        """ 
        Assumes: 
             u (3 x nt x ny x nz)
        Removes the mean of the turbsim file for the "u" component.
        """
        if icomp==0:
            if removeAllMean is True:
                self['field'] = u[icomp, :, : ,: ]-np.mean(u[icomp,:,:,:],axis=0)
            elif removeConstant is not None:
                self['field'] = u[icomp, :, : ,: ]-removeConstant
            else:
                self['field'] = u[icomp, :, : ,: ]
        else:
            self['field'] = u[icomp, :, : ,: ]

if __name__=='__main__':
    mb = MannBoxFile('mini-u_1024x32x32.bin')
#     mb = MannBoxFile('mann_bin/mini-u.bin', N=(2,4,8))
#     F1=mb['field'].ravel()
#     mb.write('mann_bin/mini-u-out.bin')
# 
#     mb2= MannBoxFile('mann_bin/mini-u-out.bin', N=(2,4,8))
#     F2=mb2['field'].ravel()
#     print(F1-F2)
