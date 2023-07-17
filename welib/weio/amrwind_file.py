"""Read AMR-Wind NETCDF file

"""
import xarray as xr
import numpy as np

class AMRWind:
    """ 
    Read a AMR-Wind output file (.nc)
    """

    def __init__(self, filename, timestep, output_frequency, **kwargs):
        self.filename   = filename
        self.amrwind_dt = timestep
        self.output_dt  = timestep * output_frequency 
            

    def read(self, group_name):        
        """
        Parameters
        ----------
        
        group_name : str,
            group name inside netcdf file that you want to read, e.g. p_slice
        """   
        
        ds = xr.open_dataset(self.filename,group=group_name)    
    
        coordinates = {"x":(0,"axial"), "y":(1,"lateral"),"z":(2,"vertical")}
        c           = {}
        for coordinate,(i,desc) in coordinates.items():
            c[coordinate] = xr.IndexVariable( 
                                             dims=[coordinate],
                                             data=np.sort(np.unique(ds['coordinates'].isel(ndim=i))), 
                                             attrs={"description":"{0} coordinate".format(desc),"units":"m"}
                                            )
        c["t"]                    = xr.IndexVariable( 
                                                     dims=["t"],
                                                     data=ds.num_time_steps*self.output_dt,
                                                     attrs={"description":"time from start of simulation","units":"s"}
                                                    )    

        self.nt = len(c["t"])
        self.nx = len(c["x"])
        self.ny = len(c["y"])
        self.nz = len(c["z"])

        coordinates = {"x":(0,"axial","u"), "y":(1,"lateral","v"),"z":(2,"vertical","w")}    
        v           = {}    
        for coordinate,(i,desc,u) in coordinates.items():        
            v[u] = xr.DataArray(np.reshape(getattr(ds,"velocity{0}".format(coordinate)).values,(self.nt,self.nx,self.ny,self.nz)), 
                                 coords=c, 
                                 dims=["t","x","y","z"],
                                 name="{0} velocity".format(desc), 
                                 attrs={"description":"velocity along {0}".format(coordinate),"units":"m/s"})

        ds = xr.Dataset(data_vars=v, coords=v[u].coords)           
        ds.attrs = {"original file":self.filename}
        
        self.data = ds
