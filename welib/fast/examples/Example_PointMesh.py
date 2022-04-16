import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local 
from welib.fast.fast_mesh import *

pm = PointMesh(2)
pm.Position[0,:] = (0,0,0)
pm.Position[1,:] = (0,0,-10)

pm.Connectivity=np.array([[0,1]])

obj = pm.toConnectedObject()


# --- 
time =np.linspace(0,10,100)
x=10*np.sin(1.5*time)
theta_y=np.sin(1.5*time)*10*np.pi/180



ms = MeshStorage(pm, time)


for it, t in enumerate(time):

    #pm.rigidBodyMotion(u=(x[it],0,0))
    pm.rigidBodyMotion(u=(x[it],0,0), theta=(0,theta_y[it],0), RefPoint=(0,0,0))
    ms.store(pm, it)


obj = pm.toConnectedObject()
print(obj)
obj = ms.toConnectedObject(pm)
print(obj)
ms.toJSON3D(pm, '_MeshMotion.json')


if __name__ == '__main__':
    pass
