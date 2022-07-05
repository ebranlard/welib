""" 
Examples to use the Graph class and generate a JSON file that can be used with viz3danim

References:
 viz3danim: https://github.com/ebranlard/viz3danim


"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local 
from welib.FEM.graph import *


def main():
    # --- Creating a simple graph (two nodes, one element), export to JSON
    g = GraphModel()
    g.addNode(Node(1, 0,0,40))
    g.addNode(Node(2, 0,0,50))
    g.addElement(Element(nodeIDs=[1,2], ID=2))
    g.toJSON('_GraphE1N2_NoModes_NoTimeSeries.json')

    # --- Add Modes
    # TODO
    #g.addModes()

    # --- Create Time Series for JSON file
    time  = np.linspace(3,120,200)
    t0    = time[0]
    T     = time[-1]-t0
    f     = 1/T
    om    = 2*np.pi*f
    x     = 0.8*(time-t0)
    nt    = len(time)
    nN    = len(g.Nodes)
    displ = np.zeros((nt,nN,3))
    rot   = np.zeros((nt,nN,3,3))

    # First time series: linear translation of node 2
    for it, t in enumerate(time):
        displ[it, 1, 2 ] = 0.8*(t-t0)
        rot[it, 0, :,: ] =np.eye(3) # NOTE: for now rotations are not used
        rot[it, 1, :,: ] =np.eye(3)
    g.addTimeSeries(time, name='transNode2', displ=displ, rot=rot)

    # Second time series: sinusoidal variation of node 2
    displ = np.zeros((nt,nN,3))
    rot   = np.zeros((nt,nN,3,3))
    for it, t in enumerate(time):
        displ[it, 0, : ] = 0
        displ[it, 1, 1 ] = 3.*np.sin(om*(t-t0))
        rot[it, 0, :,: ] =np.eye(3) # NOTE: for now rotations are not used
        rot[it, 1, :,: ] =np.eye(3)
    g.addTimeSeries(time, name='sinNode2', displ=displ, rot=rot)

    # Third time series: rotation of element
    displ = np.zeros((nt,nN,3))
    rot   = np.zeros((nt,nN,3,3))
    for it, t in enumerate(time):
        displ[it, 0, : ] = 0
        displ[it, 1, 1 ] = 10.*np.sin(om*(t-t0)) 
        displ[it, 1, 2 ] = 10.*(np.cos(om*(t-t0)) -1) 
        rot[it, 0, :,: ] =np.eye(3) # NOTE: for now rotations are not used
        rot[it, 1, :,: ] =np.eye(3)
    g.addTimeSeries(time, name='RotNode2', displ=displ, rot=rot)

    # Export to JSON
    g.toJSON('_GraphE1N2_NoModes_WithTimeSeries.json')

    return g


if __name__ == '__main__':
    g = main()
    print(g)

if __name__ == '__test__':
    g = main()

    import os
    try:
        os.remove('_GraphE1N2_NoModes_NoTimeSeries.json')
        os.remove('_GraphE1N2_NoModes_WithTimeSeries.json')
    except:
        pass
