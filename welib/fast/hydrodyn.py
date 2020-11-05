import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local 
import weio

from welib.FEM.graph import *
from welib.tools.clean_exceptions import *
# import welib.FEM.fem_model as fem


if __name__ == '__main__':
    filename='../../_data/Monopile/MT100_HD.dat'
    filename='../../_data/Monopile/TetraSpar_HydroDyn_v2.dat'

    hd = weio.FASTInputFile(filename)
    hd.write('Out.dat')

    Graph = hd.toGraph(hd)

    Graph.divideElements(3)

    print(Graph)


    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib import collections  as mc
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = fig.add_subplot(1,2,1,projection='3d')
    lines=Graph.toLines(output='coord')
    for l in lines:
    #     ax.add_line(l)
        ax.plot(l[:,0],l[:,1],l[:,2])

    ax.autoscale()
    # ax.set_xlim([-40,40])
    # ax.set_ylim([-40,40])
    # ax.set_zlim([-40,40])
    # ax.margins(0.1)
    plt.show()

# 
# if __name__ == '__main__':
#     pass
