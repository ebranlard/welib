import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local 
from welib.tools.clean_exceptions import *
from welib.plot.json3d import *
from welib.FEM.graph import * # for dummy*, ConnectedObject

def main(verbose=False, Large=False):

    # --- Large file
    if Large:
        # time1=np.arange(0, 300, 0.01)
        time1=np.arange(0, 1, 0.1)
        np.random.seed(0)
        Nodes, Conn = dummyNodesConn(10, name='random')
        Props={'shape':'cylinder','type':0, 'Diam':5}
        Dummy = ConnectedObject('dummy', Nodes, Conn, Props)
        # Mode_1 = dummyModeDisplacements(Nodes, 'y-rigid-translation', 10.0)
        TS_1   = dummyTimeSeriesDisplacements(time1, Nodes, name = 'x-rot', Ampl = 20)
        js = JSON3DFile()
        js.addObject(Dummy)
        # js.addObjectMode('dummy', Mode_1,  'mode1')
        js.addObjectTimeSeries('dummy', time1, displ=TS_1,  name='TS1')
        #print(js)
        #print('>>> Writting...')
        js.write('_Dummy.json') #, mat4=True)

    else:
        # --- Simple Small Objects
        time1=np.linspace(0, 1, 100)
        time2=np.linspace(0, 20, 100)
        NodesR1=np.array([[0,0,40],
                         [0,20,50],
                         [0,0,60]])
        ConncR1=np.array([[0,1], [0,2]])
        PropsR1={'shape':'cylinder','type':0, 'Diam':5}

        DisplR1_1 = dummyModeDisplacements(NodesR1, 'x-rigid-translation', 1.0)
        DisplR1_2 = dummyModeDisplacements(NodesR1, 'y-rigid-translation', 2.5)
        DisplR1_3 = dummyModeDisplacements(NodesR1, 'x-rigid-rotation', 10)
        TS_R1_1   = dummyTimeSeriesDisplacements(time1, NodesR1, name           = 'x-lin', Ampl = 10)
        TS_R1_2   = dummyTimeSeriesDisplacements(time1, NodesR1, name           = 'x-lin', Ampl = 10)
        TS_R1_3   = dummyTimeSeriesDisplacements(time1, NodesR1, name           = 'x-lin', Ampl = 10)
        TS_R1_4   = dummyTimeSeriesDisplacements(time1, NodesR1, name           = 'x-rot', Ampl = 10)

        Rod = ConnectedObject('rod1', NodesR1, ConncR1, PropsR1)
        print(Rod)



        NodesR2=np.array([[10,0,40],
                          [10,0,50]])
        ConncR2=np.array([[0,1]])
        PropsR2={'shape':'cylinder','type':2, 'Diam':10}
        DisplR2_1=np.array([[0,0,1],
                            [0,0,1]])
        DisplR2_2=np.array([[1,0,0],
                            [1,0,0]])

        TS_R2_1 =  dummyTimeSeriesDisplacements(time1, NodesR2, name='static', Ampl=1)
        TS_R2_2 =  dummyTimeSeriesDisplacements(time1, NodesR2, name='x-lin', Ampl=10)
        TS_R2_3 =  dummyTimeSeriesDisplacements(time1, NodesR2, name='y-lin', Ampl=10)
        TS_R2_4 =  dummyTimeSeriesDisplacements(time1, NodesR2, name='x-lin', Ampl=10)

        Rod2 = ConnectedObject('rod2', NodesR2, ConncR2, PropsR2)
        print(Rod2)

        js = JSON3DFile()
        js.addObject(Rod)
        js.addObject(Rod2)

        js.addObjectMode('rod1', DisplR1_1,  'mode1')
        js.addObjectMode('rod1', DisplR1_2,  'mode2')
        js.addObjectMode('rod1', DisplR1_3,  'mode3')
        js.addObjectMode('rod2', DisplR2_1, 'mode1')
        js.addObjectMode('rod2', DisplR2_2, 'mode2')
        js.addObjectMode('rod2', DisplR2_2, 'mode3')

        js.addObjectTimeSeries('rod1', time1, displ=TS_R1_1,  name='TS1')
        js.addObjectTimeSeries('rod1', time1, displ=TS_R1_2,  name='TS2')
        js.addObjectTimeSeries('rod1', time1, displ=TS_R1_3,  name='TS3')
        js.addObjectTimeSeries('rod1', time1, displ=TS_R1_4,  name='TS4')
        js.addObjectTimeSeries('rod2', time1, displ=TS_R2_1,  name='TS1')
        js.addObjectTimeSeries('rod2', time1, displ=TS_R2_2,  name='TS2')
        js.addObjectTimeSeries('rod2', time1, displ=TS_R2_3,  name='TS3')
        js.addObjectTimeSeries('rod2', time1, displ=TS_R2_4,  name='TS4')



    print(js)
    # print(js.Nodes)
    # print(js.Connectivity)
    # print(js.ElemProps)
    # print(js.Modes)
    # print(js.TimeSeries)
    js.write('_Rods.json')


if __name__ == '__main__':
    main()
if __name__ == '__test__':
    main(verbose=False, Large=False)
    try:
        os.delete('_Rods.json')
        os.delete('_Dummy.json')
    except:
        pass
