import numpy as np

def skew(x):
    """ Returns the skew symmetric matrix M, such that: cross(x,v) = M v 
    [ 0, -z , y]
    [ z,  0 ,-x]
    [-y,  x , 0]
    """
    x=np.asarray(x).ravel()
    return np.array([[0, -x[2], x[1]],[x[2],0,-x[0]],[-x[1],x[0],0]])


def translateLoadsJacobian(JS, r0, FS0, method='matrix', sj0=None, j_is_k=False, undisplaced=True):
    """ 
    Transfer Jacobians of loads at a source point "S" to a destination point "D"
    assuming a rigid body motion between the source and destination point

    INPUTS:
    - JS: Jacobian of loads at the source point, 6x6 array dF_S/dx_S
          where F_S  = (f_x, f_y, f_z, m_x, m_y, m_z)_S
          where dx_S = (dx, dy, dz, dt_x, dt_y, dt_z)_S (displacement and small rotations)
    - r0 = rS-rD : 3-vector from the destination point to the source point
    - FS0: 3-vector of the operating point forces at the source node

    OUTPUTS:
    - JD: Jacobian of loads at the source point, 6x6 array dF_S/dx_S

    Reference:
       See Branlard, Jonkman, Brown to be published 2022

    """
    r0til  = skew(r0)
    FS0til = skew(FS0)
    if method=='matrix':
        I3 = np.eye(3)
        Z3 = np.zeros((3,3))
        if sj0 is not None:
            # Multiple points
            sj0til  = skew(sj0)
            T1 = np.block([ [ I3   ,  Z3 ],
                            [ r0til,  I3 ]])
            T2 = np.block([ [ I3   ,-sj0til ],
                            [ Z3   ,  I3 ]])



            if j_is_k:
                T3 = np.block([ [ Z3   ,  Z3 ],
                                [ Z3   ,  FS0til.dot(r0til) ] ]) # TODO
            else:
                T3 = np.block([ [ Z3   ,  Z3 ],
                                [ Z3   ,  Z3 ] ]) # TODO
            JD = T1.dot(JS.dot(T2)) + T3
            if undisplaced and j_is_k:
                JD[3:6,0:3] += -FS0til

        else:
            T1 = np.block([ [ I3   ,  Z3 ],
                            [ r0til,  I3 ]])
            T2 = np.block([ [ I3   ,-r0til ],
                            [ Z3   ,  I3 ]])
            T3 = np.block([ [ Z3   ,  Z3 ],
                            [ Z3   ,  FS0til.dot(r0til) ] ]) 
            JD = T1.dot(JS.dot(T2)) + T3
    else:
        # component by component (for debugging)
        dFdU  = JS[0:3, 0:3]
        dMdU  = JS[3:6, 0:3]
        dFdT  = JS[0:3, 3:6]
        dMdT  = JS[3:6, 3:6]
        if sj0 is not None:
            sj0til  = skew(sj0)
            dFdUD = dFdU
            dFdTD = dFdT - dFdU.dot(sj0til)
            dMdUD = dMdU + r0til.dot(dFdU) 
            if undisplaced and j_is_k:
                dMdUD += FS0til

            #print('dMdU',dMdU[1,0],(r0til.dot(dFdU))[1,0],dMdUD[1,0])
         
            #dMdTD = dMdT - dMdU.dot(r0til) + r0til.dot(dFdT) - r0til.dot(dFdU).dot(r0til)  + FS0til.dot(r0til)
            C1 = dMdT 
            C2 =      - dMdU.dot(sj0til) 
            C3 =                        + r0til.dot(dFdT)
            C4 =                                          - r0til.dot(dFdU.dot(sj0til))
            if j_is_k:
                C5 =                                                                        + FS0til.dot(r0til)
            else:
                C5=np.zeros((3,3))
            dMdTD = np.zeros((3,3))
            dMdTD +=  C1
            dMdTD +=  C2
            dMdTD +=  C3
            dMdTD +=  C4
            dMdTD +=  C5
        else:
            raise Exception('Temporary stopped')
            dFdUD = dFdU
            dFdTD = dFdT - dFdU.dot(r0til)
            dMdUD = dMdU + r0til.dot(dFdU) 
            #dMdTD = dMdT - dMdU.dot(r0til) + r0til.dot(dFdT) - r0til.dot(dFdU).dot(r0til)  + FS0til.dot(r0til)
            C1 = dMdT 
            C2 =      - dMdU.dot(r0til) 
            C3 =                        + r0til.dot(dFdT)
            C4 =                                          - r0til.dot(dFdU.dot(r0til))
            C5 =                                                                        + FS0til.dot(r0til)
            dMdTD = np.zeros((3,3))
            dMdTD +=  C1
            dMdTD +=  C2
            dMdTD +=  C3
            dMdTD +=  C4
            dMdTD +=  C5
            #     if rref is None:
            #         C5b = r0til.dot(FS0til) # NOTE REVERSED SIGN
            #     else:
            #         r1til  = skew(r0+np.asarray(rref))
            #         C5b = r1til.dot(FS0til) # NOTE REVERSED SIGN
            #     C5b[0,0]=0 # NOTE
            #     C5b[1,1]=0 # NOTE
            #     C5b[2,2]=0 # NOTE
            #     dMdTD +=  C5b
            # print('C1\n', C1)
            # print('C2\n', C2)
            # print('C3\n', C3)
            # print('C4\n', C4)
            # print('C5\n', C5)
            #     print('C52\n',C5b)
            #     print(r0, FS0)
        JD = np.block([ [dFdUD,  dFdTD],
                        [dMdUD , dMdTD]
                      ])
    return JD


def translateLoads(FS0, MS0, r0):
    FD0 = FS0
    MD0 = MS0 + skew(r0).dot(FS0)
    return FD0, MD0





# --------------------------------------------------------------------------------}
# --- 3D springs 
# --------------------------------------------------------------------------------{
def Fspring3D(r1,r2,k,l0):
    """ return linear spring force between two points """
    r1 = np.asarray(r1).flatten()
    r2 = np.asarray(r2).flatten()
    l  = np.sqrt((r2-r1).dot(r2-r1))
    dl = (l-l0)
    return -k * dl *(r2-r1)/l

def Espring3D(r1,r2,k,l0):
    """ return linear spring energy between two points """
    r1 = np.asarray(r1).flatten()
    r2 = np.asarray(r2).flatten()
    l  = np.sqrt((r2-r1).dot(r2-r1))
    dl = (l-l0)
    return 1/2 * k * dl**2

def Fdamp3D(r1,r2,r1d,r2d,c):
    """ return linear damping force between two points """
    r1  = np.asarray(r1).flatten()
    r2  = np.asarray(r2).flatten()
    r1d = np.asarray(r1d).flatten()
    r2d = np.asarray(r2d).flatten()
    l  = np.sqrt((r2-r1).dot(r2-r1))
    e  = (r2-r1)/l
    dv = (r2d-r1d).dot(e) 
    return -c * dv * e

def FspringDamp3D(r1, r2, r1d, r2d, k, l0, c):
    """ return linear spring/damper force between two points """
    r1 = np.asarray(r1).flatten()
    r2 = np.asarray(r2).flatten()
    r1d = np.asarray(r1d).flatten()
    r2d = np.asarray(r2d).flatten()
    l  = np.sqrt((r2-r1).dot(r2-r1))
    dl = (l-l0)
    e  = (r2-r1)/l
    dv = (r2d-r1d).dot(e)
    return - (k * dl  + c * dv) * e   

