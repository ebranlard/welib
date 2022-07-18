import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation

from scipy.integrate import  solve_ivp #odeint


from welib.yams.kinetics import Fspring3D, Espring3D


GRAVITATIONAL_CONSTANT = 6.6743e-11 # [m3 kg-1 s-2] [Nm^2/kg^2]
GRAVITY                = 9.81 


class Part():
    def __init__(self, m, r0, v0=None):
        if v0 is None:
            v0=np.zeros(3)
        self.r0 = r0
        self.v0 = v0
        self.m  = m

    @property
    def r0(self):
        return self._r0

    @r0.setter
    def r0(self,r0):
        r0 = np.asarray(r0).flatten()
        self._r0 = r0
        self.r   = r0.copy()

    @property
    def v0(self):
        return self._v0

    @v0.setter
    def v0(self,v0):
        v0 = np.asarray(v0).flatten()
        self._v0 = v0
        self.v   = v0.copy()

    def state(self):
        return np.concatenate((self.r, self.v))

    def update(self, r, v):
        self.r=np.asarray(r).flatten()
        self.v=np.asarray(v).flatten()

    def gravitationalForce(self, p2, G=GRAVITATIONAL_CONSTANT):
        """ Force exerted by point 2 on current particle"""
        r1 = np.asarray(self.r)
        r2 = np.asarray(p2.r  )
        m1 = self.m
        m2 = p2.m
        dr = r1-r2
        r  = np.linalg.norm(dr)
        F2_on_self = - G*m1*m2*dr/r**3
        return F2_on_self

    def gravitationalEnergy(self, p2, G=GRAVITATIONAL_CONSTANT):
        r1 = np.asarray(self.r)
        r2 = np.asarray(p2.r  )
        m1 = self.m
        m2 = p2.m
        dr = r1-r2
        r  = np.linalg.norm(dr)
        return  - G*m1*m2/r**2

    def gravityForce(self, g=[0,0,GRAVITY]):
        return self.m*np.asarray(g)

class Point():
    def __init__(self, r0):
        self.r  = np.asarray(r0).flatten()
        self.m = 0
        self.v = np.array([0,0,0])

    def update(self, r, v):
        pass

    @property
    def r0(self):
        return self.r

    @property
    def v0(self):
        return self.v

class PartSystem(list):
    def __init__(self, elements=None, activeForces=None):
        list.__init__(self, elements)

        # --- Switch for forces
        self.activeForces=activeForces

        # --- Constant for forces
        self._gravityVect       = np.array([0,0,-GRAVITY])
        self.G                  = GRAVITATIONAL_CONSTANT

        # --- Connectivity
        self._springs=[]
        self._fixed=[]

    def __repr__(self):
        s='{} object\n'.format(type(self))
        s+=' - _gravityForce      : {}\n'.format(self._gravityForce)
        s+=' - _gravitationalForce: {}\n'.format(self._gravitationalForce)
        s+=' - g: {}\n'.format(self.g)
        s+=' - G: {}\n'.format(self.G)
        return s

    @property
    def n(self):
        return len(self)   

    @property
    def q0(self):
        """ Initial state vector (all positions and then all velocities)"""
        n = self.n
        q0 = np.zeros((n*6))
        ioff=n*3
        for ip,p in enumerate(self):
            q0[     ip*3:ip*3+3]      = p.r0
            q0[ioff+ip*3:ioff+ip*3+3] = p.v0
        return q0

    @property
    def q(self):
        """ Initial state vector (all positions and then all velocities)"""
        n = self.n
        q = np.zeros((n*6))
        ioff=n*3
        for ip,p in enumerate(self):
            q[     ip*3:ip*3+3]      = p.r
            q[ioff+ip*3:ioff+ip*3+3] = p.v
        return q

    def update(self, q):
        """ Update kinematics """
        n = int(len(q)/2) # Number of degrees of freedom
        x = q[:n]
        v = q[n:]
        for ip,p in enumerate(self):
            x_p = x[ip*3:ip*3+3]
            v_p = v[ip*3:ip*3+3]
            p.update(x_p,v_p)

    def forces(self, t=None):
        """ Compute """
        F = np.zeros(self.n*3) 
        # --- Gravitational force contribution
        if self._gravitationalForce:
            for ip,pi in enumerate(self):
                F_i=np.zeros(3)
                for jp,pj in enumerate(self):
                    if jp!=ip:
                        F_i +=pi.gravitationalForce(pj, G=self.G)
                F[ip*3:ip*3+3] += F_i
        # --- Gravity force contribution
        if self._gravityForce:
            for ip,pi in enumerate(self):
                F[ip*3:ip*3+3] += pi.gravityForce(g=self._gravityVect)
        # --- User field force contribution
        if self._fieldForce:
            for ip,pi in enumerate(self):
                F[ip*3:ip*3+3] += self._fieldForceFunction(t, pi)
        # --- Springs
        if self._springForce:
            for sp in self._springs:
                ip1 = sp['ip1']
                ip2 = sp['ip2']
                Fp = Fspring3D(self[ip1].r, self[ip2].r, sp['k'], sp['l0'])
                F[ip1*3:ip1*3+3] -= Fp
                F[ip2*3:ip2*3+3] += Fp

        return F

    def energies(self, res=None):
        """ Compute kinetic energy and potentials from misc forces when possible """
        E = {}
        E['kinetic'] = self.kineticEnergy(res)
        if self._gravitationalForce:
            E['gravitational'] = self.gravitationalEnergy(res)
        if self._gravityForce:
            E['gravity']       = self.gravityEnergy(res)
        if self._fieldForce:
            print('>>> TODO field force potential energy when possible')
        if self._springForce:
            E['spring'] = self.springEnergy(res)

        Etot = np.zeros(len(res.t))
        for _,Ek in E.items():
            Etot += Ek
        E['total'] = Etot
        return E

    def acc(self, t=None, F=None):
        if F is None:
            F = self.forces(t=t)
        a = np.zeros(F.shape)
        for ip,pi in enumerate(self):
            if isinstance(pi, Point):
                a[ip*3:ip*3+3] = 0 # Points have no mass
            else:
                a[ip*3:ip*3+3] = 1/pi.m * F[ip*3:ip*3+3]
        return a

    def dqdt(self, t, q):
        # Integrate system where states are position velocity
        #
        n=int(len(q)/2) # Number of degrees of freedom
        qdot = np.zeros(q.shape)
        v = q[n:]
        # Update kinematics
        self.update(q)
        # Kinetics
        a = self.acc(t)
        # Constraints (not pretty)
        v, a = self.kinConstraints(v,a)
        # 
        qdot[:n] = v
        qdot[n:] = a
        return qdot

    def kinConstraints(self, v,a):
        # Fixed BC
        for ip in self._fixed:
            v[ip*3:ip*3+3] = 0
            a[ip*3:ip*3+3] = 0
        return v, a


    def integrate(self, t, dqdt=None, method='RK45'):
        """
        methods: RK45, RK23, DOP852, LSODA
        """
        if dqdt is None:
            dqdt = self.dqdt
        res = solve_ivp(fun = dqdt, t_span=[t[0], t[-1]], y0 = self.q0, t_eval = t, method=method)
        return res


    def animate2Dtrajectories(self, res, XLIM=[-6,6], YLIM=[-6,6], plane='XY'):
        # Create figure for plotting
        fig, ax = plt.subplots()
        if plane=='XY':
            i1,i2=0,1
        elif plane=='XZ':
            i1,i2=0,2
        def animate(i):
            ax.clear()
            for ip,p in enumerate(self):
                x_p  = res.y[ip*3:ip*3+3,i]
                ax.plot(x_p[i1], x_p[i2], 'o')
            for sp in self._springs:
                ip1 = sp['ip1']
                ip2 = sp['ip2']
                x_p1 = res.y[ip1*3:ip1*3+3,i]
                x_p2 = res.y[ip2*3:ip2*3+3,i]
                ax.plot([x_p1[i1],x_p2[i1]], [x_p1[i2], x_p2[i2]],'k-')

            ax.set_xlim(XLIM)
            ax.set_ylim(YLIM)
            ax.set_title('{:.2f}'.format(res.t[i]))
        ani = animation.FuncAnimation(fig, animate, interval=1, frames=len(res.t))
        plt.show()

    def plot2Dtrajectories(self, res, XLIM=[-6,6], YLIM=[-6,6], plane='XY'):
        fig, ax = plt.subplots()
        if plane=='XY':
            i1,i2=0,1
        elif plane=='XZ':
            i1,i2=0,2
        for ip,p in enumerate(self):
            x_p  = res.y[ip*3:ip*3+3,:]
            ax.plot(x_p[i1,:], x_p[i2,:], '-')

        ax.set_xlim(XLIM)
        ax.set_ylim(YLIM)
        ax.set_title('')
        return fig, ax

    def plotTrajectories(self, res, comp=0):
        fig, ax = plt.subplots()
        for ip,p in enumerate(self):
            x_p  = res.y[ip*3:ip*3+3,:]
            ax.plot(res.t, x_p[comp,:], '-', label='Particle {}'.format(ip+1))
        ax.set_xlabel('Time [s]')
        ax.set_ylabel(['x','y','z'][comp] + ' [m]')

        return fig, ax

    def plotEnergy(self, res):
        E = self.energies(res)
        fig, ax = plt.subplots()
        vls=['-','--',':','-','.-',':']
        for ik,(k,v) in enumerate(E.items()):
            if k.lower()=='total':
                ax.plot(res.t, E[k], 'k-', label=k)
            else:
                ax.plot(res.t, E[k], ls=vls[ik], label=k)
        ax.legend()
        ax.set_xlabel('Time [s]')
        ax.set_ylabel('Energy [J]')

        return E, fig, ax


    # --- Average over particles
    def barycenter(self, res=None):
        if res is None:
            R = np.zeros(3)
            M = self.m
            for p in self:
                R[:] += p.r * (p.m/M)
        else:
            R = np.zeros((3,len(res.t)))
            M = self.m
            for ip,p in enumerate(self):
                x_p  = res.y[ip*3:ip*3+3,:]
                R[0:3,:] += x_p * (p.m/M)
        return R

    @property
    def m(self):
        return np.sum([p.m for p in self])

    def kineticEnergy(self, res=None):
        """ Return kinetic energy of all particles: 1/2 m v^2 """
        n = self.n
        iOff= n*3 #  three positions per particle
        if res is None:
            E = 0
            for ip,p in enumerate(self):
                v_p = p.v
                E += 1/2 * p.m * (v_p[0]**2 + v_p[1]**2 +v_p[2]**2)
        else:
            E = np.zeros(len(res.t))
            for ip,p in enumerate(self):
                v_p  = res.y[iOff+ ip*3:iOff+ ip*3+3,:]
                E += 1/2 * p.m * (v_p[0,:]**2 + v_p[1,:]**2 +v_p[2,:]**2)
        return E

    def springEnergy(self, res=None):
        """ Return potential energy of all springs: 1/2 k (l-l0)^2 """
        n = self.n
        if res is None:
            E = 0
            for sp in self._springs:
                ip1 = sp['ip1']
                ip2 = sp['ip2']
                E += Espring3D(self[ip1].r, self[ip2].r, sp['k'], sp['l0'])
        else:
            E = np.zeros(len(res.t))
            for sp in self._springs:
                ip1 = sp['ip1']
                ip2 = sp['ip2']
                x_p1 = res.y[ip1*3:ip1*3+3,:]
                x_p2 = res.y[ip2*3:ip2*3+3,:]
                l = np.sqrt((x_p2[0,:]-x_p1[0,:])**2 + (x_p2[1,:]-x_p1[1,:])**2 + (x_p2[2,:]-x_p1[2,:])**2)
                E += 1/2 * sp['k'] * (l-sp['l0'])**2
        return E

    def gravityEnergy(self, res=None):
        """ - m g.(r-r0) """
        gVect = self._gravityVect
        r0 = self[0].r0 # Reference
        if res is None:
            E = 0
            for ip,p in enumerate(self):
                dr = p.r - p.r0
                E += - p.m * g.dot(dr)
        else:
            E = np.zeros(len(res.t))
            for ip,p in enumerate(self):
                dr = res.y[ip*3:ip*3+3,:].copy()
                dr[0,:] -= p.r0[0]
                dr[1,:] -= p.r0[1]
                dr[2,:] -= p.r0[2]
                E += - p.m * gVect.dot(dr)
        return E

    def gravitationalEnergy(self, res=None):
        """ - G m1 m2 /r """
        G = self.G
        if res is None:
            E = 0
            for ip,pi in enumerate(self):
                for jp,pj in enumerate(self):
                    if jp!=ip:
                        E += pi.gravitationalEnergy(pj, G=self.G)
        else:
            E = np.zeros(len(res.t))
            n = self.n
            for ip,pi in enumerate(self):
                r_i = res.y[ip*3:ip*3+3,:]
                m_i = pi.m
                for jp,pj in enumerate(self):
                    if jp>ip:
                        r_j = res.y[jp*3:jp*3+3,:]
                        m_j = pj.m
                        r = np.sqrt((r_j[0,:]-r_i[0,:])**2 + (r_j[1,:]-r_i[1,:])**2 + (r_j[2,:]-r_i[2,:])**2)
                        E += - G*m_i*m_j/r
        return E


    # --- Connectivity / Constraints / BC
    def connect(self, p1, p2, connType, **kwargs):
        if connType=='spring':
            k = kwargs['k']
            if 'l0' not in kwargs.keys():
                l0 = np.linalg.norm(p2.r-p1.r)
            else:
                l0 = kwargs['l0']
            ip1 = self.index(p1) # NOTE: assume that indices of particle do not change
            ip2 = self.index(p2)
            d = {'p1':p1, 'p2':p2, 'ip1':ip1, 'ip2':ip2, 'k':k, 'l0':l0}
            self._springs.append(d)
        else:
            raise NotImplementedError()

    def fix(self, p1):
        ip1 = self.index(p1)
        self._fixed.append(ip1)

    # --- Force related
    @property
    def activeForces(self):
        return self._activeForces

    @activeForces.setter
    def activeForces(self, forces=None):
        """ set which forces are active, as a list of strings"""
        if forces is None:
            forces =['gravitational']
        self._activeForces=forces
        # Triggers
        self._gravitationalForce = 'gravitational' in forces
        self._gravityForce       = 'gravity' in forces
        self._fieldForce         = 'field' in forces
        self._springForce        = 'spring' in forces
        for f in forces:
            if f not in ['gravitational','gravity','field','spring']:
                raise Exception('Unknown force {}'.format(f))
    @property
    def g(self):
        return np.linalg.norm(self._gravityVect)
    @g.setter
    def g(self, g):
        self._gravityVect=np.array([0,0,-g])

    def setFieldForce(self, f):
        """ 
        f is a function with interface f(t,p) 
            t a time
            p a particle
          return a 3-vector of force at particle p
        """
        try:
            F = f(0, self[0])
        except:
            raise Exception('Force field `f` provided does not compute at t=0 and first particle')
        if len(F)!=3:
            raise Exception('Force field `f` does not return a 3-vector')

        self.activeForces+=['field']
        self._fieldForceFunction=f


if __name__ == '__main__':
    m0=100000000
#     m0=10
    r=1
    p1 = Part(m0, (0,0,0), (0.0,0,0))
#     vor = np.sqrt(G*p1.m/r)*0
    p2 = Part(m0, (r,0,0), (0,0,0))
    p3 = Part(m0, (1,1.5,0), (0,0,0))
    sys = PartSystem([p1,p2,p3], activeForces=['gravitational'])
    #sys = PartSystem([p1])
    # sys._gravitationalForce = True
    # sys._gravityForce       = False
#     sys._gravitationalForce = False
#     sys._gravityForce       = True

    print('q0',sys.q0)
    t   = np.linspace(0,760,700)
    res = sys.integrate(t)
    print(res.y.shape)
    sys.animate2Dtrajectories(res)

    pass
