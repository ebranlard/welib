from vortilib.elements.SourceEllipsoid import *
import numpy as np
import matplotlib.pyplot as plt

# Ellipse parameters
a = 3
b = 2
e = np.sqrt(a**2 - b**2)/a

# e = 1/a # HACK for now such that a*e=1
# b= a * np.sqrt(1-e**2)

k=a*e
print('a    ',a)
print('b    ',b)
print('e    ',e)
print('k    ',k)
print('zeta0',1/e)



# Coordinates of the ellipse
ne=15
zeta0 = np.ones(ne)*1/e
mu0   = np.linspace(-1,1,ne)
x0,y0,z0 = T_semielliptic_cart(mu0,zeta0,mu0*0,k)
x0e = np.linspace(-a,a,ne)
y0e=b*np.sqrt(1-(x0e/a)**2)

# mu02,zeta02,theta02 = T_cart_semielliptic(x0,y0,z0,k)
# print(theta02)


# Grid
nmu=41
nzeta=10
vMu   = np.cos(np.linspace(-np.pi,np.pi,nmu))   # full range [-1,1]
vZeta = np.cosh(np.linspace(0,2,nzeta))
# vMu   = np.linspace(-1,1,nmu)   # full range [-1,1]
# vZeta = np.unique(np.concatenate((np.linspace(1,zeta0,nzeta),np.linspace(zeta0,2*zeta0,nzeta))))
# vZeta = np.array([zeta0*1.2])

MU,ZETA = np.meshgrid(vMu, vZeta)
THETA = MU*0



X,Y,Z    = T_semielliptic_cart(MU,ZETA,THETA,k)
X1,Y1,Z1 = T_semielliptic_cart(MU,ZETA,THETA+np.pi,k)


MU2,ZETA2,THETA2 = T_cart_semielliptic(X,Y,Z,k)
print(np.min(np.abs(MU2-MU)))
print(np.max(np.abs(MU2-MU)))
print(np.min(np.abs(ZETA2-ZETA)))
print(np.max(np.abs(ZETA2-ZETA)))





plt.figure()
plt.plot(X  ,Y,'-')
# plt.plot(X1,Y1,'-')
plt.plot(X.T  ,Y.T,'-')
# plt.plot(X1.T,Y1.T,'-')
plt.plot(x0e,y0e ,'k-',lw=2)
plt.plot(x0e,-y0e,'k-',lw=2)
# plt.plot(x02,y02,'b:')
plt.plot(a,0,'kd')
plt.plot(0,b,'kd')
plt.plot(-e*a,0,'ko')
plt.plot(e*a,0,'ko')
plt.axis('equal')
plt.show()




