import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local 
from welib.stoch.distribution import *


x=np.linspace(-10,10,301)

n = gaussian_pdf(x, mu=0, sig=2)
r = rayleigh_pdf(x, sig=2)
u = uniform_pdf(x, a=-3, b=3)
g = gamma_pdf(x, alpha=1, beta=2)
e = exponential_pdf(x, beta=1)
w = weibull_pdf(x, k=2, A=4, x0=0)
l = lognormal_pdf(x, mu=0, sig=2) 


fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
ax.plot(x, n    , label='Normal')
ax.plot(x, r    , label='Rayleigh')
ax.plot(x, u    , label='Uniform')
ax.plot(x, g    , label='Gamma')
ax.plot(x, e    , label='Exponential')
ax.plot(x, w    , label='Weibull')
ax.plot(x, l    , label='Lognormal')
ax.set_xlabel('x')
ax.set_ylabel('Probability density function')
ax.legend()
ax.set_xlim([-1,10])
plt.show()



if __name__ == '__main__':
    pass
