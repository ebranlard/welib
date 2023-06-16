
from matplotlib.image import imread
import matplotlib.pyplot as plt
import numpy as np
import os
from welib.system.chaos import generatePNG # TODO find a better one
from welib.tools.svd import *
from welib.tools.regression import *

# --- Generate and read an image
tmpFile = os.path.join(os.path.dirname(__file__),'_lorenz.png')
generatePNG(tmpFile)
plt.close('all')

A = imread(tmpFile)
X = np.mean(A, -1); # Convert RGB to grayscale

plt.figure()
img = plt.imshow(X)
img.set_cmap('gray')
plt.axis('off')

U, S, VT = np.linalg.svd(X, full_matrices = False)
j = 0
for r in (5, 100, X.shape[1]):
    # Construct approximate image
    Xr = truncEvaluate(U, S, VT, r)
    #Xr = U[:,:r] @ S[0:r,:r] @ VT[:r,:]
    plt.figure()
    img = plt.imshow(Xr)
    img.set_cmap('gray')
    plt.axis('off')
    plt.title('r = ' + str(r))

plotS(S)


# --- Example 2 Linear Regression
x = 3 # True slope
a = np.arange(-2,2,0.25)
a = a.reshape(-1, 1)
b = x*a + np.random.randn(*a.shape) # Add noise
xtilde1 = regression(a, b, method='svd')
xtilde2 = regression(a, b, method='pinv')
plt.plot(a, x*a, color='k', linewidth=2, label='True line') # True relationship
plt.plot(a, b, 'x', color='r', ms = 10, label='Noisy data') # Noisy measurements
plt.plot(a,xtilde1 * a,'--',color='b',linewidth=2, label='Regression line1')
plt.plot(a,xtilde2 * a,'x', color='b',linewidth=2, label='Regression line2')
plt.xlabel('a')
plt.ylabel('b')
plt.legend()



if __name__ == '__main__':
    plt.show()

