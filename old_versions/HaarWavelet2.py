#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 12 16:21:25 2021

@author: charliesowerby
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 12 13:35:16 2021

@author: charliesowerby
"""

import numpy as np
import matplotlib.pyplot as plt 
from scipy.integrate import quad
import scipy.signal




def Haar(t, n, k):
    if n == 0:
        return np.piecewise(t, [(k <= t) & (t < k+1), (t < k) | (k+1 <= t)], [1, 0])
    else:
        return np.piecewise(t, [(k/2**n <= t) & (t < (k + 1/2)/2**n), ((k + 1/2)/2**n <= t) & (t < (k+1)/2**n), (t < (k)/2**n) | ((k+1)/2**n <= t)], [-2**(n/2), 2**(n/2), 0])

def gauss(x, a=1, sigma=1, mu=0):
    return a * np.exp(-(x - mu)**2/(2 * sigma**2))

def quadratic(t):
    return np.piecewise(t, [(t < -1) | (t > 1), (-1 <= t) & (t <= 1)], [0, lambda x: 1 - x**2])

def haarFunctionProd(x, f, n, k):
    return f(x) * Haar(x, n, k)


def haarCoeff(f, n, k):
    coeff = quad(haarFunctionProd, k/2**n, (k+1)/2**n, args=(f, n, k), limit=500)
    return coeff[0]


def haarSeries(x, f, n):
    sum = 0
    
    
    
    for nn in range(0, n+1):
        print(nn)
        for kk in range(- ((5 * 2**nn)+1) , 5 * 2**nn ):
            coeff = haarCoeff(f, nn, kk )
            sum = sum + coeff * Haar(x, nn, kk)
    return sum
    
x = np.linspace(-1, 1, 10000)
#plt.plot(x, Haar(x, 0, 0))



sinc = lambda x: np.sinc(x)
saw = lambda x: scipy.signal.sawtooth(10 * x)



func = gauss
numIterations = 5

annotation
plt.plot(x, func(x))
plt.plot(x, haarSeries(x, func, numIterations))
plt.title("Approximation of Gaussian with Haar Wavelets")
   plt.annotate(f' n =  ',
                     xy=(0.03, 0.97), xycoords='axes fraction', bbox=dict(boxstyle='round', facecolor='tan', alpha=0.3),
                     verticalalignment='top', horizontalalignment='left')
plt.show()

    


#



#plt.plot(x, Haar(x, 0, 0),"k-")







