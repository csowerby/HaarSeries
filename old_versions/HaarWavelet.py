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




import os
import numpy as np
import matplotlib.pyplot as plt
import imageio



#----------------------------------
""" Math Functions """ 



def gauss(x, a=1, sigma=1, mu=0):
    return a * np.exp(-(x - mu)**2/(2 * sigma**2))

def quadratic(t):
    return np.piecewise(t, [(t < -1) | (t > 1), (-1 <= t) & (t <= 1)], [0, lambda x: 1 - x**2])

def randomSignal(t):
    return np.sin(t) + np.random.normal(x, scale = 0.1, size = len(t))


sinc = lambda x: np.sinc(x)
saw = lambda x: scipy.signal.sawtooth(10 * x)

#---------------------------------------


#------------------------
""" Haar Functions """ 
def HaarMother(t):
    temp = np.piecewise(t, [(0 <= t) & (t < 1/2), (1/2 <= t) & (t < 1), (t < 0) | (t >= 1)], [-1, 1, 0])
    return temp


def haarScaling(t, k):
    return np.piecewise(t, [(k <= t) & (t < k+1), (t < k) | (t >= k+1)], [1, 0])

def Haar(t, n, k):
        return 2**(n/2) * HaarMother(2**n * t - k)

def haarFunctionProd(x, f, n, k):
    return f(x) * Haar(x, n, k)

def haarCoeff(f, n, k):
    coeff = quad(haarFunctionProd, k/2**n, (k+1)/2**n, args=(f, n, k), limit=500)
    return coeff[0]

def scalingCoeff(f, k):
    coeff = quad(f, k, k+1, limit=500)
    return coeff[0]


def haarSeries(x, f, n, k):
    
    func_sum = 0

    LEGACY FOR DOING ENTIRE N LEVELS AT A TIME 
    print("Running for all k values")
        
    print("Zero Frequency Component...")
    for kk in range(-3, 3):
        coeff = scalingCoeff(f, kk)
        func_sum = func_sum + coeff * haarScaling(x, kk)
        
    for nn in range(0, n):
        print(nn)
        for kk in range(- ((3 * 2**nn)) , 3 * 2**nn ):
            coeff = haarCoeff(f, nn, kk )
            func_sum = func_sum + coeff * Haar(x, nn, kk)
            
    else: 
    """
    for m in range(-1, n):
        if m == -1: 
            # Do zero frequency shit 
            for kk in range(-3, 3):
                coeff = scalingCoeff(f, kk)
                func_sum = func_sum + coeff * haarScaling(x, kk)
        else: 
            # Do non-zero frequency shit 
            for kk in range(-3 * 2**m, 3 * 2**m):
                coeff = haarCoeff(f, m, kk)
                func_sum = func_sum + coeff * Haar(x, m, kk)
        
        # At this stage do the partial sum at the n level up until k 
        
        if n == -1:
            for kk in range(-3, kk+1):
                coeff = scalingCoeff(f, kk)
                func_sum = func_sum + coeff * haarScaling(x, kk)
        else: 
            for kk in range(-3* 2**n, k+1):
                coeff = haarCoeff(f, n, kk)
                func_sum = func_sum + coeff * Haar(x, n, kk)
        
        """
    return func_sum
#-----------------------------------


def createGIF(xvalues, func, title:str, numIterations:int):
    "Writes images to current directory and compiles them into gifs"
    
    n_frames_start = 10
    
    filenames = []
    
    ## Plot the original function 

    
    plt.plot(x, func(x), "g")
    plt.title(f"Approximation of {title} with Haar Wavelet Series")
    plt.annotate("Original Function", xy=(0.97, 0.97), xycoords="axes fraction", verticalalignment='top', horizontalalignment='right')
    plt.xlim((-3, 3))
    plt.ylim((min(func(x)),1))
    filename = "original.png"
    plt.savefig(filename)
    
    plt.show()
    
    for i in range(n_frames_start):
        filenames.append(filenames)
    
    
    ## Plot Zero Frequency Components
    max_translation = 3
    
    for k in range(-max_translation, max_translation):
        pass
        
    
    ## Plot Sucessive approximations
    for n in range(-1, numIterations+1):
        for k in range(-3 * 2**n, 3 * 2**n):
            print(f'running Haar series for n = {n} and k = {k}')
            plt.plot(x, haarSeries(x, func, n, k))
            plt.xlim((-3, 3))
            plt.ylim((min(func(x)), 1))
            plt.show()
    



if __name__ == "__main__":
    """ Main Routine""" 
    
    x = np.linspace(-3, 3, 10000)
    numIterations = 5
    func = gauss
    title = "Gaussian"
    
    
    
    
    createGIF(x, func, title, numIterations)




"""


plt.plot(x, func(x))
plt.title(f"Haar Series of {funcName} Function")
plt.show()





for i in range(numIterations):
    
    plt.plot(x, haarSeries(x, func, i))
    plt.title(f'Haar Series of {funcName} Function') 
    plt.annotate(f'n={i} ',
                     xy=(0.97, 0.97), xycoords='axes fraction', bbox=dict(boxstyle='round', facecolor='white', alpha=0.3),
                     verticalalignment='top', horizontalalignment='right')

    #Create File
    filename = f'{i}.png'
    for j in range(10):
        filenames.append(filename)
    if  i == numIterations - 1:
        for j in range(20):
            filenames.append(filename)
        
    
    plt.savefig(filename)
    plt.show()
    plt.close()

print(filenames)
with imageio.get_writer(f"{funcName}.gif", mode="I") as writer: 
    for filename in filenames:
        image = imageio.imread(filename)
        writer.append_data(image)
        
#Remove files 
for filename in set(filenames):
    os.remove(filename)
    pass

#




# list comprehension 
def myFunc(x:int) -> int:
    return x + 5


#plt.plot(x, Haar(x, 0, 0),"k-")


"""




