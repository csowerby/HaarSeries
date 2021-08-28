#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 13 15:35:59 2021

@author: charliesowerby
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 13 13:56:11 2021

@author: charliesowerby
"""


import numpy as np
import matplotlib.pyplot as plt 
from scipy.integrate import quad
import scipy.signal

from typing import Callable

import random

import os
import imageio

""" NOTES FOR FURTHER IMPROVEMENT 

- Create a class that does the wavelet transform so we only need to calculate the values once 
- 2d plot for CWT

"""



#----------------------------------
""" Math Functions """ 



def gauss(x, a=1, sigma=1, mu=0):
    return a * np.exp(-(x - mu)**2/(2 * sigma**2))

def quadratic(t):
    return np.piecewise(t, [(t < -1) | (t > 1), (-1 <= t) & (t <= 1)], [0, lambda x: 1 - x**2])

def randomPoly(t, order = 3):
    func = 0
    for i in range(order + 1):
        func = func + (random.random() - 0.5) * t**i
    return func
    
        


sinc = lambda x: np.sinc(2*x)
saw = lambda x: scipy.signal.sawtooth(10 * x)

triangle = lambda x: scipy.signal.sawtooth(2 *  np.pi * (x + 0.33 ), 0.5    )

def beat(t, a=3, b=5):
    return np.cos(3 * np.pi * t / 3) * np.sin(5 * np.pi * t / 3)
#---------------------------------------


#------------------------
""" Haar Functions """ 
def HaarMother(t):
    temp = np.piecewise(t, [(0 <= t) & (t < 1/2), (1/2 <= t) & (t < 1), (t < 0) | (t >= 1)], [-1, 1, 0])
    return temp


def haarScaling(t, k):
    """ function for the identity function on a unit interval """ 
    return np.piecewise(t, [(k <= t) & (t < k+1), (t < k) | (t >= k+1)], [1, 0])

def Haar(t, n, k):
    """ Haar function Haar(n, k) """ 
    return 2**(n/2) * HaarMother(2**n * t - k)

def haarFunctionProd(x, f, n, k):
    return f(x) * Haar(x, n, k)

def haarCoeff(f, n, k):
    coeff = quad(haarFunctionProd, k/2**n, (k+1)/2**n, args=(f, n, k), limit=50)
    return coeff[0]

def scalingCoeff(f, k):
    coeff = quad(f, k, k+1, limit=500)
    return coeff[0]

def Nth_Haar_sum(x:np.ndarray, f:Callable, n:int) -> Callable: 
    """ Returns the Nth partial sum (including all k values in domain) up to and including N """ 
        
    func_sum = 0
    #Zero frequency Component
    for kk in range(-3, 3):
        coeff = scalingCoeff(f, kk)
        func_sum = func_sum + coeff * haarScaling(x, kk)
    
    # frequency starting at (n=0) going to and including n=N 
    for nn in range(0, n+1):
        for kk in range(- ((3 * 2**nn)) , 3 * 2**nn ):
            coeff = haarCoeff(f, nn, kk )
            func_sum = func_sum + coeff * Haar(x, nn, kk)
    
    return func_sum 

def Kth_Haar_sum(x:np.ndarray, f:Callable, n:int, k:int) -> Callable: 
    """ Returns the Kth partial sum of Haar(n, k) from (-2^n, k+1) at the Nth level. """ 
    func_sum = 0 
    if n == -1:
        #Zero Frequency sum
        for kk in range(-3, k+1):
            coeff = scalingCoeff(f, kk)
            func_sum = func_sum + coeff * haarScaling(x, kk)
    else: 
        # Nth frequency sum
        for kk in range(-3*2**n, k+1):
            coeff = haarCoeff(f, n, kk)
            func_sum = func_sum + coeff* Haar(x, n, kk)
            
        func_sum = func_sum + Nth_Haar_sum(x, f, n-1)    
    
    return func_sum 
    

    



#-----------------------------------------
"GIF Creation for Final Product"



def createGIF(x, func, title:str, numIterations:int):
    "Writes images to current directory and compiles them into gifs"
    
    n_frames_start = 10
    
    filenames = []
    
    ## Plot the original function 
    """
    
    plt.plot(x, func(x), "g")
    plt.title(f"Approximation of {title} with Haar Wavelet Series")
    plt.annotate("Original Function", xy=(0.97, 0.97), xycoords="axes fraction", verticalalignment='top', horizontalalignment='right')
    plt.xlim((-3, 3))
    plt.ylim(1.2 * min(func(x)), 1.2 * max(func(x)))
    filename = "images/original.png"
    plt.savefig(filename)
    
    plt.show()
    
    for i in range(n_frames_start):
        filenames.append(filename)
    """
    max_translation = 3
    ## Plot Nth Partial Sum Approximation
    """
    for N in range(-1, numIterations):
        plt.plot(x, Nth_Haar_sum(x, func, N))
        plt.title(f"Approximation of {title} with Haar Wavelet Series")
        if N == -1:
            plt.annotate("Zero Frequency", xy=(0.97, 0.97), xycoords="axes fraction", verticalalignment='top', horizontalalignment='right')
        else:
            plt.annotate(f"n={N}", xy=(0.97, 0.97), xycoords="axes fraction", verticalalignment='top', horizontalalignment='right')
        plt.xlim((-3, 3))
        plt.ylim(1.2 * min(func(x)), 1.2 * max(func(x)))
        
        
        plt.show()

    """
    ## Plot N, Kth partial sum approximation 
    for N in range(-1, numIterations):
        k_range = 3 if N==-1 else 3 * 2**N
        for k in range(-k_range, k_range):
            print(f"Running n={N}, k={k}")
            
            fig, axs = plt.subplots(1, 2, figsize=(12, 5))
            fig.suptitle(f"Approximation of {title} Function with Haar Wavelet Series")
            axs[0].plot(x, func(x), label=f"{title.lower()}(x)")
            
            new_ax=axs[0].twinx()
            
            if N == -1: 
                temp_string = f"[{k}, {k+1})"
                new_ax.plot(x, haarScaling(x, k), "orange", label=rf'$1_{{{temp_string}}}$')
                axs[1].axvspan(k, k+1, alpha=0.5, color='yellow')
                axs[0].plot([], [], "orange", label=rf'$1_{{{temp_string}}}$' )
            else:
                new_ax.plot(x, Haar(x, N, k), "orange", label=rf'Haar$_{N}^{{{k}}}$(x)')
                axs[0].plot([], [], "orange", label=rf'Haar$_{N}^{{{k}}}$(x) ')
                axs[1].axvspan(k/2**N, (k+1)/2**N, alpha=0.5, color='yellow')
                
            

            axs[1].plot(x, Kth_Haar_sum(x, func, N, k), 'g', label = f'({N+1}, {k + k_range})-th Partial Sum')
            axs[0].legend(loc="upper right")

                
            new_ax.set_ylim(-2**((numIterations-1)/2)* 1.1 , 1.5 * 2**((numIterations-1)/2) )
            axs[0].set_ylim(-1.1 * max(func(x)), 1.5 * max(func(x)))
            axs[0].set_xlim(-3, 3)
            axs[1].set_xlim(-3, 3)
            axs[1].set_ylim(-1.1 * max (func(x)), 1.5 * max(func(x)))
            
            axs[1].legend(loc="upper right")     
            
            
            plt.subplots_adjust(wspace=0.2)
            
            filename = f"images/n_{N}_k_{k}.png"
            for i in range(max([1, 4-N])):
                filenames.append(filename)
            
            plt.savefig(filename)
            plt.close()
            
            """
            OLD SINGLE PLOT 
            plt.plot(x, Kth_Haar_sum(x, func, N, k))
            plt.title(f"Approximation of {title} with Haar Wavelet Series")
            if N == -1:
                plt.annotate("Zero Frequency", xy=(0.97, 0.97), xycoords="axes fraction", verticalalignment='top', horizontalalignment='left')
            else:
                plt.annotate(f"n={N}\nk={k}", xy=(0.97, 0.97), xycoords="axes fraction", verticalalignment='top', horizontalalignment='right')
            plt.xlim((-3, 3))
            plt.ylim(1.2 * min(func(x)), 1.2 * max(func(x)))
            
            filename = f"images/n_{N}_k_{k}.png"
            filenames.append(filename)
            plt.savefig(filename)
            plt.close()
            """
        
    ## Create the GIF 
    print("Creating GIF...")
    with imageio.get_writer(f"2part_gifs/{title}_{numIterations}iterated.gif", mode="I") as writer: 
        for filename in filenames:
            image = imageio.imread(filename)
            writer.append_data(image)
        
        #Remove files 
    for filename in set(filenames):
        os.remove(filename)
            



if __name__ == "__main__":
    """ Main Routine""" 
    
    x = np.linspace(-3, 3, 10000)
    numIterations = 6
    func = triangle
    title = "Triangle"
    
    
    
    
    createGIF(x, func, title, numIterations)
    print("Done!")

