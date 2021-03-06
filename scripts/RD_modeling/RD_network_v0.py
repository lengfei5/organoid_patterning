#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 28 16:14:06 2021
@author: jingkui.wang

The caltech course http://be150.caltech.edu/2020/content/lessons/20_turing.html
was much help

Main function to assess the RD stability

Inputs of function: 
    ODE model
    reaction parameters
    diffusion molecules (0, immobile; 1, diffusible)
    diffusion rate
    
Output of function for each of given parameters (reaction, diffusion):
    parameters
    steady state (could be multiple)
    stability of steady state without diffusion
    stability of steady state with diffusion
    
"""
import numpy as np
import pandas as pd
import numba
import scipy.integrate
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
from RD_network_functions import state_plotter

import matplotlib.pyplot as plt

import biocircuits
import bokeh.io
import bokeh.plotting

import panel as pn
pn.extension()

import math
import itertools
import sympy as sym
from itertools import permutations


from RD_network_functions import *

#ID = '/Users/jiwang/workspace/imp/organoid_patterning/results/RD_topology_test/3N_testExample_python'

#%% specify the network topology or model with parameters
n = 3 # nb of node
k_length = 15 # nb of reaction parameters: 3* number of nodes (3*3) + number of interactions (6)

# define symbolic variables for Jacobian matrix 
X = sym.symbols(('x0:3'))
K = sym.symbols(('k0:15'))

# define ODE model without diffusion
def f_ode(x, t, k):
    #dRdt = np.empty(n)
    
    ## the test example from Zheng et al, 2016, Fig.S1A
    #dx0dt = 55.14*x[2]**2/(18.48**2 + x[2]**2) + 0.1 - 1.341*x[0]
    #dx1dt = 29.28*x[2]**2/(15.90**2 + x[2]**2) + 0.1 - 0.3508*x[1]
    #dx2dt = 16.17*x[0]**2/(x[0]**2 + 0.6421**2)*1.316**2/(1.316**2 + x[1]**2) + 0.1 - 1.203*x[2]
    
    ## the test example from Zheng et al, 2016, Fig.S1B (something not right in the formula)
    #dx0dt = 50.86*x[0]**2/(x[0]**2 + 0.02315**2)*17.64**2/(17.64**2 + x[1]**2) + 0.1 - 0.09367*x[0]
    #dx0dt = 50.86*x[0]**2/(x[0]**2 + 0.02315**2)*x[1]**2/(17.64**2 + x[1]**2) + 0.1 - 0.09367*x[0]
    #dx1dt = 17.43*x[0]**2/(x[0]**2 + 5.230**2)*1.038**2/(1.038**2 + x[2]**2) + 0.1 - 2.699*x[1]
    #dx1dt = 17.43*(5.230**2/(x[0]**2 + 5.230**2) * 1.038**2/(1.038**2 + x[2]**2)) + 0.1 - 2.699*x[1]
    #dx2dt = 69.57*x[2]**2/(x[2]**2 + 1.000**2)*0.02100**2/(0.02100**2 + x[1]**2) + 0.1 - 0.1503*x[2]
    
    ## only non-competitive interactions were chosen for approximiation 
    ## NT organoid phase III pattern selection: FoxA2, Noggin and BMP
    dx0dt = k[0] - k[3]*x[0] + k[6]*(x[0]**2/(x[0]**2 + k[9]**2) * k[10]**2/(x[2]**2 + k[10]**2)) # Foxa2
    dx1dt = k[1] - k[4]*x[1] + k[7]*(x[0]**2/(x[0]**2 + k[11]**2) * x[2]**2/(x[2]**2 + k[12]**2)) # Noggin
    dx2dt = k[2] - k[5]*x[2] + k[8]*(x[0]**2/(x[0]**2 + k[13]**2) * k[14]**2/(x[1]**2 + k[14]**2)) # BMP
    
    dRdt = [dx0dt, dx1dt, dx2dt]
    
    return dRdt
    
        
#%% sampling the parameters, which node is difusor and diffusion coeffs
nb_sampling = 2
binary_diffusor = [0, 1, 1]

# paramtere vector k sampling
ks = np.logspace(-1, 2.0, num=nb_sampling)
gammas = np.logspace(-2, 1, num = nb_sampling)

k_grid1 = list(itertools.product(ks, ks, ks, gammas, gammas, gammas))
k_grid2 = list(itertools.product(ks, repeat = (k_length - 6)))
k_grid = list(itertools.product(k_grid1, k_grid2))
#k_grid = itertools.combinations_with_replacement(ks, repeat = k_length)
#k_grid = pd.DataFrame(k_grid)
#k_grid.columns = K
#gamma_grid = 

# diffusion rate sampling
d_range = np.logspace(-1, 2.0, num = 10)
d_grid = list(itertools.product(d_range, repeat=2))

# initial conditions: each node has 3 initial values
nb_init = 3
x_init = np.logspace(-1, 4, nb_init)
c_init = itertools.combinations_with_replacement(x_init, n)

# time 
t_final = 1000


#%% big loop over each k parameter vector and save the result for each sampled d combination
#for i in range(len(k_grid)):
for i in range(100):
    #i = 0 # test first k_grid
    i = 0
    import time
    start_time = time.process_time()
    #k = np.asarray(k_grid[i])
    #k0 = k_grid[i]
    k = [element for tupl in k_grid[i] for element in tupl]
    k = np.asarray(k)
    
    k[0] = 0.0256364845015646
    k[1] = 0.196353074934487
    k[2] = 0.01
    k[3] = 0.312886357461524
    k[4] = 0.058476504153979
    k[5] = 0.1
    k[6] = 2.89617848764822
    k[7] = 0.317400763608912
    k[8] = 1
    k[9] = 1
    k[10] = 1
    k[11] = 1
    k[12] = 7.65768501184393
    k[13] = 0.949509951832922
    k[14] = 100
    #k[3], k[4], k[5] = 0.3, 0.5, 0.4
    #k[6], k[7], k[8] = math.log(2)/10, math.log(2)/5, math.log(2)/10
    #k[6], k[7], k[8] = 30, 50, 20
    #k[9], k[10], k[11], k[12], k[13], k[14], k[15], k[16] = 14, 3, 0.2, 5, 10, 1, 2, 5
    
    
    # %% find the steady state by integration (initial guess)
    x0 = np.random.random(1) * np.ones(n)
    #x0 = [2.3, 0.4, 1.3]
    #x0 = [0.0975, 0.0975, 0.0975]
    
    t = np.linspace(0, t_final, 200)
    sol = odeint(f_ode, x0, t, args=(k,))
    
    ## check the integration solution
    fig,ax = plt.subplots()
    ax.plot(t,sol[:,0],label='x1')
    ax.plot(t,sol[:,1],label='x2')
    ax.plot(t,sol[:,2],label='x3')
    #ax.plot(t,result[:,2],label='R0=1')
    ax.legend()
    ax.set_xlabel('t')
    ax.set_ylabel('x')
    sol[199, ]
    
    # double check if steady state is reache by considering the first 100 time points as BurnIn
    ss0 = check_BurnIn_steadyState(sol, f_ode, k, n, x0, t_final)
    
    #%% check if multiple steady states exist and save
    # define initial conditions for ODE
    ss_saved = Multi_steadyStates(ss0, c_init, f_ode, k, n)
    
    for kk in range(len(ss_saved)):
        ss = ss_saved[kk]    
        if any(ss <= 0) or any([isinstance(j, complex) for j in ss]):
            ss_saved.remove(ss)
    
    #%% calculate the eigenvalue of Jacobian matrix without and with diffusion matrix
    # some codes from https://www.sympy.org/scipy-2017-codegen-tutorial/notebooks/20-ordinary-differential-equations.html were 
    # very helpful
    f_sym = sym.Matrix(f_ode(X, None, K))
    J = f_sym.jacobian(X)
    J_func = sym.lambdify((X, K),  J)
    
    #%% eigenvalue computation with diffusion matrix to test if Turing instability 
    # disperion relation plot for specific set of parameters
    #d = [1, 10.03, 0]
    q = 2*3.15169 / np.logspace(0, 3.0, num=50) # wavenumber 
    
    # loop over stedy states
    for ii in range(len(ss_saved)):
        # ii = 0
        ss = ss_saved[ii]
        #J_inputs = J.free_symbols
        S = J_func(ss, k)
        w  =  np.linalg.eigvals(S)
        w.real.max(), S
        
        for val in d_grid: # loop diffusion matrix for each steady state
            d = np.asarray(val)
            d = [0.0, d[0], d[1]]
            
            d = [0, 1, 1.4]
            lam_real = np.empty_like(q)
            lam_im = np.empty_like(q)
            
            # loop over the wavenumber 
            for j in range(len(q)):
                #j = 1
                S2 = S - np.diag(np.multiply(d, q[j]**2))
                #wk,vk =  np.linalg.eig(S2)
                wk = np.linalg.eigvals(S2)
                lam_real[j] = wk.real.max()
                lam_im[j] = wk.imag[np.argmax(wk.real)]
                
            plt.plot(q, lam_real)
            plt.show()
            #plt.axis([0, max(q), -1, 1])
            #plt.Axes.axhline(y=0, color='r', linestyle='-', xmin = np.min(lam_real), xmax = np.max(lam_real))
            print(max(lam_real))
            
            index_max = np.argmax(lam_real) 
            lam_real_max = lam_real[index_max]
            #lam_im_max = lam_im[index_max]
            #q_max = q[index_max]
            
            # save the result, k parameter, steady state, d parameters, q values, lambda_real, lambda imaginary
            if lam_real_max >= 0:
                # save  
                turing = 0 
    
                
    print(time.process_time() - start_time, "seconds")
    

def main():

    data = "My data read from the Web"

    print(data)

    modified_data = process_data(data)

    print(modified_data)


if __name__ == "__main__":

    main()


