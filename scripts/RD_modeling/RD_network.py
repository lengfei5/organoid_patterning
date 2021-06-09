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

from RD_network_functions import *
from itertools import permutations


#%% specify parameters
ID = '/Users/jiwang/workspace/imp/organoid_patterning/results/RD_topology_test/3N_testExample_python'

n = 3 # nb of node
k_length = 15 # nb of reaction parameters

states = 0 # steady state
turing = 0 # patterning 

# define ODE model without diffusion
def f_ode(x, t, k):
    #dRdt = np.empty(n)
    
    ## the test example from Zheng et al, 2016, Fig.S1A
    #dx0dt = 55.14*x[2]**2/(18.48**2 + x[2]**2) + 0.1 - 1.341*x[0]
    #dx1dt = 29.28*x[2]**2/(15.90**2 + x[2]**2) + 0.1 - 0.3508*x[1]
    #dx2dt = 16.17*x[0]**2/(x[0]**2 + 0.6421**2)*1.316**2/(1.316**2 + x[1]**2) + 0.1 - 1.203*x[2]
    
    ## the test example from Zheng et al, 2016, Fig.S1B (something not right in the formula)
    dx0dt = 50.86*x[0]**2/(x[0]**2 + 0.02315**2)*17.64**2/(17.64**2 + x[1]**2) + 0.1 - 0.09367*x[0]
    #dx0dt = 50.86*(x[0]**2/(x[0]**2 + 0.02315**2) * x[1]**2/(17.64**2 + x[1]**2)) + 0.1 - 0.09367*x[0]
    dx1dt = 17.43*x[0]**2/(x[0]**2 + 5.230**2)*1.038**2/(1.038**2 + x[2]**2) + 0.1 - 2.699*x[1]
    #dx1dt = 17.43*(5.230**2/(x[0]**2 + 5.230**2) * 1.038**2/(1.038**2 + x[2]**2)) + 0.1 - 2.699*x[1]
    dx2dt = 69.57*x[2]**2/(x[2]**2 + 1.000**2)*0.02100**2/(0.02100**2 + x[1]**2) + 0.1 - 0.1503*x[2]
    
    #dx0dt = k[0]*(1.0/(1.0 + pow(k[9]/x[0], 2))*1.0/(1.0 + pow(x[1]/k[10], 2))) + k[3] - k[6]*x[0]
    #dx1dt = k[1]*(1.0/(1.0 + pow(k[11]/x[0], 2)) + 1.0/(1.0 + pow(k[12]/x[2], 2))) + k[4] - k[7]*x[1]
    #dx2dt = k[2]*(1.0/(1.0 + pow(k[13]/x[0], 2))*1.0/(1.0 + pow(x[1]/k[14], 2))) + k[5] - k[8]*x[2]
    
    dRdt = [dx0dt, dx1dt, dx2dt]
    
    return dRdt

#%% speicify k parameters, difusor and d
k = np.ones(k_length)
k[0], k[1], k[2] = 10, 10, 20 
k[3], k[4], k[5] = 5, 5, 10
#k[6], k[7], k[8] = math.log(2)/10, math.log(2)/5, math.log(2)/10
k[6], k[7], k[8] = 0.1, 1, 0.1
k[9], k[10], k[11], k[12], k[13], k[14] = 30, 50, 50, 70, 40, 80

k 

binary_diffusor = [0, 1, 1]
d = [0, 2, 0.5]

# define initial conditions for ODE
x_max = 20001
x_int = 10000
c_init = permutations([1, x_int, x_max]) 

for i in list(c_init): 
    print (i) 


# %% find the steady state by integration
x0 = np.random.random(1) * np.ones(n)
x0 = [2.3, 0.4, 1.3]
#x0 = [0.0975, 0.0975, 0.0975]
err_tole = 0.0000001

t_final = 1000
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

# tspan = np.linspace(0, 5, 100)
# yinit = [0, -3]
# c = [4, 3, -2, 0.5]
# # Solve differential equation
# sol = solve_ivp(lambda t, y: f(t, y, c), 

#                 [tspan[0], tspan[-1]], yinit, t_eval=tspan, rtol = 1e-5)
# # Plot states
# state_plotter(sol.t, sol.y, 1)

#%% double check if steady state is reache by considering the first 100 time points as BurnIn
ss = np.zeros(n)
ss_fluc = np.zeros(n)
nb_passThreshold = 0

for i in range(n):
    ss[i] = np.mean(sol[100:149, i])
    ss_fluc[i] = np.abs(np.mean(sol[150:199, i]) - np.mean(sol[100:149, i]))
    if ss_fluc[i] >= err_tole:
        nb_passThreshold = nb_passThreshold + 1

while nb_passThreshold > 0:
    t_final = t_final*2
    t = np.linspace(0,t_final,200)
    sol = odeint(f_ode, x0, t,args=(k,))
    ss = np.zeros(n)
    ss_fluc = np.zeros(n)
    nb_passThreshold = 0
    for i in range(n):
        ss[i] = np.mean(sol[100:149, i])
        ss_fluc[i] = np.abs(np.mean(sol[150:199, i]) - np.mean(sol[100:149, i]))
        if ss_fluc[i] >= err_tole:
            nb_passThreshold = nb_passThreshold + 1

# check if ss is negative or imaginary solution
if any(ss <= 0) or any([isinstance(j, complex) for j in ss]):
    state = 2 

# continue to check if multiple steady states exist (to do)

#%% calculate the eigenvalue of Jacobian matrix without and with diffusion matrix
# some codes from https://www.sympy.org/scipy-2017-codegen-tutorial/notebooks/20-ordinary-differential-equations.html were 
# very helpful
import sympy as sym
X = sym.symbols(('x0:3'))
K = sym.symbols(('k0:15'))

f_sym = sym.Matrix(f_ode(X, None, K))
J = f_sym.jacobian(X)
J_func = sym.lambdify((X, K),  J)

#J_inputs = J.free_symbols
S = J_func(ss, k)
w  =  np.linalg.eigvals(S)
max(w), S


#%% eigenvalue computation with diffusion matrix to test if Turing instability 
# disperion relation plot for specific set of parameters
diffusing_nodes = binary_diffusor

#d = [0.07018, 1, 0.01057]
d = [1.00, 1.303, 0]
q = np.linspace(0, 5, 50)
lam_real = np.empty_like(q)
lam_im = np.empty_like(q)

for j in range(len(q)):
    #j = 1
    S2 = S - np.diag(np.multiply(d, q[j]**2))
    #wk,vk =  np.linalg.eig(S2)
    wk = np.linalg.eigvals(S2)
    lam_real[j] = wk.real.max()
    lam_im[j] = wk.imag[np.argmax(wk.real)]
    
plt.plot(q, lam_real)
#plt.axis([0, max(q), -1, 1])
plt.axhline(y=0, color='r', linestyle='-')
plt.show()

max(lam_real)

index_max = np.argmax(lam_real) 
lam_real_max = lam_real[index_max]
lam_im_max = lam_im[index_max]
q_max = q[index_max]

#%% define turing pattern types according to the eigenvalues (to finish)
# matlab code from Scholes how to distinguish type I and II
if max(Eig_save) > threshold % Check if any positive real eigenwert exists 
        if Eig_save(length(Eig_save)) < max(Eig_save)-0.01*max(Eig_save)
            %Check if Type I
            w = 1;
        else
            %This is Type II
            w = 2;
        end
    else 
        %No Turing instability
        w = 0;
    end

if lam_real_max < 0:
    turing = 0 
else:
    if np.abs(lam_im_max) > 0:
       turing = 1
    else:
        if lam_real_max > 0 and k_max < 0.000001:
            turing = 3
        else: 
            if lam_real_max > 0 and k_max > 0.000001:
                turing = 4
        

# def dispersion_relation(k_vals, d, mu):
#     lam = np.empty_like(k_vals)
#     for i, k in enumerate(k_vals):
#         A = np.array([[1-d*k**2,          1],
#                       [-2*mu,    -mu - k**2]])
#         lam[i] = np.linalg.eigvals(A).real.max()

#     return lam
# d = 0.05
# mu = 1
# k = np.linspace(0, 10, 200)
# lam_max_real_part = dispersion_relation(k, d, mu)
# plt.plot(k, lam_max_real_part)
# plt.axis([0, 10, -5, 5])
# plt.axhline(y=0, color='r', linestyle='-')
# plt.show()

