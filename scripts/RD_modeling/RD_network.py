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
import matplotlib.pyplot as plt

import biocircuits
import bokeh.io
import bokeh.plotting

import panel as pn
pn.extension()

from RD_network_functions import *

from itertools import permutations 

# specify parameters
ID = '/Users/jiwang/workspace/imp/organoid_patterning/results/RD_topology_test/2N_testExample_python'

n = 2 # nb of node
k_length = 9 # nb of reaction parameters

# define ODE model without diffusion
from scipy.integrate import solve_ivp
from RD_network_functions import state_plotter

def f_ode(x, t, k):
    #dRdt = np.empty(n)
    dx1dt = k[0]*pow(k[1], -2)*pow(x[0], 2)*pow(1.0+pow(k[0],-2)*pow(x[0],2)+pow(k[2], -2)*pow(x[1], 2), -1)+k[5]-k[7]*x[0]  
    dx2dt = k[3]*pow(k[4], -2)*pow(x[0], 2)*pow(1.0+pow(k[4],-2)*pow(x[0],2), -1) + k[6] - k[8] * x[1] 
    
    dRdt = [dx1dt, dx2dt]
    
    return dRdt

# speicify k parameters, difusor and d
k = np.ones(k_length)
binary_diffusor = [1, 1]
d = [1, 0.1]

# define initial conditions for ODE
x_max = 20001
x_int = 10000
c_init = permutations([1, x_int, x_max]) 

for i in list(c_init): 
    print (i) 

# %% find the steady state by integration
states = 0

x0 = np.random.random(1) * np.ones(n)
err_tole = 0.000001

t_final = 1000
t = np.linspace(0,t_final,200)
sol = odeint(f_ode, x0, t,args=(k,))

## check the integration solution
fig,ax = plt.subplots()
ax.plot(t,sol[:,0],label='x1')
ax.plot(t,sol[:,1],label='x2')
#ax.plot(t,result[:,2],label='R0=1')
ax.legend()
ax.set_xlabel('t')
ax.set_ylabel('x')

# tspan = np.linspace(0, 5, 100)
# yinit = [0, -3]
# c = [4, 3, -2, 0.5]
# # Solve differential equation
# sol = solve_ivp(lambda t, y: f(t, y, c), 

#                 [tspan[0], tspan[-1]], yinit, t_eval=tspan, rtol = 1e-5)
# # Plot states
# state_plotter(sol.t, sol.y, 1)


# double check if steady state is reache by considering the first 100 time points as BurnIn
#np.nonzero(t > 500))
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

# continue to check if multiple steady states exist

#%% calculate the eigenvalue of Jacobian matrix without and with diffusion matrix
# some codes from https://www.sympy.org/scipy-2017-codegen-tutorial/notebooks/20-ordinary-differential-equations.html were 
# very helpful
import sympy as sym
X = sym.symbols(('x0:2'))
K = sym.symbols(('k0:9'))
f_sym = sym.Matrix(f_ode(X, None, K))
J = f_sym.jacobian(X)
J_func = sym.lambdify((X, K),  J)

J_inputs = J.free_symbols

Xoverlap = np.zeros(len(X))
xx = []
for i in range(len(X)):
    if X[i] in J_inputs:
        Xoverlap[i] = 1
        xx.append(ss[i])

Koverlap = np.zeros(len(K))
kk = []
for i in range(len(K)):
    if K[i] in J_inputs:
        Koverlap[i] = 1
        kk.append(k[i])

#f = sy.Function('f')
#eq = sy.Eq(f(x).diff(x, 2) - 2*f(x).diff(x) + f(x), sy.sin(x)) sy.dsolve(eq)
#y = y0, y1 = sym.symbols('y0 y1')
#mu = sym.symbols('mu')
#J = sym.Matrix(f_ode(y, None, mu)).jacobian(y)
#J_func = sym.lambdify((y, t, mu), J)
#J
#inputs = [kk, xx]
S = J_func(ss, k)

w, v  =  np.linalg.eig(S)





