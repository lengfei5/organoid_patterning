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
turing = 0 # 0 

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

# continue to check if multiple steady states exist (to do)

#%% calculate the eigenvalue of Jacobian matrix without and with diffusion matrix
# some codes from https://www.sympy.org/scipy-2017-codegen-tutorial/notebooks/20-ordinary-differential-equations.html were 
# very helpful
import sympy as sym
X = sym.symbols(('x0:2'))
K = sym.symbols(('k0:9'))
f_sym = sym.Matrix(f_ode(X, None, K))
J = f_sym.jacobian(X)
J_func = sym.lambdify((X, K),  J)

#J_inputs = J.free_symbols
S = J_func(ss, k)
w, v  =  np.linalg.eig(S)

# Xoverlap = np.zeros(len(X))
# xx = []
# for i in range(len(X)):
#     if X[i] in J_inputs:
#         Xoverlap[i] = 1
#         xx.append(ss[i])

# Koverlap = np.zeros(len(K))
# kk = []
# for i in range(len(K)):
#     if K[i] in J_inputs:
#         Koverlap[i] = 1
#         kk.append(k[i])

#f = sy.Function('f')
#eq = sy.Eq(f(x).diff(x, 2) - 2*f(x).diff(x) + f(x), sy.sin(x)) sy.dsolve(eq)
#y = y0, y1 = sym.symbols('y0 y1')
#mu = sym.symbols('mu')
#J = sym.Matrix(f_ode(y, None, mu)).jacobian(y)
#J_func = sym.lambdify((y, t, mu), J)
#J
#inputs = [kk, xx]

#%% eigenvalue computation with diffusion matrix to test if Turing instability 
# disperion relation plot for specific set of parameters
diffusing_nodes = binary_diffusor

k = np.linspace(0, 100, 500)
lam_real = np.empty_like(k)
lam_im = np.empty_like(k)

for j in range(len(k)):
    #j = 1
    S2 = S - np.diag(np.multiply(d, k[j]*k[j]))
    #wk,vk =  np.linalg.eig(S2)
    wk = np.linalg.eigvals(S2)
    lam_real[j] = wk.real.max()
    lam_im[j] = wk.imag[np.argmax(wk.real)]
    
plt.plot(k, lam_real)
plt.axis([0, 8, -5, 5])
plt.axhline(y=0, color='r', linestyle='-')
plt.show()

index_max = np.argmax(lam_real) 
lam_real_max = lam_real[index_max]
lam_im_max = lam_im[index_max]
k_max = k[index_max]


# define turing pattern types according to the eigenvalues (to finish)
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

# orignial code from http://be150.caltech.edu/2020/content/lessons/20_turing.html
# try to make panal sidebar, not used for the moment

# d_slider = pn.widgets.FloatSlider(
#     name="d", start=0.01, end=1, value=0.05, step=0.01, width=150
# )

# mu_slider = pn.widgets.FloatSlider(
#     name="μ", start=0.01, end=2, value=1.5, step=0.005, width=150
# )

# @pn.depends(d_slider.param.value, mu_slider.param.value)
# def plot_dispersion_relation(d, mu):
#     d = np.linspace(0.01, 1, 150)
#     mu = np.linspace(0.01, 2, 150)
    
#     d = 0.05
#     mu = 1
#     k = np.linspace(0, 10, 200)
#     lam_max_real_part = dispersion_relation(k, d, mu)
#     plt.plot(k, lam_max_real_part)
#     plt.axis([0, 10, -5, 5])
#     plt.axhline(y=0, color='r', linestyle='-')
#     plt.show()

#     p = bokeh.plotting.figure(
#         frame_width=350,
#         frame_height=200,
#         x_axis_label="k",
#         y_axis_label="Re[λ-max]",
#         x_range=[0, 10],
#     )
#     p.line(k, lam_max_real_part, color="black", line_width=2)

#     return p


# pn.Column(
#     pn.Row(d_slider, mu_slider), pn.Spacer(height=20), plot_dispersion_relation
# )


#%% numerical solution of RD equation 
from RD_network_functions import *

def constant_diff_coeffs(c_tuple, t, x, diff_coeffs):
    n = len(c_tuple[0])
    return tuple([diff_coeffs[i] * np.ones(n) for i in range(len(c_tuple))])


def asdm_rxn(as_tuple, t, mu):
    """
    Reaction expression for activator-substrate depletion model.

    Returns the rate of production of activator and substrate, respectively.

    r_a = a**2 * s - a
    r_s = mu * (1 - a**2 * s)
    """
    # Unpack concentrations
    a, s = as_tuple

    # Compute and return reaction rates
    a2s = a ** 2 * s
    return (a2s - a, mu * (1.0 - a2s))


# Set up steady state (using 500 grid points)
a_0 = np.ones(500)
s_0 = np.ones(500)

# Make a small perturbation to a_0
a_0 += 0.001 * np.random.rand(len(a_0))

# Time points
t = np.linspace(0.0, 1000.0, 100)

# Physical length of system
L = 20.0

# x-coordinates for plotting
x = np.linspace(0, L, len(a_0))

# Diffusion coefficients
diff_coeffs = (0.05, 1.0)

# Reaction parameter (must be a tuple of params, even though only 1 for ASDM)
rxn_params = (1.5,)

# Solve
conc = rd_solve(
    (a_0, s_0),
    t,
    L=L,
    derivs_0=0,
    derivs_L=0,
    diff_coeff_fun=constant_diff_coeffs,
    diff_coeff_params=(diff_coeffs,),
    rxn_fun=asdm_rxn,
    rxn_params=rxn_params,
)

t_point = 1000000
i = np.searchsorted(t, t_point)
i = 98

plt.plot(x, conc[0][i, :])
plt.plot(x, conc[1][i, :],
        color="tomato"
    ) 

# t_slider = pn.widgets.FloatSlider(
#     name="t", start=t[0], end=t[-1], value=t[-1], step=t[1] - t[0]
# )

@pn.depends(t_slider)
def plot_turing(t_point):
    i = np.searchsorted(t, t_point)

    p = bokeh.plotting.figure(
        frame_width=400,
        frame_height=200,
        x_axis_label="x",
        y_axis_label="a, s",
        x_range=[0, L],
        y_range=[0, np.concatenate(conc).max()*1.02],
    )

    p.line(x, conc[0][i, :], legend_label="activator", line_width=2)
    p.line(
        x,
        conc[1][i, :],
        color="tomato",
        legend_label="substrate",
        line_width=2,
    )

    return p

# pn.Column(t_slider, plot_turing)

# xy plot
# bokeh.io.show(
#     biocircuits.xyt_plot(
#         x,
#         conc,
#         t,
#         legend_names=["activator", "substrate"],
#         palette=["#1f77b4", "tomato"],
#     )
# )

import rdsolver
from bokeh.io import output_notebook
from bokeh.plotting import figure, show
#notebook_url = 'localhost:8888'

# Load a standard ASDM model
D, beta, gamma, f, f_args, homo_ss = rdsolver.models.asdm()

# Set up the space and time grid
n = (32, 32)
L = (50, 50)
t = np.linspace(0, 100000, 100)

# Initial condition and solve
c0 = rdsolver.initial_condition(uniform_conc=homo_ss, n=n, L=L)
c = rdsolver.solve(c0, t, D=D, beta=beta, gamma=gamma, f=f, f_args=f_args, L=L)

# Interpolate the solution
c_interp = rdsolver.viz.interpolate_concs(c)

#plt.plot(c_interp[0, :,:,-1], c_interp[1, :,:,-1])
plt.imshow(c_interp[:,:,:, 99])
plt.show()

# Display the final time point
bokeh.io.show(rdsolver.viz.display_single_frame(c_interp, i = -1, frame_height=300), 
              browser=None, notebook_handle=False)

bokeh.io.show(
    rdsolver.viz.display_notebook(t, c_interp, frame_height=300),
    notebook_url=notebook_url,
)



