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
import os

#import numba
import scipy.integrate
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
#from RD_network_functions import state_plotter

import matplotlib.pyplot as plt

#import biocircuits
#import bokeh.io
#import bokeh.plotting
import panel as pn
pn.extension()

#import math
import itertools
import sympy as sym
from itertools import permutations
from skopt.space import Space
from skopt.sampler import Lhs

import sys, getopt
from RD_network_functions import *

#%% utility functions
def check_BurnIn_steadyState(sol, f_ode, k, S, n, x0, t_final):
    from scipy.integrate import odeint
    err_tole = 0.00001
    ss = np.zeros(n)
    
    ss_fluc = np.zeros(n)
    nb_passThreshold = 0
    
    for i in range(n):
        ss[i] = np.mean(sol[150:199, i])
        ss_fluc[i] = np.abs(np.mean(sol[150:199, i]) - np.mean(sol[100:149, i]))
        if ss_fluc[i] >= err_tole:
            nb_passThreshold = nb_passThreshold + 1
    nb_try = 0        
    while nb_passThreshold > 0 and nb_try < 5:
        nb_try = nb_try + 1
        t_final = t_final*2
        t = np.linspace(0, t_final, 200)
        sol = odeint(f_ode, x0, t,args=(k,S))
        ss = np.zeros(n)
        ss_fluc = np.zeros(n)
        nb_passThreshold = 0
        for i in range(n):
            ss[i] = np.mean(sol[150:199, i])
            ss_fluc[i] = np.abs(np.mean(sol[150:199, i]) - np.mean(sol[100:149, i]))
            if ss_fluc[i] >= err_tole:
                nb_passThreshold = nb_passThreshold + 1
    
    return  ss

def Multi_steadyStates(ss0, c_init, f_ode, k, S, n):
    import itertools
    from scipy.integrate import odeint
    
    ss_saved = [ss0]
    
    t_final = 2000
    t = np.linspace(0, t_final, 200)
    
    #x_init = np.logspace(-2, 4, nb_init)
    #c_init = itertools.combinations_with_replacement(x_init, n)
    for val in c_init:
        x0 = np.asarray(val)
        #print(x0)
        sol2 = odeint(f_ode, x0, t,args=(k,S))
        ss2 = check_BurnIn_steadyState(sol2, f_ode, k, S, n, x0, t_final)
        
        dists = np.ones(len(ss_saved))
        for j in range(len(ss_saved)):
            dists[j] = np.linalg.norm(ss2 -ss_saved[j])
        
        if (dists > 10**(-5)).all():
            ss_saved.append(ss2)
        #x0 = (*val)
        #print(x0)
    return ss_saved

#%% specify the network topology or model with parameters
def f_ode_simple(x, t, k): # simplified version of ode for specific S matrix
    #dRdt = np.empty(n)
    
    ## NT organoid phase III pattern selection:  Noggin, BMP, FoxA2
    dx0dt = k[0]  -     x[0] + k[5]*(1.0/2.0 * 1.0/(1.0 + (k[8]/x[1])**(2.0)) * 1.0/(1.0 + (k[9]/x[2])**(2.0))) # Noggin
    dx1dt = k[1] - k[3]*x[1] + k[6]*(1.0/(1.0+(k[10]/x[0])**(-2.0)) * 1.0/(2.0) * 1.0/(1.0 + (k[11]/x[2])**(2.0))) # BMP
    dx2dt = k[2] - k[4]*x[2] + k[7]*(1.0/2.0 * 1.0/(1.0 + (k[13]/x[1])**(-2.0)) * 1.0/(1.0 + (1.0/x[2])**(2.0))) # Foxa2
    
    dRdt = [dx0dt, dx1dt, dx2dt]
    
    return dRdt

# define ODE model in general form without diffusion 
def f_ode(x, t, k, S):
    #dRdt = np.empty(n)
    
    ## NT organoid phase III pattern selection:  Noggin, BMP, FoxA2
    dx0dt = k[0]  -     x[0] + k[5]*(1.0/(1.0+(1.0/x[0])**(2.0*S.iloc[0,0])) * 1.0/(1.0 + (k[8]/x[1])**(2.0*S.iloc[1, 0])) * 1.0/(1.0 + (k[9]/x[2])**(2.0*S.iloc[2, 0]))) # Noggin
    dx1dt = k[1] - k[3]*x[1] + k[6]*(1.0/(1.0+(k[10]/x[0])**(2.0*S.iloc[0,1])) * 1.0/(1.0 + (1.0/x[1])**(2.0*S.iloc[1, 1])) * 1.0/(1.0 + (k[11]/x[2])**(2.0*S.iloc[2, 1]))) # BMP
    dx2dt = k[2] - k[4]*x[2] + k[7]*(1.0/(1.0+(k[12]/x[0])**(2.0*S.iloc[0,2])) * 1.0/(1.0 + (k[13]/x[1])**(2.0*S.iloc[1, 2])) * 1.0/(1.0 + (1.0/x[2])**(2.0*S.iloc[2, 2]))) # Foxa2
    
    dRdt = [dx0dt, dx1dt, dx2dt]
    
    return dRdt


#%% linear stability analysis for given kinetic prameter;
# save the parameter, k, steady state, x, sampled diffusion d, wavenumber q and lambda (real and imaginary)
def linear_stability_test_param(n, f_ode, k, S, t_final, c_init, X, K, d_grid, q, i, outputDir):
    
    names = [sym.symbols('k0:' + str(len(k))), 
             sym.symbols('X0:' + str(n)),
             sym.symbols('noDiffusion0:2'),
             sym.symbols('d0:' + str(len(d_grid[0]))), 
             sym.symbols('q0:' + str(len(q))), 
             sym.symbols('lambda_re0:' + str(len(q))), 
             sym.symbols('lambda_im0:' + str(len(q))), 
             sym.symbols('eigenvec_posSign0:' + str(len(q)))
             ]
    
    names = [element for tupl in names for element in tupl]
    keep = pd.DataFrame(columns=names)
    
    # start with some test
    #k = np.asarray(k_grid[i])
    #k0 = k_grid[i]
    #k[0], k[1], k[2] = 0.1, 0.1, 0.1
    #k[3], k[4], k[5] = 0.3, 0.5, 0.4
    #k[6], k[7], k[8] = math.log(2)/10, math.log(2)/5, math.log(2)/10
    #k[6], k[7], k[8] = 30, 50, 20
    #k[9], k[10], k[11], k[12], k[13], k[14], k[15], k[16] = 14, 3, 0.2, 5, 10, 1, 2, 5
    
    # %% find the steady state by integration (initial guess)
    x0 = np.random.random(1) * np.ones(n)
    #x0 = [2.3, 0.4, 1.3]
    
    #f_ode_xk = f_ode(x, None, k, S)
    t = np.linspace(0, t_final, 200)
    sol = odeint(f_ode, x0, t, args=(k,S))
    #sol_test = odeint(f_ode_xk, x0, t, args=(k,))
    # check the integration solution with plot
    # fig,ax = plt.subplots()
    # ax.plot(t,sol[:,0],label='x1')
    # ax.plot(t,sol[:,1],label='x2')
    # ax.plot(t,sol[:,2],label='x3')
    # #ax.plot(t,result[:,2],label='R0=1')
    # ax.legend()
    # ax.set_xlabel('t')
    # ax.set_ylabel('x')
    # # sol[199, ]
    
    # double check if steady state is reache by considering the first 100 time points as BurnIn
    ss0 = check_BurnIn_steadyState(sol, f_ode, k, S, n, x0, t_final)
    
    #%% check if multiple steady states exist and save
    # define initial conditions for ODE
    ss_saved = Multi_steadyStates(ss0, c_init, f_ode, k, S, n)
    for kk in range(len(ss_saved)):
        ss = ss_saved[kk]    
        if any(ss <= 0) or any([isinstance(j, complex) for j in ss]):
            ss_saved.remove(ss)
    
    if len(ss_saved) > 0: 
        #%% calculate the eigenvalue of Jacobian matrix without and with diffusion matrix
        # some codes from https://www.sympy.org/scipy-2017-codegen-tutorial/notebooks/20-ordinary-differential-equations.html 
        # were very helpful
        f_sym = sym.Matrix(f_ode(X, None, K, S))
        J = f_sym.jacobian(X)
        J_func = sym.lambdify((X, K),  J)
        
        #%% eigenvalue computation with diffusion matrix to test if Turing instability 
        # disperion relation plot for specific set of parameters
        #d = [0, 0.01057, 1]
        
        # loop over stedy states
        for ii in range(len(ss_saved)):
            # ii = 0
            ss = ss_saved[ii]
            #J_inputs = J.free_symbols
            S1 = J_func(ss, k)
            
            try:
                w  =  np.linalg.eigvals(S1)
                #w, v = np.linalg.eig(S1)
                # max(w), S
            except:
                return True
            
            eigen0 = np.zeros(2)
            eigen0[0] = w.real.max()
            eigen0[1] = w.imag[np.argmax(w.real)]
            
            if eigen0[0] < 0.0: # if stabe without diffusion; otherwise don't continue the test
                
                for val in d_grid: # loop diffusion matrix for each steady state
                
                    d = np.asarray(val)
                    #d = [0.0, d[0], d[1]]
                    
                    lam_real = np.empty_like(q)
                    lam_im = np.empty_like(q)
                    eigenvec = []
                    # loop over the wavenumber 
                    for j in range(len(q)):
                        # j = 1
                        S2 = S1 - np.diag(np.multiply(d, q[j]**2))
                        #wk,vk =  np.linalg.eig(S2)
                        try:
                            wk, vk = np.linalg.eig(S2)
                        except:
                            return True
                        
                        lam_real[j] = wk.real.max()
                        lam_im[j] = wk.imag[np.argmax(wk.real)]
                        vec_sel = vk[:, np.argmax(wk.real)]
                        eigenvec.append(str(int(vec_sel[0] >0)) + ';' + str(int(vec_sel[1] >0)) + ';' + str(int(vec_sel[2] >0)))
                    #plt.plot(q, lam_real)
                    #plt.show()
                    #plt.axis([0, max(q), -1, 1])
                    #plt.axhline(y=0, color='r', linestyle='-'
                    #print(max(lam_real))
                    index_max = np.argmax(lam_real) 
                    lam_real_max = lam_real[index_max]
                    #lam_im_max = lam_im[index_max]
                    #q_max = q[index_max]
                
                    # save the result, k parameter, steady state, d parameters, q values, lambda_real, lambda imaginary
                    if lam_real_max >= 0:  
                        arr = np.concatenate((k, ss, eigen0,  d, q, lam_real, lam_im, eigenvec))
                        keep = keep.append(pd.DataFrame(arr.reshape(1,-1), columns=list(keep)), ignore_index=True)
    
    if keep.shape[0] > 1: # for one reaction kinetic parameter save one file
        keep.to_csv(outputDir + '/linear_stability_out_' + str(i) + '.csv', index = False) # Use Tab to seperate data

#%% linear stability analysis for single parameter                
def linear_stability_singleParam(i, k_grid_log, nb_params, Index_K_unsampled, n, f_ode, S, t_final, c_init, X, K, d_grid, q):
    # i = 50
    #if i % 100 == 0:
    #print(i)
        
    ks = np.asarray(k_grid_log[i])
    ks = np.power(10.0, ks) # transform to linear scale
        
    k = np.ones(nb_params)
    
    if len(ks) < len(k) :
        index_ks = 0
        for index_k in range(nb_params):
            if index_k not in Index_K_unsampled:
                k[index_k] = ks[index_ks]
                index_ks = index_ks + 1
    
        # test if the parameter assignment correct             
        for index_kns in Index_K_unsampled:
            if k[index_kns] > 1.0 or k[index_kns] < 1.0:
                print(" non smpled parameters assignment not correct  !")
                os._exit(1)
        
    linear_stability_test_param(n, f_ode, k, S, t_final, c_init, X, K, d_grid, q, i)

#%% numerical solution
# constant diffusion coeffients
def constant_diff_coeffs(c_tuple, t, x, diff_coeffs):
    n = len(c_tuple[0])
    return tuple([diff_coeffs[i] * np.ones(n) for i in range(len(c_tuple))])

# deinfe the reaction equations
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
    
def fode_3N2M_rxn(as_tuple, t, k):
    #dRdt = np.empty(n)
    
    a, s, f = as_tuple
    #x[0], x[1], x[2] = as_tuple
    ## the test example from Zheng et al, 2016, Fig.S1A
    #dx0dt = 55.14*x[2]**2/(18.48**2 + x[2]**2) + 0.1 - 1.341*x[0]
    #dx1dt = 29.28*x[2]**2/(15.90**2 + x[2]**2) + 0.1 - 0.3508*x[1]
    #dx2dt = 16.17*x[0]**2/(x[0]**2 + 0.6421**2)*1.316**2/(1.316**2 + x[1]**2) + 0.1 - 1.203*x[2]
    
    dx0dt = 55.14*f**2/(18.48**2 + f**2) + 0.1 - 1.341*a
    dx1dt = 29.28*f**2/(15.90**2 + f**2) + 0.1 - 0.3508*s
    dx2dt = 16.17*a**2/(a**2 + 0.6421**2)*1.316**2/(1.316**2 + s**2) + 0.1 - 1.203*f
    
        
    return (dx0dt, dx1dt, dx2dt)


def RD_numericalSolver():
    
    print('RD numerical solution 1D and 2D')
    
    ## here we start the example from http://be150.caltech.edu/2020/content/lessons/20_turing.html
    
    
    # Set up intial condition (using 500 grid points)
    a_0 = np.ones(500)
    s_0 = np.ones(500)
    f_0 = np.ones(500)
    #a_0
    # Make a small perturbation to a_0 by adding noise
    a_0 += 0.1 * np.random.rand(len(a_0))
    s_0 += 0.1 * np.random.rand(len(s_0))
    #f_0 += 0.1 * np.random.rand(len(f_0))
    
    # Time points
    t = np.linspace(0.0, 1000.0, 100)
    # Physical length of system
    L = 20
    
    # x-coordinates for plotting
    x = np.linspace(0, L, len(a_0))
    
    # Diffusion coefficients
    #diff_coeffs = (0.05, 1.0, 0)
    diff_coeffs = (0.07018, 1, 0.01057)
    # Reaction parameter (must be a tuple of params, even though only 1 for ASDM)
    k = np.ones(15)
    k[0], k[1], k[2] = 10, 10, 20 
    k[3], k[4], k[5] = 5, 5, 10
    #k[6], k[7], k[8] = math.log(2)/10, math.log(2)/5, math.log(2)/10
    k[6], k[7], k[8] = 0.1, 1, 0.1
    k[9], k[10], k[11], k[12], k[13], k[14] = 30, 50, 50, 70, 40, 80
    rxn_params = (k,)
    
    # periodic boundary condtion
    periodic_bc = False
    if periodic_bc:
        a_0[-1] = a_0[0]
        s_0[-1] = s_0[0]
        f_0[-1] = f_0[0]
    
    # Solve numeerically the RD with no-flux boundary condition
    conc = rd_solve((a_0, s_0, f_0),
        t,
        L=L,
        derivs_0=0,
        derivs_L=0,
        periodic_bc = periodic_bc,
        diff_coeff_fun=constant_diff_coeffs,
        diff_coeff_params=(diff_coeffs,),
        rxn_fun=fode_3N2M_rxn,
        rxn_params=rxn_params,
        mxstep=5000
    )
    
    #t_point = 1000000
    #i = np.searchsorted(t, t_point)
    i = len(t) - 1
    plt.plot(x, conc[0][i, :])
    plt.plot(x, conc[1][i, :], color="orange")
    plt.plot(x, conc[2][i, :], color="green")
    
        
    
#%% main function
def main(argv):
    
    modelfile = ''
    paramfile = ''
    try:
        opts, args = getopt.getopt(argv,"hm:p:",["model=", "param="])
    except getopt.GetoptError:
        print('python RD_numericalSolution_perturbation.py -m <modelfile> -p <paramfile>')
        sys.exit(2)
      
    for opt, arg in opts:
        if opt == '-h':
            print('python RD_numericalSolution_perturbation.py -m <modelfile> -p <paramfile>')
            sys.exit()
        elif opt in ("-m", "--modelfile"):
            modelfile = arg
        elif opt in ('-p', '--paramfile'):
            paramfile = arg
    
    
    #inputfile = '3N2M_topology_enumerate/Model_21.csv'
    #modelfile = '/Users/jiwang/workspace/imp/organoid_patterning/results/RD_topology_screening/topology_screening_3N2M_v2/topology_summary_selection_D.larger.1_lambda.neg.max.q_phase/table_params/Model_16.csv'
    #paramfile = '/Users/jiwang/workspace/imp/organoid_patterning/results/RD_topology_screening/topology_screening_3N2M_v2/topology_summary_selection_D.larger.1_lambda.neg.max.q_phase/table_params/params_saved_Model_16.csv'
    print('model file is -- ', modelfile)
    print('param file is -- ', paramfile)
    
    sys.exit()
    
    outputDir = os.path.basename()
    outputDir = './RD_numSolution_perturbation/' + outputDir.rsplit('.', 1)[0]
    print('Output directory is ', outputDir)
    
    try:
        os.makedirs(outputDir)
    except FileExistsError:
        # directory already exists
        pass
    
    # total number for parameter sampling 
    nb_sampling_parameters = 1000 # reaction parameters
    
    nb_sampling_diffusion = 50 # diffusion rate
    nb_sampling_init = 3 # nb of initial condtion sampled
    
    q = 2*3.14159 / np.logspace(-2, 3.0, num=nb_sampling_diffusion) # wavenumber
    
    n = 3 # nb of node
    nb_params = 14
    #binary_diffusor = [1, 1, 0]
    
    print('--  main function starts --')
    
    import time
    start_time = time.process_time()
    
    # read the network topology
    S = pd.read_csv(inputfile, index_col=0) 
    
    if S.shape[0] != 3 or S.shape[1] != 3:
        print("Required 3x3 matrix for network topology !")
        os._exit(1)
    
    # define symbolic variables for Jacobian matrix 
    X = sym.symbols(('x0:' + str(n)))
    K = sym.symbols(('k0:' + str(nb_params)))
    
    ## keep a record of unsampled parameters
    #K_total = [None] * nb_params
    #for index_par in range(nb_params):
    #    K_total[index_par] = 'k' + str(index_par)
        
    Index_K_unsampled = []
    k_length = 8 # nb of reaction parameters: 3* number of nodes (3*3) + number of interactions (6)
    for index_j in range(3): 
        for index_i in range(3):
            #print(S.iloc[index_j, index_i])
            if index_i != index_j:
                if np.abs(S.iloc[index_j, index_i]) > 0: 
                    k_length = k_length + 1
                else:
                    if index_i == 0 and index_j == 1:
                        Index_K_unsampled.append(8)
                    elif index_i ==0 and index_j == 2:
                        Index_K_unsampled.append(9)
                    elif index_i == 1 and index_j == 0:
                        Index_K_unsampled.append(10)
                    elif index_i == 1 and index_j == 2:
                        Index_K_unsampled.append(11)
                    elif index_i == 2 and index_j == 0:
                        Index_K_unsampled.append(12)
                    elif index_i == 2 and index_j == 1:
                        Index_K_unsampled.append(13)
    if len(Index_K_unsampled) + k_length != nb_params:
        print(" nb of sampled parameters not correct  !")
        os._exit(1)
                            
    #%% sampling the parameters, which node is difusor and diffusion coeffs
    
    ## lhs sampling for parameter
    np.random.seed(123)
    
    list_dimensions = [(-2, 2)]*k_length
    space = Space(list_dimensions)
    
    lhs = Lhs(criterion="maximin", iterations=1000)
    k_grid_log = lhs.generate(space.dimensions, nb_sampling_parameters)
    
    # diffusion rate sampling
    d_range = np.logspace(-3, 3.0, num = nb_sampling_diffusion)
    d_grid = list(itertools.product(np.ones(1),  d_range, np.zeros(1)))
    
    # initial conditions: each node has 3 initial values
    x_init = np.logspace(-1, 4, nb_sampling_init)
    c_init = itertools.combinations_with_replacement(x_init, n)
    
    # time 
    t_final = 1000
    
    print(time.process_time() - start_time, "seconds to set up parameters ")
    
    #%% big loop over each k parameter vector and save the result for each sampled d combination
    start_time = time.process_time()
    
    # try to parallize the for loop
    #from joblib import Parallel, delayed
    #import multiprocessing
    # #pool_obj = multiprocessing.Pool()
    # pool = multiprocessing.Pool()
    # args = ((foo, bar, foobar, baz) 
    #     for foo in range(3) 
    #     for bar in range(5) 
    #     for baz in range(4) 
    #     for foobar in range(10))
    # pool.starmap(calculation, args)
    # pool.close()
    # pool.join()
    # pool_obj.map(sumall, range(0, len(k_grid_log)))
    
    #for i in range(len(k_grid)):
    for i in range(len(k_grid_log)):
        
        if i % 100 == 0 and i > 0:
            print(i)
        ks = np.asarray(k_grid_log[i])
        ks = np.power(10.0, ks) # transform to linear scale
        
        k = np.ones(nb_params)
        
        if len(ks) < len(k) :
            index_ks = 0
            for index_k in range(nb_params):
                if index_k not in Index_K_unsampled:
                    k[index_k] = ks[index_ks]
                    index_ks = index_ks + 1
    
            # test if the parameter assignment correct             
            for index_kns in Index_K_unsampled:
                if k[index_kns] > 1.0 or k[index_kns] < 1.0:
                    print(" non smpled parameters assignment not correct  !")
                    os._exit(1)
        
        linear_stability_test_param(n, f_ode, k, S, t_final, c_init, X, K, d_grid, q, i, outputDir)
        #Parallel(n_jobs=2)(delayed(linear_stability_singleParam)(i, k_grid_log, nb_params, Index_K_unsampled, n, f_ode, S, t_final, c_init, X, K, d_grid, q) for i in range(len(k_grid_log)))
    
    print(time.process_time() - start_time, "seconds for for loop")
    
    
if __name__ == "__main__":

    main(sys.argv[1:])


