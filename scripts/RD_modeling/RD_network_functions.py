#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 29 11:44:52 2021

@author: jingkui.wang
"""
import numpy as np
import pandas as pd
import numba
import scipy.integrate
import biocircuits
import bokeh.io
import bokeh.plotting

import panel as pn
import matplotlib.pyplot as plt

def dc_dt(
    c,
    t,
    x,
    derivs_0,
    derivs_L,
    diff_coeff_fun,
    diff_coeff_params,
    rxn_fun,
    rxn_params,
    n_species,
    h,
):
    """
    Time derivative of concentrations in an R-D system
    for constant flux BCs.

    Parameters
    ----------
    c : ndarray, shape (n_species * n_gridpoints)
        The concentration of the chemical species interleaved in a
        a NumPy array.  The interleaving allows us to take advantage
        of the banded structure of the Jacobian when using the
        Hindmarsh algorithm for integrating in time.
    t : float
        Time.
    derivs_0 : ndarray, shape (n_species)
        derivs_0[i] is the value of the diffusive flux,
        D dc_i/dx, at x = 0, the leftmost boundary of the domain of x.
    derivs_L : ndarray, shape (n_species)
        derivs_0[i] is the value of the diffusive flux,
        D dc_i/dx, at x = L, the rightmost boundary of the domain of x.
    diff_coeff_fun : function
        Function of the form diff_coeff_fun(c_tuple, t, x, *diff_coeff_params).
        Returns an tuple where entry i is a NumPy array containing
        the diffusion coefficient of species i at the grid points.
        c_tuple[i] is a NumPy array containing the concentrations of
        species i at the grid poitns.
    diff_coeff_params : arbitrary
        Tuple of parameters to be passed into diff_coeff_fun.
    rxn_fun : function
        Function of the form rxn_fun(c_tuple, t, *rxn_params).
        Returns an tuple where entry i is a NumPy array containing
        the net rate of production of species i by chemical reaction
        at the grid points.  c_tuple[i] is a NumPy array containing
        the concentrations of species i at the grid poitns.
    rxn_params : arbitrary
        Tuple of parameters to be passed into rxn_fun.
    n_species : int
        Number of chemical species.
    h : float
        Grid spacing (assumed to be constant)

    Returns
    -------
    dc_dt : ndarray, shape (n_species * n_gridpoints)
        The time derivatives of the concentrations of the chemical
        species at the grid points interleaved in a NumPy array.
    """
    # Tuple of concentrations
    c_tuple = tuple([c[i::n_species] for i in range(n_species)])

    # Compute diffusion coefficients
    D_tuple = diff_coeff_fun(c_tuple, t, x, *diff_coeff_params)

    # Compute reaction terms
    rxn_tuple = rxn_fun(c_tuple, t, *rxn_params)

    # Return array
    conc_deriv = np.empty_like(c)

    # Convenient array for storing concentrations
    da_dt = np.empty(len(c_tuple[0]))

    # Useful to have square of grid spacing around
    h2 = h ** 2

    # Compute diffusion terms (central differencing w/ Neumann BCs)
    for i in range(n_species):
        # View of concentrations and diffusion coeff. for convenience
        a = np.copy(c_tuple[i])
        D = np.copy(D_tuple[i])

        # Time derivative at left boundary
        da_dt[0] = D[0] / h2 * 2 * (a[1] - a[0] - h * derivs_0[i])

        # First derivatives of D and a
        dD_dx = (D[2:] - D[:-2]) / (2 * h)
        da_dx = (a[2:] - a[:-2]) / (2 * h)

        # Time derivative for middle grid points
        da_dt[1:-1] = D[1:-1] * np.diff(a, 2) / h2 + dD_dx * da_dx

        # Time derivative at left boundary
        da_dt[-1] = D[-1] / h2 * 2 * (a[-2] - a[-1] + h * derivs_L[i])

        # Store in output array with reaction terms
        conc_deriv[i::n_species] = da_dt + rxn_tuple[i]

    return conc_deriv


def rd_solve(
    c_0_tuple,
    t,
    L=1,
    derivs_0=0,
    derivs_L=0,
    diff_coeff_fun=None,
    diff_coeff_params=(),
    rxn_fun=None,
    rxn_params=(),
    rtol=1.49012e-8,
    atol=1.49012e-8,
):
    """
    Parameters
    ----------
    c_0_tuple : tuple
        c_0_tuple[i] is a NumPy array of length n_gridpoints with the
        initial concentrations of chemical species i at the grid points.
    t : ndarray
        An array of time points for which the solution is desired.
    L : float
        Total length of the x-domain.
    derivs_0 : ndarray, shape (n_species)
        derivs_0[i] is the value of dc_i/dx at x = 0.
    derivs_L : ndarray, shape (n_species)
        derivs_L[i] is the value of dc_i/dx at x = L, the rightmost
        boundary of the domain of x.
    diff_coeff_fun : function
        Function of the form diff_coeff_fun(c_tuple, x, t, *diff_coeff_params).
        Returns an tuple where entry i is a NumPy array containing
        the diffusion coefficient of species i at the grid points.
        c_tuple[i] is a NumPy array containing the concentrations of
        species i at the grid poitns.
    diff_coeff_params : arbitrary
        Tuple of parameters to be passed into diff_coeff_fun.
    rxn_fun : function
        Function of the form rxn_fun(c_tuple, t, *rxn_params).
        Returns an tuple where entry i is a NumPy array containing
        the net rate of production of species i by chemical reaction
        at the grid points.  c_tuple[i] is a NumPy array containing
        the concentrations of species i at the grid poitns.
    rxn_params : arbitrary
        Tuple of parameters to be passed into rxn_fun.
    rtol : float
        Relative tolerance for solver.  Default os odeint's default.
    atol : float
        Absolute tolerance for solver.  Default os odeint's default.

    Returns
    -------
    c_tuple : tuple
        c_tuple[i] is a NumPy array of shape (len(t), n_gridpoints)
        with the initial concentrations of chemical species i at
        the grid points over time.

    Notes
    -----
    .. When intergrating for long times near a steady state, you
       may need to lower the absolute tolerance (atol) because the
       solution does not change much over time and it may be difficult
       for the solver to maintain tight tolerances.
    """
    # Number of grid points
    n_gridpoints = len(c_0_tuple[0])

    # Number of chemical species
    n_species = len(c_0_tuple)

    # Grid spacing
    h = L / (n_gridpoints - 1)

    # Grid points
    x = np.linspace(0, L, n_gridpoints)

    # Set up boundary conditions
    if np.isscalar(derivs_0):
        derivs_0 = np.array(n_species * [derivs_0])
    if np.isscalar(derivs_L):
        derivs_L = np.array(n_species * [derivs_L])

    # Set up parameters to be passed in to dc_dt
    params = (
        x,
        derivs_0,
        derivs_L,
        diff_coeff_fun,
        diff_coeff_params,
        rxn_fun,
        rxn_params,
        n_species,
        h,
    )

    # Set up initial condition
    c0 = np.empty(n_species * n_gridpoints)
    for i in range(n_species):
        c0[i::n_species] = c_0_tuple[i]

    # Solve using odeint, taking advantage of banded structure
    c = scipy.integrate.odeint(
        dc_dt,
        c0,
        t,
        args=params,
        ml=n_species,
        mu=n_species,
        rtol=rtol,
        atol=atol,
    )

    return tuple([c[:, i::n_species] for i in range(n_species)])


def state_plotter(times, states, fig_num):
    num_states = np.shape(states)[0]
    num_cols = int(np.ceil(np.sqrt(num_states)))
    num_rows = int(np.ceil(num_states / num_cols))
    plt.figure(fig_num)
    plt.clf()
    fig, ax = plt.subplots(num_rows, num_cols, num=fig_num, clear=True,
                         squeeze=False)
    for n in range(num_states):
        row = n // num_cols
        col = n % num_cols
        ax[row][col].plot(times, states[n], 'k.:')
        ax[row][col].set(xlabel='Time',
                         ylabel='$y_{:0.0f}(t)$'.format(n),
                         title='$y_{:0.0f}(t)$ vs. Time'.format(n))
        
    for n in range(num_states, num_rows * num_cols):
        fig.delaxes(ax[n // num_cols][n % num_cols])

    fig.tight_layout()

    return fig, ax

def check_BurnIn_steadyState(sol, f_ode, k, n, x0, t_final):
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
        sol = odeint(f_ode, x0, t,args=(k,))
        ss = np.zeros(n)
        ss_fluc = np.zeros(n)
        nb_passThreshold = 0
        for i in range(n):
            ss[i] = np.mean(sol[150:199, i])
            ss_fluc[i] = np.abs(np.mean(sol[150:199, i]) - np.mean(sol[100:149, i]))
            if ss_fluc[i] >= err_tole:
                nb_passThreshold = nb_passThreshold + 1
    
    return  ss