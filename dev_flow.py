#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: jd
contact: jamesduv@umich.edu
affiliation: University of Michigan, Department of Aerospace Eng., CASLAB
"""


import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt

from dev_thermo import get_gas_viscosity

def comp_flow(design_vars, flow_bc, therm_out, geom_out, mesh_out, mesh_param):
    '''Flow model where transport properties are treated as constant along engine axis'''

    z_inner     = mesh_out['zones']['inner']['zc']
    r_inner     = mesh_out['zones']['inner']['rc']
    Ac_inner    = mesh_out['zones']['inner']['Ac']
    dz_vec      = mesh_out['zones']['global']['dz_vec']
    rc          = design_vars['rc']

    rco         = mesh_out['zones']['global']['rco']
    rci         = mesh_out['zones']['global']['rci']
    beta_vec    = mesh_out['zones']['global']['beta_vec']
    Ac_channel  = mesh_out['zones']['global']['Ac_channel']
    nz_glob     = Ac_channel.shape[0]
    rho_f       = flow_bc['rho_fuel']
    mu_f        = flow_bc['mu_fuel']

    nz_c    = mesh_param['nz_c']
    nz_nc   = mesh_param['nz_nc']
    nz_nd   = mesh_param['nz_nd']
    row_start_nc    = mesh_out['zones']['global']['row_start_nc']
    row_t           = mesh_out['zones']['global']['row_t']
    row_start_nd    = mesh_out['zones']['global']['row_start_nd']
    beta_nc     = design_vars['beta_nc']
    beta_nd     = design_vars['beta_nd']

    Ac_c        = np.pi * rc**2
    Tc          = therm_out['Tc']
    pc          = therm_out['pc']
    rhoc        = therm_out['rhoc']
    R           = therm_out['R_sp']
    gamma       = therm_out['gamma']
    At          = geom_out['At']
    mu          = therm_out['mu']
    mdot        = therm_out['mdot_tot']
    mdot_f      = design_vars['mdot_f']
    mdot_cool   = design_vars['f_cool'] * mdot_f

    flow_sol = {}

    ### combustion chamber - constant flow properties
    uc          = mdot / (rhoc * Ac_c)
    ac          = np.sqrt(gamma * R * Tc)
    Mc          = uc / ac
    nz_c        = mesh_param['nz_c']
    Re_c        = 2 * rhoc * uc * rc / mu

    Tc_vec      = np.repeat([Tc], repeats = nz_c)
    pc_vec      = np.repeat([pc], repeats = nz_c)
    rhoc_vec    = np.repeat([rhoc], repeats = nz_c)
    uc_vec      = np.repeat([uc], repeats = nz_c)
    Mc_vec      = np.repeat([Mc], repeats = nz_c)
    ac_vec      = np.repeat([ac], repeats = nz_c)
    Rec_vec     = np.repeat([Re_c], repeats = nz_c)
    # gamma_c_vec     = np.repeat([gamma_c], repeats = nz_c)

    flow_sol['c'] = {
        'T'     : Tc_vec,
        'p'     : pc_vec,
        'rho'   : rhoc_vec,
        'u'     : uc_vec,
        'M'     : Mc_vec,
        'a'     : ac_vec,
        'Re'    : Rec_vec
        }

    ### compute stagnation conditions from chamber properties
    common  = (0.5 * (gamma-1))*Mc**2
    T0      = Tc * (1 + common)
    p0      = pc * (1 + common)**(gamma / (gamma-1))
    rho0    = p0 / (R * T0)

    flow_sol['stagnation'] = {
        'T'     : T0,
        'p'     : p0,
        'rho'   : rho0
        }

    ### converging nozzle section
    # set first value to chamber properties, due to repeated node
    nz_nc       = mesh_param['nz_nc']
    Tnc         = np.zeros(nz_nc)
    Tnc[0]      = Tc
    pnc         = np.zeros(nz_nc)
    pnc[0]      = pc
    rhonc       = np.zeros(nz_nc)
    rhonc[0]    = rhoc
    unc         = np.zeros(nz_nc)
    unc[0]      = uc
    Mnc         = np.zeros(nz_nc)
    Mnc[0]      = Mc
    anc         = np.zeros(nz_nc)
    anc[0]      = ac
    Renc        = np.zeros(nz_nc)
    Renc[0]     = Re_c

    # fma         = mach_area
    fma         = mach_area_sq
    fmap        = mach_area_sq_fp
    offset      = nz_c
    sol_nc      = []
    xtol        = 1e-12
    for ii in np.arange(1, nz_nc):
        Azcur       = Ac_inner[offset + ii]
        args        = (Azcur, At, gamma)

        # solution without providing gradient
        #bracket     = (0.05, 1)
        # sol_cur     = opt.root_scalar(fma, args = args, bracket=bracket, xtol = xtol)

        x0          = Mnc[ii-1]
        sol_cur     = opt.root_scalar(fma, args = args, x0 = x0,
                                      xtol = xtol, fprime = fmap, method='newton' )
        # sol_cur     = opt.root_scalar(fma, args = args, xtol = xtol,
        #                               fprime=fmap, x0=x0)
        if sol_cur.root <= 0:
            Mnc[ii] = Mnc[ii-1]
        else:
            Mnc[ii]     = sol_cur.root

        sol_nc.append(sol_cur)

        Tnc[ii]     = T_isen(T0, Mnc[ii], gamma)
        pnc[ii]     = p_isen(p0, Mnc[ii], gamma)
        rhonc[ii]   = rho_isen(rho0, Mnc[ii], gamma)
        anc[ii]     = np.sqrt(gamma * R * Tnc[ii])
        unc[ii]     = Mnc[ii] * anc[ii]
        # unc[ii]     = mdot / (rhonc[ii] * Azcur)
        # anc[ii]     = unc[ii] / Mnc[ii]
        # Tnc[ii]     = (anc[ii]**2) / (gamma * R)
        Renc[ii]    = (2 * rhonc[ii] * unc[ii] * r_inner[offset + ii]) / mu

    flow_sol['nc'] = {
        'T'     : Tnc,
        'p'     : pnc,
        'rho'   : rhonc,
        'u'     : unc,
        'M'     : Mnc,
        'a'     : anc,
        'Re'    : Renc,
        'sol_ls'    : sol_nc
        }

    flow_sol['t'] = {
        'T'     : Tnc[-1],
        'p'     : pnc[-1],
        'rho'   : rhonc[-1],
        'u'     : unc[-1],
        'M'     : Mnc[-1],
        'a'     : anc[-1],
        'Re'    : Renc[-1]
        }

    ### diverging nozzle section
    nz_nd       = mesh_param['nz_nd']
    Tnd         = np.zeros(nz_nd)
    pnd         = np.zeros(nz_nd)
    rhond       = np.zeros(nz_nd)
    und         = np.zeros(nz_nd)
    Mnd         = np.zeros(nz_nd)
    andd        = np.zeros(nz_nd)
    Rend        = np.zeros(nz_nd)

    fma         = mach_area
    offset      = nz_c + nz_nc
    sol_nd      = []
    xtol        = 1e-12

    for ii in np.arange(nz_nd):
        Azcur       = Ac_inner[offset + ii]
        args        = (Azcur, At, gamma)
        x0          = 2
        # x1          = 3
        # sol_cur     = opt.root_scalar(fma, args = args, x0=x0, x1=x1, xtol = xtol)
        sol_cur     = opt.root_scalar(fma, args = args, xtol = xtol,
                                      fprime=fmap, x0=x0, method='newton')

        # sol_cur     = opt.root_scalar(fma, args = args, bracket=bracket,
                                      # xtol = xtol, fprime = fmap )
        sol_nd.append(sol_cur)

        Mnd[ii]     = sol_cur.root
        pnd[ii]     = p_isen(p0, Mnd[ii], gamma)
        rhond[ii]   = rho_isen(rho0, Mnd[ii], gamma)
        Tnd[ii]     = T_isen(T0, Mnd[ii], gamma)
        andd[ii]    = np.sqrt(gamma * R * Tnd[ii])
        und[ii]     = Mnd[ii] * andd[ii]

        # und[ii]     = mdot / (rhond[ii] * Azcur)
        # andd[ii]    = und[ii] / Mnd[ii]
        # Tnd[ii]     = (andd[ii]**2) / (gamma * R)
        Rend[ii]    = (2 * rhond[ii] * und[ii] * r_inner[offset + ii]) / mu

    flow_sol['nd'] = {
        'T'     : Tnd,
        'p'     : pnd,
        'rho'   : rhond,
        'u'     : und,
        'M'     : Mnd,
        'a'     : andd,
        'Re'    : Rend,
        'sol_ls'    : sol_nd
        }

    #concatenate global vectors
    T_glob      = np.concatenate((Tc_vec, Tnc[1:], Tnd[1:]), axis=0)
    p_glob      = np.concatenate((pc_vec, pnc[1:], pnd[1:]), axis=0)
    rho_glob    = np.concatenate((rhoc_vec, rhonc[1:], rhond[1:]), axis=0)
    u_glob      = np.concatenate((uc_vec, unc[1:], und[1:]), axis=0)
    M_glob      = np.concatenate((Mc_vec, Mnc[1:], Mnd[1:]), axis=0)
    a_glob      = np.concatenate((ac_vec, anc[1:], andd[1:]), axis=0)
    Re_glob     = np.concatenate((Rec_vec, Renc[1:], Rend[1:]), axis=0)

    flow_sol['global'] = {
        'T'     : T_glob,
        'p'     : p_glob,
        'rho'   : rho_glob,
        'u'     : u_glob,
        'M'     : M_glob,
        'a'     : a_glob,
        'Re'    : Re_glob,
        }

    ### $
    u_ch    = mdot_cool / (rho_f * Ac_channel)
    hd      = 2 * (rco - rci)   #methane hydraulic diameter: length scale for  Re
    Re_ch   = (rho_f * u_ch * hd) / mu_f

    # compute pressure vs z in cooling channel
    ff_vec  = np.zeros(nz_glob)
    dp_ch   = np.zeros(nz_glob)

    f_loss = 0.03
    for ii in np.arange(nz_glob):
        # ff_vec[ii] = friction_factor(Re_ch[ii])
        L       = dz_vec[ii] / np.cos(beta_vec[ii] * np.pi / 180)
        dp_ch[ii] = 0.5 * f_loss * rho_f * u_ch[ii]**2 * (L/hd[ii])
    dp_ch_tot = dp_ch.sum()

    # set cooling channel outlet pressure equal to
    #1.5x chamber pressure
    pch_out = 1.5 * pc

    pch_vec = np.zeros(nz_glob)
    for ii in np.arange(nz_glob):
        if ii == 0:
            pch_vec[ii] = pch_out + dp_ch[ii]
        else:
            pch_vec[ii] = pch_vec[ii-1] + dp_ch[ii]

    flow_sol['channel'] = {
        'u'             : u_ch,
        'hydraulic_d'   : hd,
        'Re'            : Re_ch,
        'dp_vec'        : dp_ch,
        'dp_tot'        : dp_ch_tot,
        'p'             : pch_vec
        }
    return flow_sol

def mach_area(x, Az, At, gamma):
    '''Mach-number area relation relative to throat conditions, written 
    in residual form'''

    Arat = Az / At

    t1      = 1 / x
    common  = ((gamma-1)/2)
    num     = 1 + (common * x**2)
    tfrac   = ( num / (1 + common))**((gamma+1) / (gamma-1))
    t2      = np.sqrt(tfrac)
    # t2      = tfrac
    f       = Arat - (t1 * t2)
    return f

def mach_area_sq(x, Az, At, gamma):
    Arat = (Az / At)**2

    t1      = 1/ (x**2)
    common  = ((gamma-1)/2)
    num     = 1 + (common * x**2)
    t2   = ( num / (1 + common))**((gamma+1) / (gamma-1))
    # t2      = np.sqrt(tfrac)
    # t2      = tfrac
    f       = Arat - (t1 * t2)
    return f

def mach_area_sq_fp(x, Az, At, gamma):
    r = (gamma+1)/(gamma-1)
    common  = (0.5 * (gamma-1))
    T1      = (1 + (common * x**2)) / (1 + common)

    a1  = (2/(x**3)) * T1**r
    a2  = (r * T1**(r-1) / (x**2)) * ((2 * x * common)/(1+common))
    fp = a1 -a2
    return fp


def p_isen(p0, M, gamma):
    '''Compute the pressure given stagnation pressure and local Mach number
    using isentropic relation'''

    t1  = 1 + (0.5 * (gamma-1) * M**2)
    p   = p0 * (t1**(-gamma / (gamma - 1)))
    return p

def rho_isen(rho0, M, gamma):
    '''Compute the desnity given the stagnation density and local Mach number
    using isentropic relation'''

    t1  = 1 + (0.5 * (gamma-1) * M**2)
    rho     = rho0 * (t1**(-1 / (gamma-1)))
    return rho

def T_isen(T0, M, gamma):
    '''Compute the temperature given stagnation temperature and local Mach number
    using isentropic relation'''

    t1  = 1 + (0.5 * (gamma-1) * M**2)
    T   = T0 / t1
    return T

def friction_factor(Re):
    f = (0.79 * np.log(Re) - 1.64)**(-2)
    return f


def test_ma():
    At      = 0.123
    Az      = 0.27
    gamma   = 1.199
    nx      = 100
    xvec    = np.linspace(0.05, 1, nx)
    res     = np.zeros(nx)
    for ii in np.arange(nx):
        res[ii] = mach_area(xvec[ii], Az, At, gamma)

    fig = plt.figure()
    plt.plot(xvec, res, '-ko')
    plt.plot([0, 1], [0,0], '--r')

if __name__ == '__main__':
    test_ma()






