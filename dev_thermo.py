#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: jd
contact: jamesduv@umich.edu
affiliation: University of Michigan, Department of Aerospace Eng., CASLAB
"""
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt


def comp_thermochem(design_vars, flow_bc):
    ''' Thermochemical model component '''

    mdot_f  = design_vars['mdot_f']
    equiv   = design_vars['phi']

    #Initialize canter gas Solution object, use gri30 methane-air mechanism
    gas     = ct.Solution('gri30.xml', 'gri30_mix')
    PATM    = ct.one_atm

    #boundary/initial conditions
    T0          = flow_bc['Tfuel']
    pinf        = flow_bc['p_inf'] * PATM
    pc          = design_vars['pc'] * PATM

    #fuel-oxidizer mass ratio at stoichiometric conditions
    fo_stoich   = get_fostoich_methane()
    mdot_ox     = (mdot_f / equiv) / fo_stoich
    mdot_tot    = mdot_f + mdot_ox

    Y_fuel      = mdot_f / mdot_tot
    Y_ox        = mdot_ox / mdot_tot
    composition     = {'CH4':Y_fuel, 'O2':Y_ox}

    #set composition
    gas.Y       = composition
    gas.TP      = T0, pc

    #combustion at constant pressure
    gas.equilibrate('HP')
    pc      = gas.P
    Tc      = gas.T
    rhoc    = gas.density
    gamma   = gas.cp_mass / gas.cv_mass
    R       = ct.gas_constant
    MW      = gas.mean_molecular_weight
    R_sp    = R / MW

    #compute gas viscosity
    mumix           = get_gas_viscosity(gas)
    thermal_cond    = gas.thermal_conductivity
    Pr              = mumix * gas.cp_mass / thermal_cond

    # gas.TP = 2000, pc
    # k2 =gas.thermal_conductivity
    # mumix2      = get_gas_viscosity(gas)
    # gamma2      = gas.cp_mass / gas.cv_mass

    thermo_out = {
        'pc'        : pc,
        'Tc'        : Tc,
        'rhoc'      : rhoc,
        'gamma'     : gamma,
        'R'         : R,
        'R_sp'      : R_sp,
        'mdot_tot'  : mdot_tot,
        'mdot_ox'   : mdot_ox,
        'gas'       : gas,
        'mu'        : mumix,    #dynamic viscosity
        'k'         : thermal_cond,
        'Pr'        : Pr
        }

    return thermo_out

def get_gas_viscosity(gas):
    X       = gas.X
    MW      = gas.molecular_weights
    visc    = gas.species_viscosities
    nX      = X.shape[0]

    mumix   = 0
    num_ls  = []
    den_ls  = []
    for ii in np.arange(nX):
        num     = X[ii] * visc[ii]
        Mi      = MW[ii]
        mui     = visc[ii]
        den     = 0
        for jj in np.arange(nX):
            Mj      = MW[jj]
            Xj      = X[jj]
            muj     = visc[jj]
            curnum  = (1 + (mui/muj)**(0.5) * (Mj/Mi)**(0.25))**2
            curden  = (8 + 8 * (Mi/Mj))**(0.5)
            phi_ij  = curnum / curden
            den += Xj * phi_ij
        cur_term = num / den
        mumix += cur_term
    return mumix


def dev_thermochem():

    #get gas1 for thermochemical properties
    gas = ct.Solution('gri30.xml')
    PATM    = 101325


    rc  = 0.1
    lc  = 0.3
    expansion_ratio = 25 #nozzle expansion ratio


    #thermochemical parameters
    mdot_f      = 2.
    equiv       = 1.    #equivalence ratio
    pc          = 10 * PATM

    #boundary/initial conditions
    T0          = 100 #initial mixture temperature
    pinf        = 0.5 * PATM

    #fuel-oxidizer mass ratio at stoichiometric conditions
    fo_stoich   = get_fostoich_methane()
    mdot_ox     = (mdot_f / equiv) / fo_stoich
    mdot_tot    = mdot_f + mdot_ox
    composition     = {'CH4':mdot_f, 'O2':mdot_ox}

    #set initial gas composition and state
    gas.Y   = composition
    gas.TP = T0, pc

    #combustion at constant pressure
    gas.equilibrate('HP')
    Tc  = gas.T
    rhoc = gas.density
    gamma = gas.cp_mass / gas.cv_mass


def get_fostoich_methane():
    mwo2    = 2 * 15.999
    mwc     = 12.0107
    mwh     = 1.00784
    mwch4   = mwc + 4 * mwh

    fo_stoich   = mwch4 / (2 * mwo2)
    return fo_stoich


if __name__ == '__main__':
    dev_thermochem()
