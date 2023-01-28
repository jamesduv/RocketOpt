#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: jd
contact: jamesduv@umich.edu
affiliation: University of Michigan, Department of Aerospace Eng., CASLAB
"""


'''
Brainstorming XDSM

Component 1: Thermochemical solver
    Inputs:
        Shared design vars:
            pc  : combustion chamber pressure
            phi : equivalence ratio

        Local inputs:

        Outputs:
            Full thermochemical state post combustion
            pc, Tc, gamma, R,  rhoc, mdot_tot, mdot_ox

Component 2: Geometry generator
    Inputs:
        Design vars
            rc      : combustion chamber radius
            mdot_f  : fuel mass flow rate
            phi     : equivalence ratio
        Computed values:
            pc, Tc, gamma, R,  rhoc, mdot_tot, mdot_ox

    Outputs:
        geomteric design vars


Component 3: Mesh generator
    Inputs:
        geomteric design vars

    Outputs:
        Coordinates for all mesh points


Component 4: Flow model
    Inputs:
        geomteric design vars
        mesh coordinates

    Outputs:
        T, p, rho, M at all z coordinates
        

'''
import sys, os, pickle
import numpy as np
import scipy.optimize as opt

from dev_thermo import comp_thermochem
from dev_geom import comp_geom
from dev_mesh import comp_mesh
from dev_flow import comp_flow
from dev_heat import comp_heat
from dev_stress import comp_stress
from dev_mass import comp_mass
from dev_performance import comp_performance

import plots


def dev_opt():
    global XHIST
    XHIST = []

    save_dir = os.path.join('model_output')
    if not os.path.exists(save_dir):
        os.mkdir(save_dir)

    x0  = np.array([0.35, 0.35, 20, 50, 15, 20, 1.0, 15, 0.01, 0.001,  0.001, 0.1])

    Ispmin = 250
    Thrustmin = 200e+03
    cIsp = {
        'type'  : 'ineq',
        'fun'   : constr_Isp,
        'args'  : (Ispmin,)
            }

    cThrust = {
        'type'  : 'ineq',
        'fun'   : constr_thrust,
        'args'  : (Thrustmin,)
            }


    cStress = {
        'type'  : 'ineq',
        'fun'   : constr_stress,
               }
    con = [cIsp, cThrust]
    con2 = [cIsp, cThrust, cStress]

    #original - working for Isp constraint
    # bounds = [(0.5, 0.51),          #chamber radius
    #           (0.4, 0.41),         #chamber length
    #           (20, 21),      #chamber pressure
    #           (50, 51),     #beta_nc
    #           (15, 16),     #beta_nd
    #           (20, 21),      #fuel mass flow rate
    #           (0.8, 1.2),   #equivalence ratio
    #           (15,20),       #nozzle expansion ratio
    #           (0.001, 0.0015),     #inner wall thickness
    #           (0.005, 0.0051),    #channel thickness
    #           (0.001, 0.0015),   #shell thickness)
    #           (0.05, 0.5)]     #fuel cooling  fraction

    bounds = [(0.02, 0.51),          #chamber radius
          (0.1, 0.41),         #chamber length
          (5, 40),      #chamber pressure
          (50, 55),     #beta_nc
          (15, 16),     #beta_nd
          (10, 40),      #fuel mass flow rate
          (0.8, 1.2),   #equivalence ratio
          (2, 30),       #nozzle expansion ratio
          (0.005, 0.05),     #inner wall thickness
          (0.001, 0.05),    #channel thickness
          (0.001, 0.0015),   #shell thickness)
          (0.01, 0.5)]     #fuel cooling  fraction
    fmin = fobj

    res = opt.minimize(fobj, x0, method='SLSQP', bounds=bounds, constraints=con2, callback = callback)

    xopt = res['x']
    therm_out, geom_out, mesh_out, flow_out, heat_out, mass_out, stress_out, perf_out = fmodel(xopt)
    material_prop, mesh_param, flow_bc = get_prop_param_bc()
    del therm_out['gas']
    ##Update to add thermal properties etc
    data = {
        'therm_out'     : therm_out,
        'geom_out'      : geom_out,
        'mesh_out'      : mesh_out,
        'flow_out'      : flow_out,
        'heat_out'      : heat_out,
        'mass_out'      : mass_out,
        'stress_out'    : stress_out,
        'perf_out'      : perf_out,
        'xopt'          : xopt,
        'opt_res'       : res,
        'bounds'        : bounds,
        'constraints'   : con,
        'xhist'         : XHIST.copy(),
        'x0'            : x0,
        'material_prop' : material_prop,
        'mesh_param'    : mesh_param,
        'flow_bc'       : flow_bc
        }

    fn_save = os.path.join(save_dir, 'fisp_stress_f{:1.3e}_isp{:1.3e}.pickle'.format(Thrustmin, Ispmin))
    pickle.dump(data, open(fn_save, 'wb'))

    return res, therm_out, geom_out, mesh_out, flow_out, heat_out, mass_out, stress_out, perf_out

# def dev_opt():
#     ## works but not very interesting
#     x0  = np.array([0.5, 0.4, 20, 50, 15, 20, 1.0, 15, 0.001, 0.005,  0.001, 0.1])

#     bounds = [(0.5, 0.51),          #chamber radius
#               (0.4, 0.41),         #chamber length
#               (20, 21),      #chamber pressure
#               (50, 51),     #beta_nc
#               (15, 16),     #beta_nd
#               (20, 21),      #fuel mass flow rate
#               (0.8, 1.2),   #equivalence ratio
#               (15,20),       #nozzle expansion ratio
#               (0.001, 0.0015),     #inner wall thickness
#               (0.005, 0.0051),    #channel thickness
#               (0.001, 0.0015),   #shell thickness)
#               (0.05, 0.5)]     #fuel cooling  fraction
#     fmin = fobj
#     res = opt.minimize(fobj, x0, method='Powell', bounds=bounds)
#     return res

def callback(x):
    global XHIST
    XHIST.append(x)

def fobj(x):

    design_vars = {
        'rc'        : x[0],      # radius, combustion chamber [m]
        'Lc'        : x[1],      # length, combustion chamber [m]
        'pc'        : x[2],       # pressure, combustion chamber [atm]
        'beta_nc'   : x[3],       # angle, converging nozzle [deg]
        'beta_nd'   : x[4],       # angle, diverging nozzle [deg]
        'mdot_f'    : x[6],      # fuel mass flow rate, [kg/s]
        'phi'       : x[6],      # equivalence ratio CH4-O2 [-]
        'er'        : x[7],       # nozzle expansion ratio
        'tw'        : x[8],     # wall thickness
        'tc'        : x[9],    # channel thickness
        'ts'        : x[10],     # shell thickness
        'f_cool'    : x[11]       #fraction of fuel to run through cooling loop
    }

    material_prop, mesh_param, flow_bc = get_prop_param_bc()

    cp_cool     = flow_bc['Pr_fuel'] * flow_bc['k_fuel'] / flow_bc['mu_fuel']
    flow_bc['cp_fuel'] = cp_cool



    therm_out   = comp_thermochem(design_vars, flow_bc)

    geom_out    = comp_geom(design_vars, therm_out)

    # mesh_out    = comp_mesh(design_vars, geom_out, mesh_param)

    mass_out    = comp_mass(design_vars, material_prop, geom_out)

    # flow_out    = comp_flow(design_vars, flow_bc, therm_out, geom_out, mesh_out, mesh_param)

    # perf_out    = comp_performance(design_vars, flow_bc, geom_out, therm_out, flow_out)

    # heat_out    = comp_heat(design_vars, flow_bc, material_prop, mesh_param,
                            # therm_out, geom_out, mesh_out, flow_out)

    # stress_out  = comp_stress(design_vars, mesh_param, material_prop, mesh_out, flow_out, heat_out)

    mass = mass_out['mtot']

    return mass

def fmodel(x):
    design_vars = {
        'rc'        : x[0],      # radius, combustion chamber [m]
        'Lc'        : x[1],      # length, combustion chamber [m]
        'pc'        : x[2],       # pressure, combustion chamber [atm]
        'beta_nc'   : x[3],       # angle, converging nozzle [deg]
        'beta_nd'   : x[4],       # angle, diverging nozzle [deg]
        'mdot_f'    : x[6],      # fuel mass flow rate, [kg/s]
        'phi'       : x[6],      # equivalence ratio CH4-O2 [-]
        'er'        : x[7],       # nozzle expansion ratio
        'tw'        : x[8],     # wall thickness
        'tc'        : x[9],    # channel thickness
        'ts'        : x[10],     # shell thickness
        'f_cool'    : x[11]       #fraction of fuel to run through cooling loop
    }

    material_prop, mesh_param, flow_bc = get_prop_param_bc()
    cp_cool     = flow_bc['Pr_fuel'] * flow_bc['k_fuel'] / flow_bc['mu_fuel']
    flow_bc['cp_fuel'] = cp_cool

    therm_out   = comp_thermochem(design_vars, flow_bc)

    geom_out    = comp_geom(design_vars, therm_out)

    mesh_out    = comp_mesh(design_vars, geom_out, mesh_param)

    mass_out    = comp_mass(design_vars, material_prop, geom_out)

    flow_out    = comp_flow(design_vars, flow_bc, therm_out, geom_out, mesh_out, mesh_param)

    perf_out    = comp_performance(design_vars, flow_bc, geom_out, therm_out, flow_out)

    heat_out    = comp_heat(design_vars, flow_bc, material_prop, mesh_param,
                            therm_out, geom_out, mesh_out, flow_out)

    stress_out  = comp_stress(design_vars, mesh_param, material_prop, mesh_out, flow_out, heat_out)

    return therm_out, geom_out, mesh_out, flow_out, heat_out, mass_out, stress_out, perf_out


def constr_Isp(x, Ispmin):
    design_vars = {
        'rc'        : x[0],      # radius, combustion chamber [m]
        'Lc'        : x[1],      # length, combustion chamber [m]
        'pc'        : x[2],       # pressure, combustion chamber [atm]
        'beta_nc'   : x[3],       # angle, converging nozzle [deg]
        'beta_nd'   : x[4],       # angle, diverging nozzle [deg]
        'mdot_f'    : x[6],      # fuel mass flow rate, [kg/s]
        'phi'       : x[6],      # equivalence ratio CH4-O2 [-]
        'er'        : x[7],       # nozzle expansion ratio
        'tw'        : x[8],     # wall thickness
        'tc'        : x[9],    # channel thickness
        'ts'        : x[10],     # shell thickness
        'f_cool'    : x[11]       #fraction of fuel to run through cooling loop
    }

    material_prop, mesh_param, flow_bc = get_prop_param_bc()
    cp_cool     = flow_bc['Pr_fuel'] * flow_bc['k_fuel'] / flow_bc['mu_fuel']
    flow_bc['cp_fuel'] = cp_cool

    therm_out   = comp_thermochem(design_vars, flow_bc)

    geom_out    = comp_geom(design_vars, therm_out)

    mesh_out    = comp_mesh(design_vars, geom_out, mesh_param)

    # mass_out    = comp_mass(design_vars, material_prop, geom_out)

    flow_out    = comp_flow(design_vars, flow_bc, therm_out, geom_out, mesh_out, mesh_param)

    perf_out    = comp_performance(design_vars, flow_bc, geom_out, therm_out, flow_out)

    # heat_out    = comp_heat(design_vars, flow_bc, material_prop, mesh_param,
    #                         therm_out, geom_out, mesh_out, flow_out)

    # stress_out  = comp_stress(design_vars, mesh_param, material_prop, mesh_out, flow_out, heat_out)

    Isp = perf_out['Isp']
    val = Isp - Ispmin
    return val

def constr_thrust(x, thrustmin):
    design_vars = {
        'rc'        : x[0],      # radius, combustion chamber [m]
        'Lc'        : x[1],      # length, combustion chamber [m]
        'pc'        : x[2],       # pressure, combustion chamber [atm]
        'beta_nc'   : x[3],       # angle, converging nozzle [deg]
        'beta_nd'   : x[4],       # angle, diverging nozzle [deg]
        'mdot_f'    : x[6],      # fuel mass flow rate, [kg/s]
        'phi'       : x[6],      # equivalence ratio CH4-O2 [-]
        'er'        : x[7],       # nozzle expansion ratio
        'tw'        : x[8],     # wall thickness
        'tc'        : x[9],    # channel thickness
        'ts'        : x[10],     # shell thickness
        'f_cool'    : x[11]       #fraction of fuel to run through cooling loop
    }

    material_prop, mesh_param, flow_bc = get_prop_param_bc()
    cp_cool     = flow_bc['Pr_fuel'] * flow_bc['k_fuel'] / flow_bc['mu_fuel']
    flow_bc['cp_fuel'] = cp_cool

    therm_out   = comp_thermochem(design_vars, flow_bc)

    geom_out    = comp_geom(design_vars, therm_out)

    mesh_out    = comp_mesh(design_vars, geom_out, mesh_param)

    # mass_out    = comp_mass(design_vars, material_prop, geom_out)

    flow_out    = comp_flow(design_vars, flow_bc, therm_out, geom_out, mesh_out, mesh_param)

    perf_out    = comp_performance(design_vars, flow_bc, geom_out, therm_out, flow_out)

    # heat_out    = comp_heat(design_vars, flow_bc, material_prop, mesh_param,
    #                         therm_out, geom_out, mesh_out, flow_out)

    # stress_out  = comp_stress(design_vars, mesh_param, material_prop, mesh_out, flow_out, heat_out)

    thrust = perf_out['F']
    val = thrust - thrustmin
    return val

def constr_stress(x):
    design_vars = {
        'rc'        : x[0],      # radius, combustion chamber [m]
        'Lc'        : x[1],      # length, combustion chamber [m]
        'pc'        : x[2],       # pressure, combustion chamber [atm]
        'beta_nc'   : x[3],       # angle, converging nozzle [deg]
        'beta_nd'   : x[4],       # angle, diverging nozzle [deg]
        'mdot_f'    : x[6],      # fuel mass flow rate, [kg/s]
        'phi'       : x[6],      # equivalence ratio CH4-O2 [-]
        'er'        : x[7],       # nozzle expansion ratio
        'tw'        : x[8],     # wall thickness
        'tc'        : x[9],    # channel thickness
        'ts'        : x[10],     # shell thickness
        'f_cool'    : x[11]       #fraction of fuel to run through cooling loop
    }

    material_prop, mesh_param, flow_bc = get_prop_param_bc()

    cp_cool     = flow_bc['Pr_fuel'] * flow_bc['k_fuel'] / flow_bc['mu_fuel']
    flow_bc['cp_fuel'] = cp_cool

    therm_out   = comp_thermochem(design_vars, flow_bc)

    geom_out    = comp_geom(design_vars, therm_out)

    mesh_out    = comp_mesh(design_vars, geom_out, mesh_param)

    # mass_out    = comp_mass(design_vars, material_prop, geom_out)

    flow_out    = comp_flow(design_vars, flow_bc, therm_out, geom_out, mesh_out, mesh_param)

    perf_out    = comp_performance(design_vars, flow_bc, geom_out, therm_out, flow_out)

    heat_out    = comp_heat(design_vars, flow_bc, material_prop, mesh_param,
                            therm_out, geom_out, mesh_out, flow_out)

    stress_out  = comp_stress(design_vars, mesh_param, material_prop, mesh_out, flow_out, heat_out)


    ### original method. Failed to find descent direction in line search, likely
    #too stiff or discontinuous
    # stot = abs(stress_out['s_tot'])
    # smax    = stot.max()

    # syield = material_prop['s_yield']
    # val = syield - smax

    ### new method: measure the l2 norm of the difference between
    nz_glob     = mesh_out['zones']['global']['zc'].shape[0]

    sy_vec = np.repeat([material_prop['s_yield']], repeats= nz_glob)
    s_tot = stress_out['s_tot']
    s_tot_abs = abs(s_tot)
    dif = sy_vec - s_tot_abs
    val = dif.sum()
    return val


def get_prop_param_bc():
    mesh_param = {
        'nr'        :   10,
        'nz_c'      :   50,
        'nz_nc'     :   50,
        'nz_t'      :   3,
        'nz_nd'     :   100,
        }

    material_prop = {
        'k_w'   : 45,     #thermal conductivity, wall
        'k_s'   : 45,     # thermal conductivity, shell
        'E_w'   : 200E+09,  #youngs modulus steel
        'v_w'   : 0.29,      # poisson's ratio
        'cte_w' : 11e-06,    #coefficient of thermal expansion
        's_yield' : 350E+06,     #steel yield strength, Pa
        'rho'   : 8E+03     #steel density
        }

    flow_bc = {
        'Tfuel'     : 200,          #Kelvin, boiling point of ~110K, pre combustion temperature
        'p_inf'      : 0.1,          #ambient pressure, atm
        'T_inf'     : 200,          #K
        'rho_fuel'  : 260,          # kg/m^3 (liq 440, supercrit ~260)
        'k_fuel'    : 40e-03,     #methane thermal conductivity, sliquid ~200e-03, supercritical ~40e-03-80e-03
        'Pr_fuel'   : 0.86,         #(0.86)fuel Prandtl number, estimate
        'mu_fuel'   : 15e-06,      #fuel dynamic viscosity - supercritical
        'h_inf'     : 10,
        'Tcool'     : 90        #coolant inlet temperature
    }
    return material_prop, mesh_param, flow_bc



def dev_model():

    design_vars = {
        'rc'        : 0.2,      # radius, combustion chamber [m]
        'Lc'        : 0.4,      # length, combustion chamber [m]
        'pc'        : 50,       # pressure, combustion chamber [atm]
        'beta_nc'   : 60,       # angle, converging nozzle [deg]
        'beta_nd'   : 15,       # angle, diverging nozzle [deg]
        'mdot_f'    : 20,       # fuel mass flow rate, [kg/s]
        'phi'       : 0.6,      # equivalence ratio CH4-O2 [-]
        'er'        : 15,       # nozzle expansion ratio
        'tw'        : 0.01,     # wall thickness
        'tc'        : 0.005,    # channel thickness
        'ts'        : 0.003,     # shell thickness
        'f_cool'    : 0.1       #fraction of fuel to run through cooling loop
    }

    design_vars = {
        'rc'        : 0.51,      # radius, combustion chamber [m]
        'Lc'        : 0.4,      # length, combustion chamber [m]
        'pc'        : 21,       # pressure, combustion chamber [atm]
        'beta_nc'   : 50,       # angle, converging nozzle [deg]
        'beta_nd'   : 15,       # angle, diverging nozzle [deg]
        'mdot_f'    : 20,       # fuel mass flow rate, [kg/s]
        'phi'       : 1.2,      # equivalence ratio CH4-O2 [-]
        'er'        : 20,       # nozzle expansion ratio
        'tw'        : 0.001,     # wall thickness
        'tc'        : 0.005,    # channel thickness
        'ts'        : 0.001,     # shell thickness
        'f_cool'    : 0.1       #fraction of fuel to run through cooling loop
    }

    flow_bc = {
        'Tfuel'     : 200,          #Kelvin, boiling point of ~110K, pre combustion temperature
        'p_inf'      : 0.1,          #ambient pressure, atm
        'T_inf'     : 200,          #K
        'rho_fuel'  : 260,          # kg/m^3 (liq 440, supercrit ~260)
        'k_fuel'    : 40e-03,     #methane thermal conductivity, sliquid ~200e-03, supercritical ~40e-03-80e-03
        'Pr_fuel'   : 0.86,         #(0.86)fuel Prandtl number, estimate
        'mu_fuel'   : 15e-06,      #fuel dynamic viscosity - supercritical
        'h_inf'     : 10,
        'Tcool'     : 90        #coolant inlet temperature
    }
    cp_cool     = flow_bc['Pr_fuel'] * flow_bc['k_fuel'] / flow_bc['mu_fuel']
    flow_bc['cp_fuel'] = cp_cool

    mesh_param = {
        'nr'        :   10,
        'nz_c'      :   50,
        'nz_nc'     :   50,
        'nz_t'      :   3,
        'nz_nd'     :   100,
        }

    material_prop = {
        'k_w'   : 45,     #thermal conductivity, wall
        'k_s'   : 45,     # thermal conductivity, shell
        'E_w'   : 200E+09,  #youngs modulus steel
        'v_w'   : 0.29,      # poisson's ratio
        'cte_w' : 11e-06,    #coefficient of thermal expansion
        's_yield' : 350E+06,     #steel yield strength, Pa
        'rho'   : 8E+03     #steel density
        }

    therm_out   = comp_thermochem(design_vars, flow_bc)

    geom_out    = comp_geom(design_vars, therm_out)

    mesh_out    = comp_mesh(design_vars, geom_out, mesh_param)

    mass_out    = comp_mass(design_vars, material_prop, geom_out)

    flow_out    = comp_flow(design_vars, flow_bc, therm_out, geom_out, mesh_out, mesh_param)

    perf_out    = comp_performance(design_vars, flow_bc, geom_out, therm_out, flow_out)

    heat_out    = comp_heat(design_vars, flow_bc, material_prop, mesh_param,
                            therm_out, geom_out, mesh_out, flow_out)

    stress_out  = comp_stress(design_vars, mesh_param, material_prop, mesh_out, flow_out, heat_out)

    plots.plot_flow_vs_z(flow_out, mesh_out)
    plots.plot_stress(stress_out, mesh_out, material_prop)
    plots.plot_temp_contourf(heat_out, mesh_out, mesh_param)

    return therm_out, geom_out, mesh_out, flow_out, heat_out, mass_out, stress_out, perf_out


if __name__ == '__main__':
    therm_out, geom_out, mesh_out, flow_out, heat_out, mass_out, stress_out, perf_out = dev_model()
    # res, therm_out, geom_out, mesh_out, flow_out, heat_out, mass_out, stress_out, perf_out = dev_opt()