import sys
import os
import pickle
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

def explore():
    fn_read = os.path.join('model_output_norm', 'thrust3.000e+05_mass.pickle')
    data = pickle.load(open(fn_read, 'rb'))
    limits = data['limits']
    xopt = data['opt_res']['x']
    xopt[0] = 1
    model_out = fmodel(xopt, limits)
    plots.plot_temp_contourf(model_out['heat_out'], model_out['mesh_out'], model_out['mesh_param'])
    plots.plot_flow_vs_z(model_out['flow_out'], model_out['mesh_out'])
    plots.plot_stress(model_out['stress_out'], model_out['mesh_out'], model_out['material_prop'])

def callback(x):
    global XHIST
    XHIST.append(x)

def dev_opt_norm(is_save_res = True):
    global XHIST
    XHIST = []

    save_dir = os.path.join('model_output_norm')
    if not os.path.exists(save_dir):
        os.mkdir(save_dir)

    ndv = 12
    x0 = np.ones(ndv) * 0.5

    bounds = [  (0, 1),          #chamber radius
                (0, 1),           #chamber length
                (0, 1),      #chamber pressure, atm
                (0, 1),     #beta_nc
                (0, 1),     #beta_nd
                (0, 1),      #fuel mass flow rate
                (0, 1),   #equivalence ratio
                (0, 1),       #nozzle expansion ratio
                (0, 1),     #inner wall thickness
                (0, 1),    #channel thickness
                (0, 1),   #shell thickness)
                (0, 1)]     #fuel cooling  fraction
    gamma   = 1.2
    R       = 370
    Tc      = 3600
    pc      = 80 * 101325
    mdot_o  = 20 / 1.2 / 0.25
    mdot    =  20 + mdot_o

    limits = {
        'rc'        : (0.3, 1.0),         # radius, combustion chamber [m]
        'Lc'        : (0.3, 0.40),          # length, combustion chamber [m]
        'pc'        : (100, 350),              # pressure, combustion chamber [atm]
        'beta_nc'   : (30, 60),             # angle, converging nozzle [deg]
        'beta_nd'   : (5, 25),             # angle, diverging nozzle [deg]
        'mdot_f'    : (2, 300),             # fuel mass flow rate, [kg/s]
        'phi'       : (0.6, 1.4),           # equivalence ratio CH4-O2 [-]
        'er'        : (10, 40),              # nozzle expansion ratio
        'tw'        : (0.005, 0.05),        # wall thickness
        'tc'        : (0.001, 0.05),        # channel thickness
        'ts'        : (0.001, 0.0015),      # shell thickness
        'f_cool'    : (0.001, 0.5)           #fraction of fuel to run through cooling loop
    }
    Ispmin = 290
    Thrustmin = 1500e+03
    Tmax  = 700
    cIsp = {
        'type'  : 'ineq',
        'fun'   : constr_Isp,
        'args'  : (limits, Ispmin,)
            }

    cThrust = {
        'type'  : 'ineq',
        'fun'   : constr_thrust,
        'args'  : (limits, Thrustmin,)
            }
    cStress = {
        'type'  : 'ineq',
        'fun'   : constr_stress,
        'args'  : (limits,)
                }
    cTemp = {
        'type'  : 'ineq',
        'fun'   : constr_temp,
        'args'  : (limits, Tmax)
                }

    con = [cThrust, cIsp, cStress, cTemp]

    args = (limits,)
    fobj = fobj_mass

    res = opt.minimize(fobj, x0, args = args, method='SLSQP', bounds=bounds,
                       constraints = con, callback = callback)

    xopt = res['x']

    model_out = fmodel(xopt, limits)
    plots.plot_temp_contourf(model_out['heat_out'], model_out['mesh_out'], model_out['mesh_param'])
    plots.plot_flow_vs_z(model_out['flow_out'], model_out['mesh_out'])
    plots.plot_stress(model_out['stress_out'], model_out['mesh_out'], model_out['material_prop'])
    xopt_sc = design_vars_vec2dict(xopt, limits)

    #delet non-pickleable cantera gas object from model output
    del model_out['therm_out']['gas']

    data = {
        'model_out'     : model_out,
        'xopt_scaled'   : xopt_sc,
        'opt_res'       : res,
        'x0'            : x0,
        'bounds'        : bounds,
        'limits'        : limits,
        'xhist'         : XHIST.copy()
        }
    if is_save_res:
        # fn_save = os.path.join(save_dir, 'isp{:1.0f}_thrust{:1.3e}_stress_Temp_mass.pickle'.format(Ispmin, Thrustmin))
        # fn_save = os.path.join(save_dir, 'thrust{:1.3e}_mass_higherplim.pickle'.format(Thrustmin))
        fn_save = os.path.join(save_dir, 'raptor_higherpc.pickle'.format(Thrustmin))
        pickle.dump(data, open(fn_save, 'wb'))

    return xopt_sc, res, model_out, limits

def opt_unc():
    '''Unconstrained mass optimization'''
    global XHIST
    XHIST = []

    save_dir = os.path.join('model_output_norm')
    if not os.path.exists(save_dir):
        os.mkdir(save_dir)

    ndv = 12
    x0 = np.ones(ndv) * 0.5

    bounds = [(0, 1),          #chamber radius
          (0, 1),           #chamber length
          (0, 1),      #chamber pressure, atm
          (0, 1),     #beta_nc
          (0, 1),     #beta_nd
          (0, 1),      #fuel mass flow rate
          (0, 1),   #equivalence ratio
          (0, 1),       #nozzle expansion ratio
          (0, 1),     #inner wall thickness
          (0, 1),    #channel thickness
          (0, 1),   #shell thickness)
          (0, 1)]     #fuel cooling  fraction

    limits = {
        'rc'        : (0.1, 0.5),         # radius, combustion chamber [m]
        'Lc'        : (0.1, 0.4),          # length, combustion chamber [m]
        'pc'        : (5, 40),              # pressure, combustion chamber [atm]
        'beta_nc'   : (30, 60),             # angle, converging nozzle [deg]
        'beta_nd'   : (15, 25),             # angle, diverging nozzle [deg]
        'mdot_f'    : (10, 11),             # fuel mass flow rate, [kg/s]
        'phi'       : (0.8, 1.2),           # equivalence ratio CH4-O2 [-]
        'er'        : (2, 30),              # nozzle expansion ratio
        'tw'        : (0.005, 0.05),        # wall thickness
        'tc'        : (0.001, 0.05),        # channel thickness
        'ts'        : (0.001, 0.0015),      # shell thickness
        'f_cool'    : (0.01, 0.5)           #fraction of fuel to run through cooling loop
    }

    args = (limits,)
    fobj = fobj_mass
    res = opt.minimize(fobj, x0, args = args, method='SLSQP', bounds=bounds, callback = callback)
    xopt = res['x']
    model_out = fmodel(xopt, limits)
    xopt_sc = design_vars_vec2dict(xopt, limits)

    #delet non-pickleable cantera gas object from model output
    del model_out['therm_out']['gas']

    data = {
        'model_out'     : model_out,
        'xopt_scaled'   : xopt_sc,
        'opt_res'       : res,
        'x0'            : x0,
        'bounds'        : bounds,
        'limits'        : limits,
        'xhist'         : XHIST.copy()
        }

    fn_save = os.path.join(save_dir, 'unconstrained_mass.pickle')
    pickle.dump(data, open(fn_save, 'wb'))

    return xopt_sc, res, model_out, limits

def design_vars_vec2dict(x, limits):
    '''Convert the normalized input vector to a dictionary of scaled, dimensional values
    
    Args:
        x (ndarray)     : vector of normalized (0-1) design variables
        limits (dict)   : named limit tuples, (minval, maxval) for each entry of the design variables, 
            to recover dimensional quantities
    Returns:
        phys_vals (dict)    : the design variables in dimensional quantities
    '''

    norm = {
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

    phys_vals = {}
    for k,v in norm.items():
        xmin, xmax = limits[k]
        phys_vals[k] = xmin + (v * (xmax - xmin))

    return phys_vals

def fobj_mass(mu, limits):
    '''Objective function - mass
    
    Args:
        x (ndarray)     : vector of normalized (0-1) design variables
        limits (dict)   : named limits, xmin, xmax for each entry of the design variables, to recover dimensional quantities

    Retuns:
        mass (float)    : the total structural mass
    '''

    design_vars = design_vars_vec2dict(mu, limits)

    material_prop, mesh_param, flow_bc = get_prop_param_bc()

    therm_out   = comp_thermochem(design_vars, flow_bc)

    geom_out    = comp_geom(design_vars, therm_out)

    mass_out    = comp_mass(design_vars, material_prop, geom_out)

    mass = mass_out['mtot']

    return mass

def constr_Isp(x, limits, Ispmin):
    '''Isp constraint function'''

    design_vars = design_vars_vec2dict(x, limits)
    material_prop, mesh_param, flow_bc = get_prop_param_bc()

    therm_out   = comp_thermochem(design_vars, flow_bc)

    geom_out    = comp_geom(design_vars, therm_out)

    mesh_out    = comp_mesh(design_vars, geom_out, mesh_param)

    flow_out    = comp_flow(design_vars, flow_bc, therm_out, geom_out, mesh_out, mesh_param)

    perf_out    = comp_performance(design_vars, flow_bc, geom_out, therm_out, flow_out)

    Isp = perf_out['Isp']
    val = Isp - Ispmin
    return val

def constr_thrust(x, limits, Thrustmin):
    '''Isp constraint function'''

    design_vars = design_vars_vec2dict(x, limits)
    material_prop, mesh_param, flow_bc = get_prop_param_bc()

    therm_out   = comp_thermochem(design_vars, flow_bc)

    geom_out    = comp_geom(design_vars, therm_out)

    mesh_out    = comp_mesh(design_vars, geom_out, mesh_param)

    flow_out    = comp_flow(design_vars, flow_bc, therm_out, geom_out, mesh_out, mesh_param)

    perf_out    = comp_performance(design_vars, flow_bc, geom_out, therm_out, flow_out)

    Thrust = perf_out['F']
    val = Thrust - Thrustmin
    return val

def constr_stress(x, limits):
    '''Stress constraint function'''

    design_vars = design_vars_vec2dict(x, limits)
    material_prop, mesh_param, flow_bc = get_prop_param_bc()

    therm_out   = comp_thermochem(design_vars, flow_bc)

    geom_out    = comp_geom(design_vars, therm_out)

    mesh_out    = comp_mesh(design_vars, geom_out, mesh_param)

    # mass_out    = comp_mass(design_vars, material_prop, geom_out)

    flow_out    = comp_flow(design_vars, flow_bc, therm_out, geom_out, mesh_out, mesh_param)

    perf_out    = comp_performance(design_vars, flow_bc, geom_out, therm_out, flow_out)

    heat_out    = comp_heat(design_vars, flow_bc, material_prop, mesh_param,
                            therm_out, geom_out, mesh_out, flow_out)

    stress_out  = comp_stress(design_vars, mesh_param, material_prop, mesh_out, flow_out, heat_out)


    stot = abs(stress_out['s_tot'])
    smax    = stot.max()

    syield = material_prop['s_yield']
    val = syield - smax

    return val

def constr_temp(x, limits, Tmax):
    '''Stress constraint function'''

    design_vars = design_vars_vec2dict(x, limits)
    material_prop, mesh_param, flow_bc = get_prop_param_bc()

    therm_out   = comp_thermochem(design_vars, flow_bc)

    geom_out    = comp_geom(design_vars, therm_out)

    mesh_out    = comp_mesh(design_vars, geom_out, mesh_param)

    # mass_out    = comp_mass(design_vars, material_prop, geom_out)

    flow_out    = comp_flow(design_vars, flow_bc, therm_out, geom_out, mesh_out, mesh_param)

    perf_out    = comp_performance(design_vars, flow_bc, geom_out, therm_out, flow_out)

    heat_out    = comp_heat(design_vars, flow_bc, material_prop, mesh_param,
                            therm_out, geom_out, mesh_out, flow_out)

    # stress_out  = comp_stress(design_vars, mesh_param, material_prop, mesh_out, flow_out, heat_out)

    Tw = heat_out['Tw']
    val = Tmax - Tw.max()

    return val

def fmodel(x, limits):
    '''Call the full model, given a vector of normalized inputs and a dictionary
    giving the physical bounds for each variable'''

    design_vars = design_vars_vec2dict(x, limits)
    material_prop, mesh_param, flow_bc = get_prop_param_bc()

    therm_out   = comp_thermochem(design_vars, flow_bc)

    geom_out    = comp_geom(design_vars, therm_out)

    mesh_out    = comp_mesh(design_vars, geom_out, mesh_param)

    mass_out    = comp_mass(design_vars, material_prop, geom_out)

    flow_out    = comp_flow(design_vars, flow_bc, therm_out, geom_out, mesh_out, mesh_param)

    perf_out    = comp_performance(design_vars, flow_bc, geom_out, therm_out, flow_out)

    heat_out    = comp_heat(design_vars, flow_bc, material_prop, mesh_param,
                            therm_out, geom_out, mesh_out, flow_out)

    stress_out  = comp_stress(design_vars, mesh_param, material_prop, mesh_out, flow_out, heat_out)

    res = {
        'therm_out'     : therm_out,
        'geom_out'      : geom_out,
        'mesh_out'      : mesh_out,
        'mass_out'      : mass_out,
        'flow_out'      : flow_out,
        'perf_out'      : perf_out,
        'heat_out'      : heat_out,
        'stress_out'    : stress_out,
        'material_prop' : material_prop,
        'mesh_param'    : mesh_param,
        'flow_bc'       : flow_bc
        }

    return res

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
    cp_cool     = flow_bc['Pr_fuel'] * flow_bc['k_fuel'] / flow_bc['mu_fuel']
    flow_bc['cp_fuel'] = cp_cool
    return material_prop, mesh_param, flow_bc

def test():
    bounds = [(0,1),]*12
    test = 1

if __name__ == '__main__':
    test()
    # xopt, res, model_out, limits = dev_opt_norm()
    # xopt, res, model_out, limits = opt_unc()
    # explore()