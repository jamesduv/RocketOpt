import numpy as np

def comp_performance(design_vars, flow_bc, geom_out, therm_out, flow_out):

    mdot    = therm_out['mdot_tot']
    Ae      = geom_out['Ae']
    At      = geom_out['At']
    beta_nd = design_vars['beta_nd']
    p_gas   = flow_out['global']['p']
    p_inf   = flow_bc['p_inf']
    u_gas   = flow_out['global']['u']

    ONEATM  = 101325
    p_inf   = p_inf * ONEATM

    #ideal values without conical nozzle correction
    F_ideal = ( mdot * u_gas[-1]) + (p_gas[-1] - p_inf) * Ae
    c_star  = p_gas[0] * At / mdot
    Isp_ideal = F_ideal / (mdot * 9.81)
    Cf_ideal = F_ideal / (At * p_gas[0])

    #values with conical nozzle correction
    lam     = 0.5 * (1 + np.cos(beta_nd * np.pi / 180))
    F_c     = (lam * mdot * u_gas[-1]) + (p_gas[-1] - p_inf) * Ae
    Isp_c   = F_c / (mdot * 9.81)
    Cf_c    = F_c / (At * p_gas[0])

    perf_out = {
        'F_ideal'   : F_ideal,
        'F'         : F_c,
        'Isp_ideal'     : Isp_ideal,
        'Isp'       : Isp_c,
        'Cf_ideal'  : Cf_ideal,
        'Cf'        : Cf_c,
        'lambda'    : lam,
        'c_star'    : c_star
        }

    return perf_out
