import numpy as np

def comp_stress(design_vars, mesh_param, material_prop, mesh_out, flow_out, heat_out):

    nr          = mesh_param['nr']
    p_gas       = flow_out['global']['p']
    Tw          = heat_out['Tw']
    p_ch        = flow_out['channel']['p']
    beta_vec    = mesh_out['zones']['global']['beta_vec']
    cte         = material_prop['cte_w']
    E           = material_prop['E_w']
    poisson     = material_prop['v_w']
    t_w         = design_vars['tw']

    rc_glob     = mesh_out['zones']['global']['rc']
    rc_wi       = rc_glob[:,0]
    rc_wo       = rc_glob[:,(nr-1)]
    rc_mean     = 0.5 * (rc_wi + rc_wo)

    ### stress from thermal expansion
    #varies linearly through thickness from -s_cte to +s_cte
    deltaT = Tw[:,-1] - Tw[:,0]
    s_cte   = 2 * cte * E * deltaT / (1 - poisson)

    ### hoop stress
    s_hoop = ((p_gas - p_ch) * rc_mean) / (t_w * np.cos(beta_vec * np.pi/180))
    s_tot   = s_cte + s_hoop

    stress_out = {
        's_hoop'    : s_hoop,
        's_cte'     : s_cte,
        's_tot'     : s_tot
        }

    return stress_out





