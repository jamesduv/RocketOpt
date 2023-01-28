#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: jd
contact: jamesduv@umich.edu
affiliation: University of Michigan, Department of Aerospace Eng., CASLAB
"""

import numpy as np

def comp_geom(design_vars, thermo_out):
    '''Geometry generator component'''

    gamma   = thermo_out['gamma']
    R      = thermo_out['R_sp']
    Tc     = thermo_out['Tc']
    pc     = thermo_out['pc']
    mdot   = thermo_out['mdot_tot']
    er     = design_vars['er']
    rc     = design_vars['rc']
    beta_nc    = design_vars['beta_nc']
    beta_nd    = design_vars['beta_nd']

    #compute the throat area and radius
    At1     = np.sqrt(gamma * R * Tc)
    At2     = np.sqrt( ( 2 / (gamma + 1))**((gamma+1) / (gamma-1))  )
    At3     = mdot / (pc * gamma)
    At      = (At1 * At3) / At2
    rt      = np.sqrt( At / np.pi)

    #exit area
    Ae      = er * At
    re      = np.sqrt( Ae / np.pi)

    #nozzle convergent/divergent section lengths, beta_nc is in deg
    Lnc     = (rc - rt) / np.tan(beta_nc * np.pi / 180)
    Lnd     = (re - rt) / np.tan(beta_nd * np.pi / 180)

    geom_out = {
        'At'        : At,
        'rt'        : rt,
        'Ae'        : Ae,
        're'        : re,
        'Lnc'       : Lnc,
        'Lnd'       : Lnd}

    return geom_out

