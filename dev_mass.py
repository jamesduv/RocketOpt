#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: jd
contact: jamesduv@umich.edu
affiliation: University of Michigan, Department of Aerospace Eng., CASLAB
"""

import numpy as np
import matplotlib.pyplot as plt

def comp_mass(design_vars, material_prop, geom_out):
    rho     = material_prop['rho']
    rc      = design_vars['rc']
    tw      = design_vars['tw']
    Lc      = design_vars['Lc']
    beta_nc = design_vars['beta_nc']
    beta_nd = design_vars['beta_nd']

    rt      = geom_out['rt']
    re      = geom_out['re']
    Lnc     = geom_out['Lnc']
    Lnd     = geom_out['Lnd']

    #combustion chamber mass
    mc      = rho * Lc * np.pi * (2 * rc * tw + tw**2)

    mnc     = ((rho * np.pi * (rc - rt) * tw) / (np.tan(beta_nc * np.pi / 180))) * (rc + rt + tw)

    mnd     = ((rho * np.pi * (re - rt) * tw) / (np.tan(beta_nd * np.pi / 180))) * (re + rt + tw)

    mtot = mc + mnc + mnd

    mass_out = {
        'mc'    : mc,
        'mnc'   : mnc,
        'mnd'   : mnd,
        'mtot'  : mtot
        }

    return mass_out

def cone_frustum_check():
    L       = 1
    r2i     = 0.5
    r1i     = 0.25
    t       = 0.01

    r2o     = r2i + t
    r1o     = r1i + t

    Vo = (np.pi / 3) * L * (r1o**2 + r2o**2 + r1o*r2o)
    Vi = (np.pi / 3) * L * (r1i**2 + r2i**2 + r1i*r2i)

    Vtot_1 = Vo - Vi

    beta = np.arctan((r2i-r1i) / L)

    Vtot_2 = ((np.pi * (r2i - r1i) * t) / np.tan(beta)) * (r1i+r2i+t)
    dif = Vtot_1 - Vtot_2
    test = 1


if __name__ == '__main__':
    cone_frustum_check()