#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: jd
contact: jamesduv@umich.edu
affiliation: University of Michigan, Department of Aerospace Eng., CASLAB
"""
import numpy as np
import matplotlib.pyplot as plt

def comp_heat(design_vars, flow_bc, material_prop, mesh_param, therm_out, geom_out, mesh_out, flow_out):

    ### Extract geometry info
    rc          = design_vars['rc']
    nr          = mesh_param['nr']
    r_inner     = mesh_out['zones']['inner']['rc']
    z_inner     = mesh_out['zones']['inner']['zc']
    Ac_inner    = mesh_out['zones']['inner']['Ac']
    beta_nc     = design_vars['beta_nc']
    beta_nd     = design_vars['beta_nd']
    dz_vec      = mesh_out['zones']['global']['dz_vec']
    Tcool       = flow_bc['Tcool']
    cp_fuel     = flow_bc['cp_fuel']

    mdot_f      = design_vars['mdot_f']
    mdot_cool   = design_vars['f_cool'] * mdot_f

    ### Extract transport properties
    Pr          = therm_out['Pr']
    k_gas       = therm_out['k']
    Pr_fuel     = flow_bc['Pr_fuel']
    k_fuel      = flow_bc['k_fuel']
    T_inf       = flow_bc['T_inf']
    h_inf       = flow_bc['h_inf']

    #Extract solid conductivities
    k_w     = material_prop['k_w']
    k_s     = material_prop['k_s']

    ### Extract global mesh (reshaped into nz x nr)
    rc_glob     = mesh_out['zones']['global']['rc']
    zc_glob     = mesh_out['zones']['global']['zc']
    nodes_glob  = mesh_out['zones']['global']['nodes']
    nz_glob, nr_glob    = zc_glob.shape
    n_glob      = nz_glob * nr_glob

    row_start_nc    = mesh_out['zones']['global']['row_start_nc']
    row_t           = mesh_out['zones']['global']['row_t']
    row_start_nd    = mesh_out['zones']['global']['row_start_nd']

    dz_c        = zc_glob[1,0] - zc_glob[0,0]
    dz_nc       = zc_glob[row_start_nc+1,0] - zc_glob[row_start_nc,0]
    dz_nd       = zc_glob[row_start_nd+1,0] - zc_glob[row_start_nd,0]

    ### Extract global T, Re profiles (and channel)
    T_glob      = flow_out['global']['T']
    Re_glob     = flow_out['global']['Re']
    Re_ch       = flow_out['channel']['Re']
    dh_ch       = flow_out['channel']['hydraulic_d']
    Acc     = mesh_out['zones']['global']['Ac_channel']
    rco     = mesh_out['zones']['global']['rco']
    rci     = mesh_out['zones']['global']['rci']
    rso     = mesh_out['zones']['global']['rso']

    ### compute heat transfer coeffs and thermal resistances between gas flow
    # and inner wall along full engine length
    Rgas        = np.zeros(nz_glob)
    ff          = np.zeros(nz_glob)
    Nu_gas      = np.zeros(nz_glob)
    h_gas       = np.zeros(nz_glob)
    A_gas        = np.zeros(nz_glob)

    R_chi       = np.zeros(nz_glob)
    R_cho       = np.zeros(nz_glob)
    ff_ch       = np.zeros(nz_glob)
    Nu_ch       = np.zeros(nz_glob)
    h_ch        = np.zeros(nz_glob)
    A_chi        = np.zeros(nz_glob)
    A_cho        = np.zeros(nz_glob)

    A_sho       = np.zeros(nz_glob)
    R_sho       = np.zeros(nz_glob)

    Nu_gas2 = np.zeros(nz_glob)
    h_gas2 = np.zeros(nz_glob)

    # combustion chamber
    for ii in np.arange(row_start_nc):
        rcur    = rc_glob[ii,0]
        ff[ii], Nu_gas[ii], h_gas[ii] = gnielinski_correlation(Re_glob[ii], Pr, rcur, k_gas)

        Nu_gas2[ii], h_gas2[ii] = simple_corr(Re_glob[ii], Pr, rcur, k_gas)
        rch_cur = 0.5 * dh_ch[ii]   #use half of hydraulic diameter as is hard coded to use radius
        ff_ch[ii], Nu_ch[ii], h_ch[ii] = gnielinski_correlation(Re_ch[ii], Pr_fuel, rch_cur, k_fuel)

        #first node, only half dz_c in height
        if ii == 0:
            A_gas[ii]   = np.pi * rcur * dz_c
            A_chi[ii]   = np.pi * rci[ii] * dz_c
            A_cho[ii]   = np.pi * rco[ii] * dz_c
            A_sho[ii]   = np.pi * rso[ii] * dz_c
        elif ii == (row_start_nc-1):
            # gas
            A_gas[ii] = area_h_cyl_conv(r1 = rcur, rdown = rc_glob[ii+1,0],
                                     dz_straight = dz_c, dz_slant = dz_nc)

            # channel inner
            A_chi[ii] = area_h_cyl_conv(r1 = rci[ii], rdown = rci[ii+1],
                                     dz_straight = dz_c, dz_slant = dz_nc)

            # channel outer
            A_cho[ii] = area_h_cyl_conv(r1 = rco[ii], rdown = rco[ii+1],
                                     dz_straight = dz_c, dz_slant = dz_nc)

            # shell outer
            A_sho[ii] = area_h_cyl_conv(r1 = rso[ii], rdown = rso[ii+1],
                                     dz_straight = dz_c, dz_slant = dz_nc)
        else:
            A_gas[ii]   = 2 * np.pi * rcur * dz_c
            A_chi[ii]   = 2 * np.pi * rci[ii] * dz_c
            A_cho[ii]   = 2 * np.pi * rco[ii] * dz_c
            A_sho[ii]   = 2 * np.pi * rso[ii] * dz_c

        Rgas[ii]    = 1 / (h_gas[ii] * A_gas[ii])
        R_chi[ii]   = 1 / (h_ch[ii] * A_chi[ii])
        R_cho[ii]   = 1 / (h_ch[ii] * A_cho[ii])
        R_sho[ii]   = 1 / (h_inf * A_sho[ii])

    #nozzle converging section
    for ii in np.arange(row_start_nc, row_start_nd, 1):
        rcur    = rc_glob[ii,0]
        ff[ii], Nu_gas[ii], h_gas[ii] = gnielinski_correlation(Re_glob[ii], Pr, rcur, k_gas)

        Nu_gas2[ii], h_gas2[ii] = simple_corr(Re_glob[ii], Pr, rcur, k_gas)

        rch_cur = 0.5 * dh_ch[ii]   #use half of hydraulic diameter as is hard coded to use radius
        ff_ch[ii], Nu_ch[ii], h_ch[ii] = gnielinski_correlation(Re_ch[ii], Pr_fuel, rch_cur, k_fuel)

        #throat
        if ii == (row_start_nd-1):

            # gas side
            A_gas[ii] = area_h_throat(r1 = rcur, rup = rc_glob[ii-1,0], dz_up = dz_nc,
                                    rdown = rc_glob[ii+1,0], dz_down = dz_nd)

            #channel inner
            A_chi[ii] = area_h_throat(r1 = rci[ii], rup = rci[ii-1], dz_up = dz_nc,
                                    rdown = rci[ii+1], dz_down = dz_nd)

            #channel outer
            A_cho[ii] = area_h_throat(r1 = rco[ii], rup = rco[ii-1], dz_up = dz_nc,
                                    rdown = rco[ii+1], dz_down = dz_nd)

            # shell outer
            A_sho[ii] = area_h_throat(r1 = rso[ii], rup = rso[ii-1], dz_up = dz_nc,
                                    rdown = rso[ii+1], dz_down = dz_nd)
        else:
            # gas
            A_gas[ii] = area_h_slant(rcur = rcur, rup = rc_glob[ii-1,0],
                                     rdown = rc_glob[ii+1,0], dz = dz_nc)

            # Channel inner
            A_chi[ii] = area_h_slant(rcur = rci[ii], rup = rci[ii-1],
                                     rdown = rci[ii+1], dz = dz_nc)

            # Channel outer
            A_cho[ii] = area_h_slant(rcur = rco[ii], rup = rco[ii-1],
                                     rdown = rco[ii+1], dz = dz_nc)

            # shell outer
            A_sho[ii] = area_h_slant(rcur = rso[ii], rup = rso[ii-1],
                                     rdown = rso[ii+1], dz = dz_nc)

        Rgas[ii]    = 1 / (h_gas[ii] * A_gas[ii])
        R_chi[ii]   = 1 / (h_ch[ii] * A_chi[ii])
        R_cho[ii]   = 1 / (h_ch[ii] * A_cho[ii])
        R_sho[ii]   = 1 / (h_inf * A_sho[ii])

    #nozzle diverging section
    for ii in np.arange(row_start_nd, nz_glob):
        rcur    = rc_glob[ii,0]
        ff[ii], Nu_gas[ii], h_gas[ii] = gnielinski_correlation(Re_glob[ii], Pr, rcur, k_gas)

        Nu_gas2[ii], h_gas2[ii] = simple_corr(Re_glob[ii], Pr, rcur, k_gas)

        rch_cur = 0.5 * dh_ch[ii]   #use half of hydraulic diameter as is hard coded to use radius
        ff_ch[ii], Nu_ch[ii], h_ch[ii] = gnielinski_correlation(Re_ch[ii], Pr_fuel, rch_cur, k_fuel)

        #last node, end of chamber
        if ii == (nz_glob-1):
            # gas
            A_gas[ii] = area_h_slant_exit(r1 = rcur, rup = rc_glob[ii-1,0],
                                          dz = dz_nd)

            # Channel inner
            A_chi[ii] = area_h_slant_exit(r1 = rci[ii], rup = rci[ii-1],
                                          dz = dz_nd)

            # Channel outer
            A_cho[ii] = area_h_slant_exit(r1 = rco[ii], rup = rco[ii-1],
                                          dz = dz_nd)

            # shell outer
            A_sho[ii] = area_h_slant_exit(r1 = rso[ii], rup = rso[ii-1],
                                          dz = dz_nd)
        else:
            # gas
            A_gas[ii] = area_h_slant(rcur = rcur, rup = rc_glob[ii-1,0],
                                     rdown = rc_glob[ii+1,0], dz = dz_nd)

            # Channel inner
            A_chi[ii] = area_h_slant(rcur = rci[ii], rup = rci[ii-1],
                                     rdown = rci[ii+1], dz = dz_nd)

            # Channel outer
            A_cho[ii] = area_h_slant(rcur = rco[ii], rup = rco[ii-1],
                                     rdown = rco[ii+1], dz = dz_nd)

            # shell outer
            A_sho[ii] = area_h_slant(rcur = rso[ii], rup = rso[ii-1],
                                     rdown = rso[ii+1], dz = dz_nd)

        Rgas[ii]    = 1 / (h_gas[ii] * A_gas[ii])
        R_chi[ii]   = 1 / (h_ch[ii] * A_chi[ii])
        R_cho[ii]   = 1 / (h_ch[ii] * A_cho[ii])
        R_sho[ii]   = 1 / (h_inf * A_sho[ii])


    heat_out = {
        'R_gas'     : Rgas,
        'R_chi'     : R_chi,
        'R_cho'     : R_cho,
        'R_sho'     : R_sho,
        'ff_gas'    : ff,
        'Nu_gas'    : Nu_gas,
        'h_gas'     : h_gas,
        'ff_ch'     : ff_ch,
        'Nu_ch'     : Nu_ch,
        'h_ch'      : h_ch,
        'A_gas'     : A_gas,
        'A_chi'     : A_chi,
        'A_cho'     : A_cho,
        'A_sho'     : A_sho
        }
    ### compute resistances to radial heat flow in inner wall
    r_log  = np.zeros((nz_glob, nr-1))

    # ln(r2/r1)
    for jj in np.arange(nr-1):
        r_log[:,jj] = np.log(rc_glob[:,jj+1] / rc_glob[:,jj])

    dz_mat = np.expand_dims(dz_vec, axis=1)
    dz_mat = np.repeat(dz_mat, repeats=(nr-1), axis=1)

    #radial conduction resistances, inner wall
    Rr_w    = r_log / (2 * np.pi * dz_mat * k_w)

    ### compute resistance to radial heat flow in outer shell
    # ln(r2/r1)
    r_log_sh  = np.log(rso/rco)

    #radial conduction resistances
    Rr_sh    = r_log_sh / (2 * np.pi * dz_vec * k_s)

    heat_out['Rr_w']    = Rr_w
    heat_out['Rr_sh']   = Rr_sh

    ### Compute the equivalent resistance along each radial track, compute heat rate

    Rtot    = Rgas + R_cho + R_chi + R_sho + np.sum(Rr_w, axis=1) + Rr_sh
    dT_vec  = T_glob - np.repeat([T_inf], repeats = nz_glob)
    qr      = dT_vec / Rtot
    qr_flux     = qr / A_gas

    heat_out['R_tot']   = Rtot
    heat_out['qr']      = qr

    #thermal resistance between coolant and outer wall
    R_co = R_cho + Rr_sh + R_sho
    T_chan = np.zeros(nz_glob)
    T_chan[-1] = Tcool

    #succesively compute the coolant temperature, starting at the nozzle exit
    for ii in np.arange(nz_glob-2, -1, step=-1):
        num = qr[ii] + (mdot_cool * cp_fuel * T_chan[ii+1]) + (T_inf / R_co[ii])
        den = (mdot_cool * cp_fuel) + (1/R_co[ii])
        T_chan[ii] = num/den

    # thermal resistance between inner wall and coolant
    # compute the temperature at the inner wall
    R_wc    = np.sum(Rr_w, axis=1) + R_chi
    Twi     = qr * R_wc + T_chan

    #compute the temperature at all points in the inner wall
    Tw = np.zeros((nz_glob, nr))
    Tw[:,0] = Twi

    for ii in np.arange(1, nr,1):
        Tw[:,ii] = Tw[:,ii-1] - (qr * Rr_w[:,ii-1])

    # thermal resistance between coolant and Tinf
    qo = (T_chan - np.repeat([T_inf], repeats=nz_glob)) / R_co

    T_si = T_chan - (qo * R_cho)
    T_so = T_si - (qo * Rr_sh)

    Tw_glob         = np.zeros((nz_glob, (nr+3)))
    Tw_glob[:,:nr]  = Tw.copy()
    Tw_glob[:,nr]   = T_chan
    Tw_glob[:,nr+1] = T_si
    Tw_glob[:,-1]   = T_so

    heat_out['Tall_glob']   = Tw_glob
    heat_out['Tw']          = Tw
    heat_out['T_chan']      = T_chan
    heat_out['T_si']        = T_si
    heat_out['T_so']        = T_so
    heat_out['qr_flux_iw']  = qr_flux

    # fig = plt.figure()
    # plt.plot(zc_glob[:,0], qr_flux, 'ko')
    # plt.title('heat flux vs z')

    # fig = plt.figure()
    # plt.plot(zc_glob[:,0], Twi, 'ko')
    # plt.title('inner wall temp')

    # fig = plt.figure()
    # plt.contourf(zc_glob[:,:nr], rc_glob[:,:nr], Tw, levels=30, cmap=plt.cm.jet)
    # plt.colorbar()
    # plt.axis('equal')
    return heat_out

def comp_heat_orig(design_vars, flow_bc, material_prop, mesh_param, therm_out, geom_out, mesh_out, flow_out):

    ### Extract geometry info
    rc          = design_vars['rc']
    nr          = mesh_param['nr']
    r_inner     = mesh_out['zones']['inner']['rc']
    z_inner     = mesh_out['zones']['inner']['zc']
    Ac_inner    = mesh_out['zones']['inner']['Ac']
    beta_nc     = design_vars['beta_nc']
    beta_nd     = design_vars['beta_nd']
    dz_vec      = mesh_out['zones']['global']['dz_vec']

    ### Extract transport properties
    Pr          = therm_out['Pr']
    k_gas       = therm_out['k']
    Pr_fuel     = flow_bc['Pr_fuel']
    k_fuel      = flow_bc['k_fuel']
    T_inf       = flow_bc['T_inf']
    h_inf       = flow_bc['h_inf']

    #Extract solid conductivities
    k_w     = material_prop['k_w']
    k_s     = material_prop['k_s']

    ### Extract global mesh (reshaped into nz x nr)
    rc_glob     = mesh_out['zones']['global']['rc']
    zc_glob     = mesh_out['zones']['global']['zc']
    nodes_glob  = mesh_out['zones']['global']['nodes']
    nz_glob, nr_glob    = zc_glob.shape
    n_glob      = nz_glob * nr_glob

    row_start_nc    = mesh_out['zones']['global']['row_start_nc']
    row_t           = mesh_out['zones']['global']['row_t']
    row_start_nd    = mesh_out['zones']['global']['row_start_nd']

    dz_c        = zc_glob[1,0] - zc_glob[0,0]
    dz_nc       = zc_glob[row_start_nc+1,0] - zc_glob[row_start_nc,0]
    dz_nd       = zc_glob[row_start_nd+1,0] - zc_glob[row_start_nd,0]

    ### Extract global T, Re profiles (and channel)
    T_glob      = flow_out['global']['T']
    Re_glob     = flow_out['global']['Re']
    Re_ch       = flow_out['channel']['Re']
    dh_ch       = flow_out['channel']['hydraulic_d']
    Acc     = mesh_out['zones']['global']['Ac_channel']
    rco     = mesh_out['zones']['global']['rco']
    rci     = mesh_out['zones']['global']['rci']
    rso     = mesh_out['zones']['global']['rso']

    ### compute heat transfer coeffs and thermal resistances between gas flow
    # and inner wall along full engine length
    Rgas        = np.zeros(nz_glob)
    ff          = np.zeros(nz_glob)
    Nu_gas      = np.zeros(nz_glob)
    h_gas       = np.zeros(nz_glob)
    A_gas        = np.zeros(nz_glob)

    R_chi       = np.zeros(nz_glob)
    R_cho       = np.zeros(nz_glob)
    ff_ch       = np.zeros(nz_glob)
    Nu_ch       = np.zeros(nz_glob)
    h_ch        = np.zeros(nz_glob)
    A_chi        = np.zeros(nz_glob)
    A_cho        = np.zeros(nz_glob)

    A_sho       = np.zeros(nz_glob)
    R_sho       = np.zeros(nz_glob)

    Nu_gas2 = np.zeros(nz_glob)
    h_gas2 = np.zeros(nz_glob)

    # combustion chamber
    for ii in np.arange(row_start_nc):
        rcur    = rc_glob[ii,0]
        ff[ii], Nu_gas[ii], h_gas[ii] = gnielinski_correlation(Re_glob[ii], Pr, rcur, k_gas)

        Nu_gas2[ii], h_gas2[ii] = simple_corr(Re_glob[ii], Pr, rcur, k_gas)
        rch_cur = 0.5 * dh_ch[ii]   #use half of hydraulic diameter as is hard coded to use radius
        ff_ch[ii], Nu_ch[ii], h_ch[ii] = gnielinski_correlation(Re_ch[ii], Pr_fuel, rch_cur, k_fuel)

        #first node, only half dz_c in height
        if ii == 0:
            A_gas[ii]   = np.pi * rcur * dz_c
            A_chi[ii]   = np.pi * rci[ii] * dz_c
            A_cho[ii]   = np.pi * rco[ii] * dz_c
            A_sho[ii]   = np.pi * rso[ii] * dz_c
        elif ii == (row_start_nc-1):
            # gas
            A_gas[ii] = area_h_cyl_conv(r1 = rcur, rdown = rc_glob[ii+1,0],
                                     dz_straight = dz_c, dz_slant = dz_nc)

            # channel inner
            A_chi[ii] = area_h_cyl_conv(r1 = rci[ii], rdown = rci[ii+1],
                                     dz_straight = dz_c, dz_slant = dz_nc)

            # channel outer
            A_cho[ii] = area_h_cyl_conv(r1 = rco[ii], rdown = rco[ii+1],
                                     dz_straight = dz_c, dz_slant = dz_nc)

            # shell outer
            A_sho[ii] = area_h_cyl_conv(r1 = rso[ii], rdown = rso[ii+1],
                                     dz_straight = dz_c, dz_slant = dz_nc)
        else:
            A_gas[ii]   = 2 * np.pi * rcur * dz_c
            A_chi[ii]   = 2 * np.pi * rci[ii] * dz_c
            A_cho[ii]   = 2 * np.pi * rco[ii] * dz_c
            A_sho[ii]   = 2 * np.pi * rso[ii] * dz_c

        Rgas[ii]    = 1 / (h_gas[ii] * A_gas[ii])
        R_chi[ii]   = 1 / (h_ch[ii] * A_chi[ii])
        R_cho[ii]   = 1 / (h_ch[ii] * A_cho[ii])
        R_sho[ii]   = 1 / (h_inf * A_sho[ii])

    #nozzle converging section
    for ii in np.arange(row_start_nc, row_start_nd, 1):
        rcur    = rc_glob[ii,0]
        ff[ii], Nu_gas[ii], h_gas[ii] = gnielinski_correlation(Re_glob[ii], Pr, rcur, k_gas)

        Nu_gas2[ii], h_gas2[ii] = simple_corr(Re_glob[ii], Pr, rcur, k_gas)

        rch_cur = 0.5 * dh_ch[ii]   #use half of hydraulic diameter as is hard coded to use radius
        ff_ch[ii], Nu_ch[ii], h_ch[ii] = gnielinski_correlation(Re_ch[ii], Pr_fuel, rch_cur, k_fuel)

        #throat
        if ii == (row_start_nd-1):

            # gas side
            A_gas[ii] = area_h_throat(r1 = rcur, rup = rc_glob[ii-1,0], dz_up = dz_nc,
                                    rdown = rc_glob[ii+1,0], dz_down = dz_nd)

            #channel inner
            A_chi[ii] = area_h_throat(r1 = rci[ii], rup = rci[ii-1], dz_up = dz_nc,
                                    rdown = rci[ii+1], dz_down = dz_nd)

            #channel outer
            A_cho[ii] = area_h_throat(r1 = rco[ii], rup = rco[ii-1], dz_up = dz_nc,
                                    rdown = rco[ii+1], dz_down = dz_nd)

            # shell outer
            A_sho[ii] = area_h_throat(r1 = rso[ii], rup = rso[ii-1], dz_up = dz_nc,
                                    rdown = rso[ii+1], dz_down = dz_nd)
        else:
            # gas
            A_gas[ii] = area_h_slant(rcur = rcur, rup = rc_glob[ii-1,0],
                                     rdown = rc_glob[ii+1,0], dz = dz_nc)

            # Channel inner
            A_chi[ii] = area_h_slant(rcur = rci[ii], rup = rci[ii-1],
                                     rdown = rci[ii+1], dz = dz_nc)

            # Channel outer
            A_cho[ii] = area_h_slant(rcur = rco[ii], rup = rco[ii-1],
                                     rdown = rco[ii+1], dz = dz_nc)

            # shell outer
            A_sho[ii] = area_h_slant(rcur = rso[ii], rup = rso[ii-1],
                                     rdown = rso[ii+1], dz = dz_nc)

        Rgas[ii]    = 1 / (h_gas[ii] * A_gas[ii])
        R_chi[ii]   = 1 / (h_ch[ii] * A_chi[ii])
        R_cho[ii]   = 1 / (h_ch[ii] * A_cho[ii])
        R_sho[ii]   = 1 / (h_inf * A_sho[ii])

    #nozzle diverging section
    for ii in np.arange(row_start_nd, nz_glob):
        rcur    = rc_glob[ii,0]
        ff[ii], Nu_gas[ii], h_gas[ii] = gnielinski_correlation(Re_glob[ii], Pr, rcur, k_gas)

        Nu_gas2[ii], h_gas2[ii] = simple_corr(Re_glob[ii], Pr, rcur, k_gas)

        rch_cur = 0.5 * dh_ch[ii]   #use half of hydraulic diameter as is hard coded to use radius
        ff_ch[ii], Nu_ch[ii], h_ch[ii] = gnielinski_correlation(Re_ch[ii], Pr_fuel, rch_cur, k_fuel)

        #last node, end of chamber
        if ii == (nz_glob-1):
            # gas
            A_gas[ii] = area_h_slant_exit(r1 = rcur, rup = rc_glob[ii-1,0],
                                          dz = dz_nd)

            # Channel inner
            A_chi[ii] = area_h_slant_exit(r1 = rci[ii], rup = rci[ii-1],
                                          dz = dz_nd)

            # Channel outer
            A_cho[ii] = area_h_slant_exit(r1 = rco[ii], rup = rco[ii-1],
                                          dz = dz_nd)

            # shell outer
            A_sho[ii] = area_h_slant_exit(r1 = rso[ii], rup = rso[ii-1],
                                          dz = dz_nd)
        else:
            # gas
            A_gas[ii] = area_h_slant(rcur = rcur, rup = rc_glob[ii-1,0],
                                     rdown = rc_glob[ii+1,0], dz = dz_nd)

            # Channel inner
            A_chi[ii] = area_h_slant(rcur = rci[ii], rup = rci[ii-1],
                                     rdown = rci[ii+1], dz = dz_nd)

            # Channel outer
            A_cho[ii] = area_h_slant(rcur = rco[ii], rup = rco[ii-1],
                                     rdown = rco[ii+1], dz = dz_nd)

            # shell outer
            A_sho[ii] = area_h_slant(rcur = rso[ii], rup = rso[ii-1],
                                     rdown = rso[ii+1], dz = dz_nd)

        Rgas[ii]    = 1 / (h_gas[ii] * A_gas[ii])
        R_chi[ii]   = 1 / (h_ch[ii] * A_chi[ii])
        R_cho[ii]   = 1 / (h_ch[ii] * A_cho[ii])
        R_sho[ii]   = 1 / (h_inf * A_sho[ii])


    heat_out = {
        'R_gas'     : Rgas,
        'R_chi'     : R_chi,
        'R_cho'     : R_cho,
        'R_sho'     : R_sho,
        'ff_gas'    : ff,
        'Nu_gas'    : Nu_gas,
        'h_gas'     : h_gas,
        'ff_ch'     : ff_ch,
        'Nu_ch'     : Nu_ch,
        'h_ch'      : h_ch,
        'A_gas'     : A_gas,
        'A_chi'     : A_chi,
        'A_cho'     : A_cho,
        'A_sho'     : A_sho
        }
    ### compute resistances to radial heat flow in inner wall
    r_log  = np.zeros((nz_glob, nr-1))

    # ln(r2/r1)
    for jj in np.arange(nr-1):
        r_log[:,jj] = np.log(rc_glob[:,jj+1] / rc_glob[:,jj])

    dz_mat = np.expand_dims(dz_vec, axis=1)
    dz_mat = np.repeat(dz_mat, repeats=(nr-1), axis=1)

    #radial conduction resistances, inner wall
    Rr_w    = r_log / (2 * np.pi * dz_mat * k_w)

    ### compute resistance to radial heat flow in outer shell
    # ln(r2/r1)
    r_log_sh  = np.log(rso/rco)

    #radial conduction resistances
    Rr_sh    = r_log_sh / (2 * np.pi * dz_vec * k_s)

    heat_out['Rr_w']    = Rr_w
    heat_out['Rr_sh']   = Rr_sh

    ### Compute the equivalent resistance along each radial track, compute heat rate

    Rgas    = Rgas*1000
    Rtot    = Rgas + R_cho + R_chi + R_sho + np.sum(Rr_w, axis=1) + Rr_sh
    dT_vec  = T_glob - np.repeat([T_inf], repeats = nz_glob)
    qr      = dT_vec / Rtot
    qr_flux     = qr / A_gas

    heat_out['R_tot']   = Rtot
    heat_out['qr']      = qr

    Twi     = T_glob - (qr * Rgas)
    Tw = np.zeros((nz_glob, nr))
    Tw[:,0] = Twi

    for ii in np.arange(1, nr,1):
        Tw[:,ii] = Tw[:,ii-1] - (qr * Rr_w[:,ii-1])

    #temperature on channel CL
    R_chan = Rgas + np.sum(Rr_w, axis=1) + R_chi
    T_chan = T_glob - (qr * R_chan)

    R_si = Rgas + np.sum(Rr_w, axis=1) + R_chi + R_cho
    T_si = T_glob - (qr * R_si)

    R_so = Rgas + np.sum(Rr_w, axis=1) + R_chi + R_cho + Rr_sh
    T_so = T_glob - (qr * R_so )

    Tw_glob         = np.zeros((nz_glob, (nr+3)))
    Tw_glob[:,:nr]  = Tw.copy()
    Tw_glob[:,nr]   = T_chan
    Tw_glob[:,nr+1] = T_si
    Tw_glob[:,-1]   = T_so

    heat_out['Tall_glob']   = Tw_glob
    heat_out['Tw']          = Tw
    heat_out['T_chan']      = T_chan
    heat_out['T_si']        = T_si
    heat_out['T_so']        = T_so

    fig = plt.figure()
    plt.plot(zc_glob[:,0], qr_flux, 'ko')
    plt.title('heat flux vs z')

    fig = plt.figure()
    plt.plot(zc_glob[:,0], Twi, 'ko')
    plt.title('inner wall temp')

    fig = plt.figure()
    plt.contourf(zc_glob[:,:nr], rc_glob[:,:nr], Tw, levels=30, cmap=plt.cm.jet)
    plt.colorbar()
    plt.axis('equal')
    return heat_out

def gnielinski_correlation(Re, Pr, r, k):
    '''Compute the friction factor, Nusselt number, and
    convective heat transfer coefficient using Gnielinski correlation'''

    f       = (0.79 * np.log(Re) - 1.64)**(-2)
    num     = (f/8) * (Re - 1000) * Pr
    den     = 1 + (12.7 * (f/8)**(0.5) * ((Pr**(2/3))-1))
    Nu      = num / den
    h       = (Nu * k) / (2 * r)
    return f, Nu, h

def simple_corr(Re, Pr, r, k):
    Nu = 0.026 * Re**0.8 * Pr**0.4
    h       = (Nu * k) / (2 * r)
    return Nu, h


def area_h_cyl_conv(r1, rdown, dz_straight, dz_slant):
    '''Compute surface area for node at interface between
    combustion chamber and converging nozzle section'''
    r2      = 0.5 * (r1 + rdown)
    h       = dz_slant / 2
    s       = np.sqrt((r1-r2)**2 + h**2)
    Aslant  = np.pi * (r1 + r2) * s
    Ausual  = np.pi * r1 * dz_straight
    A       = Aslant + Ausual

    return A

def area_h_throat(r1, rup, dz_up, rdown, dz_down):
    '''Compute surface area for node at throat'''

    r2          = 0.5 * (r1 + rup)
    h           = 0.5 * dz_up
    s           = np.sqrt((r1-r2)**2 + h**2)
    Aslant1     = np.pi * (r1 + r2) * s

    r2          = 0.5 * (r1 + rdown)
    h           = 0.5 * dz_down
    s           = np.sqrt((r1-r2)**2 + h**2)
    Aslant2     = np.pi * (r1 + r2) * s
    A = Aslant1 + Aslant2
    return A

def area_h_slant(rcur, rup, rdown, dz):
    '''Compute surface area for node on cone frustum'''

    r1      = 0.5 * (rcur + rup)
    r2      = 0.5 * (rcur + rdown)
    h       = dz
    s       = np.sqrt((r1-r2)**2 + h**2)
    Aslant  = np.pi * (r1 + r2) * s
    return Aslant

def area_h_slant_exit(r1, rup, dz):
    '''Compute surface area at nozzle exit node'''

    r2      = 0.5 * (r1 + rup)
    h       = 0.5 * dz
    s       = np.sqrt((r1-r2)**2 + h**2)
    Aslant  = np.pi * (r1 + r2) * s
    return Aslant





