#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: jd
contact: jamesduv@umich.edu
affiliation: University of Michigan, Department of Aerospace Eng., CASLAB
"""

import numpy as np
import matplotlib.pyplot as plt

def comp_mesh(design_vars, geom_out, mesh_param):
    '''Generate the meshes for the models'''

    nr      = mesh_param['nr']
    nz_c    = mesh_param['nz_c']
    nz_nc   = mesh_param['nz_nc']
    nz_nd   = mesh_param['nz_nd']
    nr_global   = nr + 3    #wall mesh + cooling channel temp + 2 shell nodes
    tc      = design_vars['tc']
    ts      = design_vars['ts']

    Lc      = design_vars['Lc']
    rc      = design_vars['rc']
    Lnc     = geom_out['Lnc']
    Lnd     = geom_out['Lnd']
    rt      = geom_out['rt']
    re      = geom_out['re']
    tw      = design_vars['tw']
    beta_nc = design_vars['beta_nc']
    beta_nd = design_vars['beta_nd']
    rco     = rc + tw

    #compute the z coordinates of the interfaces
    #combustion chamber start and stop
    zc_start    = 0
    zc_end      = Lc

    #converging section start and stop
    znc_start   = Lc
    znc_end     = Lc + Lnc

    #throat
    zt          = znc_end

    #diverging section
    znd_start   = zt
    znd_end     = zt + Lnd

    mesh_info   = {
        'zc_start'  : zc_start,
        'zc_end'    : zc_end,
        'znc_start' : znc_start,
        'znc_end'   : znc_end,
        'zt'        : zt,
        'znd_start' : znd_start,
        'znd_end'   : znd_end
        }

    mesh_zones = {}

    ### Combustion chamber
    zvec_c  = np.linspace(0, Lc, nz_c )
    rvec_c  = np.linspace(rc, rco, nr)
    dr      = rvec_c[1] - rvec_c[0]

    r_inner     = np.repeat(np.array([rc]), repeats = nz_c)

    ZC_C, RC_C  = np.meshgrid(zvec_c, rvec_c)

    zc_c = np.reshape(ZC_C, (nz_c*nr,1), order='C')
    rc_c = np.reshape(RC_C, (nz_c*nr,1), order='C')

    #node numbering for combustion chamber wall mesh only
    nodes_c  = np.arange(nz_c*nr)

    #global node numbering, including cooling channel, 2 nodes on outer shell
    nodes_c_glob  = np.arange(nz_c*nr_global)

    mesh_zones['c'] = {}
    mesh_zones['c']['zc'] = zc_c
    mesh_zones['c']['rc'] = rc_c
    mesh_zones['c']['nodes'] = nodes_c

    # fig = plt.figure()
    # plt.plot(rc_c, zc_c, 'rx')

    #track the z coordinates of the inner wall - leave repeats
    z_inner  = zvec_c.copy()

    ### Converging nozzle section
    zvec_nc     = np.linspace(0, Lnc, nz_nc)
    rvec_nc     = rc - zvec_nc * np.tan(beta_nc * np.pi / 180)

    r_inner     = np.concatenate((r_inner, rvec_nc.copy()))

    #make z-vector global
    zvec_nc     = zvec_nc + znc_start
    zvec_nc1    = np.expand_dims(zvec_nc, axis=1)
    zvec_full   = np.repeat(zvec_nc1, repeats=nr, axis=1)
    zvec_full   = np.reshape(zvec_full, (nz_nc*nr,), order='C')

    z_inner     = np.concatenate((z_inner, zvec_nc), axis=0)

    rc_old  = rvec_nc.copy()
    for ii in np.arange(nr-1):
        rc_new = rc_old + dr
        rvec_nc = np.concatenate((rvec_nc, rc_new.copy()), axis=0)
        rc_old = rc_new.copy()
    # plt.plot(rvec_nc, zvec_full, 'ko', alpha=0.5)

    zc_nc   = zvec_full.copy()
    rc_nc   = rvec_nc.copy()

    #node numbering, converging section only
    nodes_nc    = np.arange(nz_nc*nr)

    #global node numbering, including cooling channel, 2 nodes on outer shell
    nodes_nc_glob  = np.arange(nz_nc*nr_global)

    mesh_zones['nc'] = {}
    mesh_zones['nc']['zc'] = zc_nc
    mesh_zones['nc']['rc'] = rc_nc
    mesh_zones['nc']['nodes'] = nodes_nc

    ### diverging nozzle section
    zvec_nd     = np.linspace(0, Lnd, nz_nd)
    rvec_nd     = rt + zvec_nd * np.tan(beta_nd * np.pi / 180)

    r_inner     = np.concatenate((r_inner, rvec_nd), axis=0)
    # make z vector global
    zvec_nd     = znd_start + zvec_nd

    zvec_nd1    = np.expand_dims(zvec_nd, axis=1)
    zvec_full   = np.repeat(zvec_nd1, repeats=nr, axis=1)
    zvec_full   = np.reshape(zvec_full, (nz_nd*nr,), order='C')

    z_inner     = np.concatenate((z_inner, zvec_nd), axis=0)

    rc_old  = rvec_nd.copy()
    for ii in np.arange(nr-1):
        rc_new = rc_old + dr
        rvec_nd = np.concatenate((rvec_nd, rc_new.copy()), axis=0)
        rc_old = rc_new.copy()
    # plt.plot(rvec_nd, zvec_full, 'bx')

    zc_nd       = zvec_full.copy()
    rc_nd       = rvec_nd.copy()
    nodes_nd    = np.arange(nz_nd*nr)

    mesh_zones['nd'] = {}
    mesh_zones['nd']['zc'] = zc_nd
    mesh_zones['nd']['rc'] = rc_nd
    mesh_zones['nd']['nodes'] = nodes_nd

    # plt.plot( [rt, rt], [0, znd_end], '--k')
    # plt.plot( [re, re], [0, znd_end], '--b')
    # plt.axis('equal')

    mesh_zones['inner'] = {}
    mesh_zones['inner']['zc'] = z_inner
    mesh_zones['inner']['rc'] = r_inner
    mesh_zones['inner']['Ac'] = np.pi * r_inner**2


    ### build global mesh
    rc1 = np.reshape(rc_c, (nz_c, nr), order='F')
    zc1 = np.reshape(zc_c, (nz_c, nr), order='F')

    rc2 = np.reshape(rc_nc, (nz_nc, nr), order='F')
    zc2 = np.reshape(zc_nc, (nz_nc, nr), order='C')

    rc3 = np.reshape(rc_nd, (nz_nd, nr), order='F')
    zc3 = np.reshape(zc_nd, (nz_nd, nr), order='C')

    rc_glob = np.concatenate((rc1, rc2[1:,:], rc3[1:,:]), axis=0)
    zc_glob = np.concatenate((zc1, zc2[1:,:], zc3[1:,:]), axis=0)

    nz_glob     = zc_glob.shape[0]
    nr_glob     = nr + 3

    # add point for cooling channel centerline
    tc_rep  = np.repeat([0.5*tc], repeats = nz_glob )
    cvec    = rc_glob[:,-1] + tc_rep
    cvec    = np.expand_dims(cvec, axis=1)
    rc_glob     = np.concatenate((rc_glob, cvec), axis=1)

    #add point for shell inner surface
    cvec        = rc_glob[:,-1] + tc_rep
    cvec        = np.expand_dims(cvec, axis=1)
    rc_glob     = np.concatenate((rc_glob, cvec), axis=1)

    #add point for shell outer surface
    ts_rep  = np.repeat([ts], repeats = nz_glob )
    cvec    = rc_glob[:,-1] + ts_rep
    cvec    = np.expand_dims(cvec, axis=1)
    rc_glob     = np.concatenate((rc_glob, cvec), axis=1)

    zvec_glob   = zc_glob[:,[0]]
    zc_glob     = np.concatenate((zc_glob, zvec_glob, zvec_glob, zvec_glob), axis=1)

    #reshape global coordinates, use C (row wise) reshaping
    nglob   = nz_glob*nr_glob
    nodes_glob  = np.arange(nglob)
    nodes_glob  = np.reshape(nodes_glob, (nz_glob, nr_glob), order='C')
    # rc_glob = np.reshape(rc_glob, (nglob,), order='C')
    # zc_glob = np.reshape(zc_glob, (nglob,), order='C')

    mesh_zones['global'] = {}
    mesh_zones['global']['rc'] = rc_glob
    mesh_zones['global']['zc'] = zc_glob
    mesh_zones['global']['nodes'] = nodes_glob
    nz_chamber      = nz_c
    row_start_nc    = nz_c
    row_t           = row_start_nc + nz_nc - 2
    row_start_nd    = row_t + 1

    zc_check = np.reshape(zc_glob, (nz_glob, nr_glob), order='C')
    rc_check    = np.reshape(rc_glob, (nz_glob, nr_glob), order='C')

    mesh_zones['global']['row_start_nc']    = row_start_nc
    mesh_zones['global']['row_t']           = row_t
    mesh_zones['global']['row_start_nd']    = row_start_nd

    #build vector with local dz
    dz_c        = zc_glob[1,0] - zc_glob[0,0]
    dz_nc       = zc_glob[row_start_nc+1,0] - zc_glob[row_start_nc,0]
    dz_nd       = zc_glob[row_start_nd+1,0] - zc_glob[row_start_nd,0]
    dz_vec     = np.zeros(nz_glob)
    dz_vec[:row_start_nc] = np.repeat([dz_c], repeats = (row_start_nc))
    dz_vec[0] = 0.5 * dz_c
    dz_vec[row_start_nc-1] = 0.5 * (dz_c + dz_nc)

    nrep    = row_start_nd - row_start_nc
    dz_vec[row_start_nc:row_start_nd] = np.repeat([dz_nc], repeats=nrep)
    dz_vec[row_t] = 0.5 * (dz_nc + dz_nd)

    nrep = nz_glob - row_start_nd
    dz_vec[row_start_nd:] = np.repeat([dz_nd], repeats=nrep)
    dz_vec[-1] = 0.5 * dz_nd

    mesh_zones['global']['dz_vec'] = dz_vec
    mesh_info['dz_c']   = dz_c
    mesh_info['dz_nc']  = dz_nc
    mesh_info['dz_nd']  = dz_nd

    # build global beta vector
    beta_glob = np.zeros(nz_glob)
    nrep    = row_start_nd - row_start_nc
    beta_glob[row_start_nc:row_start_nd] = np.repeat([beta_nc], repeats=nrep)

    beta_glob[row_t] = 0

    nrep = nz_glob - row_start_nd
    beta_glob[row_start_nd:] = np.repeat([beta_nd], repeats=nrep)
    mesh_zones['global']['beta_vec'] = beta_glob

    # Compute the cross sectional area vs z of the cooling channel
    rci         = rc_glob[:,-1]  #radius, channel, inner
    rco         = rci + np.repeat([tc], repeats = nz_glob)
    rso         = rco + np.repeat([ts], repeats = nz_glob)
    Ac_channel  = np.pi * (rco**2 - rci**2)

    mesh_zones['global']['Ac_channel']  = Ac_channel
    mesh_zones['global']['rci']         = rci
    mesh_zones['global']['rco']         = rco
    mesh_zones['global']['rso']         = rso

    mesh_out = {
        'info'     : mesh_info,
        'zones'    : mesh_zones
        }

    return mesh_out









