#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: jd
contact: jamesduv@umich.edu
affiliation: University of Michigan, Department of Aerospace Eng., CASLAB
"""
import os, pickle
import numpy as np
import matplotlib.pyplot as plt

def load_test():

    read_dir    = 'model_output_norm'
    fn       = os.path.join(read_dir, 'isp325_mass.pickle')
    d200    = pickle.load(open(fn, 'rb'))
    return d200

def compare_all_constr():
    plot_opt = get_plot_opt()
    read_dir    = 'model_output_norm'


    fn1     = os.path.join(read_dir, 'isp325_thrust3.000e+05_stress_mass.pickle')
    fn2     = os.path.join(read_dir, 'isp325_thrust3.000e+05_stress_Temp_mass.pickle')

    d200      = pickle.load(open(fn1, 'rb'))
    d250      = pickle.load(open(fn2, 'rb'))
    # d300      = pickle.load(open(fn3, 'rb'))
    # d400      = pickle.load(open(fn4, 'rb'))

    mesh200 = d200['model_out']['mesh_out']
    mesh250 = d250['model_out']['mesh_out']

    mesh_param = d250['model_out']['mesh_param']
    nr      = mesh_param['nr']
    # mesh300 = d300['model_out']['mesh_out']
    # mesh400 = d400['model_out']['mesh_out']

    zc200   = mesh200['zones']['global']['zc']
    zc250   = mesh250['zones']['global']['zc']
    # zc300   = mesh300['zones']['global']['zc']
    # zc400   = mesh400['zones']['global']['zc']

    rc200   = mesh200['zones']['global']['rc']
    rc250   = mesh250['zones']['global']['rc']

    heat200 = d200['model_out']['heat_out']
    heat250 = d250['model_out']['heat_out']
    # heat300 = d300['model_out']['heat_out']
    # heat400 = d400['model_out']['heat_out']

    stress200 = d200['model_out']['stress_out']
    stress250 = d250['model_out']['stress_out']
    # stress300 = d300['model_out']['stress_out']
    # stress400 = d400['model_out']['stress_out']
    material_prop   = d200['model_out']['material_prop']
    sy  = material_prop['s_yield']

    s200    = abs( stress200['s_tot'])
    s250    = abs( stress250['s_tot'])
    # s300    = abs( stress300['s_tot'])
    # s400    = abs( stress400['s_tot'])

    st200    = abs( stress200['s_cte'])
    st250    = abs( stress250['s_cte'])

    sp200    = abs( stress200['s_hoop'])
    sp250    = abs( stress250['s_hoop'])

    Tw200   = heat200['Tw']
    Tw250   = heat250['Tw']
    # Tw300   = heat300['Tw']
    # Tw400   = heat400['Tw']

    Tc200   = heat200['T_chan']
    Tc250   = heat250['T_chan']
    # Tc300   = heat300['T_chan']
    # Tc400   = heat400['T_chan']

    qf200   = heat200['qr_flux_iw']
    qf250   = heat250['qr_flux_iw']
    # qf300   = heat300['qr_flux_iw']
    # qf400   = heat400['qr_flux_iw']

    flow200 = d200['model_out']['flow_out']
    flow250 = d250['model_out']['flow_out']
    # flow300 = d300['model_out']['flow_out']
    # flow400 = d400['model_out']['flow_out']

    Tg200   = flow200['global']['T']
    Tg250   = flow250['global']['T']
    # Tg300   = flow300['global']['T']
    # Tg400   = flow400['global']['T']

    M200   = flow200['global']['M']
    M250   = flow250['global']['M']
    # M300   = flow300['global']['M']
    # M400   = flow400['global']['M']

    u200   = flow200['global']['u']
    u250   = flow250['global']['u']
    # u300   = flow300['global']['u']
    # u400   = flow400['global']['u']

    p200   = flow200['global']['p']
    p250   = flow250['global']['p']
    # p300   = flow300['global']['p']
    # p400   = flow400['global']['p']
    fig = plt.figure(figsize=(20,10))
    ax1 = fig.add_subplot(221)
    plt.contourf(zc200[:,:nr], rc200[:,:nr], Tw200, levels=30)
    ax1.axis('equal')
    plt.grid()
    plt.xlabel('$z$ [m]', **plot_opt['axis_font'])
    plt.ylabel('$r$ [m]', **plot_opt['axis_font'])
    plt.colorbar()

    ax2 = fig.add_subplot(222)
    plt.contourf(zc250[:,:nr], rc250[:,:nr], Tw250, levels=30)
    ax2.axis('equal')
    plt.grid()
    plt.xlabel('$z$ [m]', **plot_opt['axis_font'])
    plt.ylabel('$r$ [m]', **plot_opt['axis_font'])
    plt.colorbar()

    ax3 = fig.add_subplot(223, sharex=ax1)
    plt.plot(zc200[:,0], st200, '-r', label='s therm.')
    plt.plot(zc200[:,0], sp200, '-b', label='s press.')
    plt.plot(zc200[:,0], s200, '-k', label='s tot')
    # plt.plot(zc250[:,0], s250, '-r', label=r'With stress cons.')
    # plt.plot(zc300[:,0], s300, '-g', label='340s')
    # plt.plot(zc400[:,0], s400, '-b', label='350s')
    xl = plt.xlim()
    plt.plot(xl, [sy, sy], '--k')
    plt.grid()
    plt.xlabel('$z$ [m]', **plot_opt['axis_font'])
    plt.ylabel('$\sigma$ [Pa] ', **plot_opt['axis_font'])
    plt.legend(prop=plot_opt['cb_font'])
    # fig.tight_layout()

    ax4 = fig.add_subplot(224, sharex=ax2)
    plt.plot(zc250[:,0], st250, '-r', label='s therm.')
    plt.plot(zc250[:,0], sp250, '-b', label='s press.')
    plt.plot(zc250[:,0], s250, '-k', label='s tot')
    # plt.plot(zc250[:,0], s250, '-r', label=r'With stress cons.')
    # plt.plot(zc300[:,0], s300, '-g', label='340s')
    # plt.plot(zc400[:,0], s400, '-b', label='350s')
    xl = plt.xlim()
    plt.plot(xl, [sy, sy], '--k')
    plt.grid()
    plt.xlabel('$z$ [m]', **plot_opt['axis_font'])
    plt.ylabel('$\sigma$ [Pa] ', **plot_opt['axis_font'])
    plt.legend(prop=plot_opt['cb_font'])
    fig.tight_layout()

    fn_save = os.path.join('tex',r'thrust_isp_stress_temp_constr_contourf_compare.png')
    plt.savefig(fn_save, **plot_opt['save'])
    plt.close()

    return

    fig = plt.figure(figsize=(20,12))
    ax1 = fig.add_subplot(231)
    plt.plot(zc200[:,0], Tg200, '-k', label=r'No temp. cons.')
    plt.plot(zc250[:,0], Tg250, '-r', label=r'With temp. cons.')
    # plt.plot(zc300[:,0], Tg300, '-g', label='340s')
    # plt.plot(zc400[:,0], Tg400, '-b', label='350s')
    plt.grid()
    plt.xlabel('$z$ [m]', **plot_opt['axis_font'])
    plt.ylabel('$T_{gas}$ [K]', **plot_opt['axis_font'])
    plt.legend(prop=plot_opt['cb_font'])

    ax2 = fig.add_subplot(233)
    plt.plot(zc200[:,0], qf200, '-k', label=r'No temp. cons.')
    plt.plot(zc250[:,0], qf250, '-r', label=r'With temp. cons.')
    # plt.plot(zc300[:,0], qf300, '-g', label='340s')
    # plt.plot(zc400[:,0], qf400, '-b', label='350s')
    plt.grid()
    plt.xlabel('$z$ [m]', **plot_opt['axis_font'])
    plt.ylabel('$q`` [W/m^{2}]$  ', **plot_opt['axis_font'])
    plt.legend(prop=plot_opt['cb_font'])


    ax3 = fig.add_subplot(232)
    plt.plot(zc200[:,0], u200, '-k', label=r'No temp. cons.')
    plt.plot(zc250[:,0], u250, '-r', label=r'With temp. cons.')
    # plt.plot(zc300[:,0], u300, '-g', label='340s')
    # plt.plot(zc400[:,0], u400, '-b', label='350s')
    plt.grid()
    plt.xlabel('$z$ [m]', **plot_opt['axis_font'])
    plt.ylabel('$u$ [m/s] ', **plot_opt['axis_font'])
    plt.legend(prop=plot_opt['cb_font'])
    sy  = material_prop['s_yield']

    ax4 = fig.add_subplot(234)
    plt.plot(zc200[:,0], Tw200[:,0], '-k', label=r'No temp. cons.')
    plt.plot(zc250[:,0], Tw250[:,0], '-r', label=r'With temp cons.')
    # plt.plot(zc300[:,0], Tw300[:,0], '-g', label='340s')
    # plt.plot(zc400[:,0], Tw400[:,0], '-b', label='350s')
    xlim = plt.xlim()
    plt.plot(xlim, [700, 700], '--k')
    plt.grid()
    plt.xlabel('$z$ [m]', **plot_opt['axis_font'])
    plt.ylabel('$T_{w,i}$ [K]', **plot_opt['axis_font'])
    plt.legend(prop=plot_opt['cb_font'])

    ax5 = fig.add_subplot(235)
    plt.plot(zc200[:,0], Tc200, '-k', label=r'No temp. cons.')
    plt.plot(zc250[:,0], Tc250, '-r', label=r'With temp. cons.')
    # plt.plot(zc300[:,0], Tc300, '-g', label='340s')
    # plt.plot(zc400[:,0], Tc400, '-b', label='350s')
    plt.grid()
    plt.xlabel('$z$ [m]', **plot_opt['axis_font'])
    plt.ylabel('$T_{c}$ [K]', **plot_opt['axis_font'])
    plt.legend(prop=plot_opt['cb_font'])

    ax6 = fig.add_subplot(236)
    plt.plot(zc200[:,0], s200, '-k', label=r'No temp. cons.')
    plt.plot(zc250[:,0], s250, '-r', label=r'With tmep. cons.')
    # plt.plot(zc300[:,0], s300, '-g', label='340s')
    # plt.plot(zc400[:,0], s400, '-b', label='350s')
    xl = plt.xlim()
    plt.plot(xl, [sy, sy], '--k')
    plt.grid()
    plt.xlabel('$z$ [m]', **plot_opt['axis_font'])
    plt.ylabel('$\sigma$ [Pa] ', **plot_opt['axis_font'])
    plt.legend(prop=plot_opt['cb_font'])
    fig.tight_layout()

    fn_save = os.path.join('tex',r'thrust_isp_stress_temp_constr_compare.png')
    plt.savefig(fn_save, **plot_opt['save'])
    plt.close()

def compare_isp_thrust_stress():
    plot_opt = get_plot_opt()
    read_dir    = 'model_output_norm'

    fn1     = os.path.join(read_dir, 'isp325_thrust3.000e+05_mass.pickle')
    fn2     = os.path.join(read_dir, 'isp325_thrust3.000e+05_stress_mass.pickle')
    # fn3     = os.path.join(read_dir, 'isp340_thrust3.000e+05_mass.pickle')
    # fn4     = os.path.join(read_dir, 'isp350_thrust3.000e+05_mass.pickle')

    d200      = pickle.load(open(fn1, 'rb'))
    d250      = pickle.load(open(fn2, 'rb'))
    # d300      = pickle.load(open(fn3, 'rb'))
    # d400      = pickle.load(open(fn4, 'rb'))

    mesh200 = d200['model_out']['mesh_out']
    mesh250 = d250['model_out']['mesh_out']

    mesh_param = d250['model_out']['mesh_param']
    nr      = mesh_param['nr']
    # mesh300 = d300['model_out']['mesh_out']
    # mesh400 = d400['model_out']['mesh_out']

    zc200   = mesh200['zones']['global']['zc']
    zc250   = mesh250['zones']['global']['zc']
    # zc300   = mesh300['zones']['global']['zc']
    # zc400   = mesh400['zones']['global']['zc']

    rc200   = mesh200['zones']['global']['rc']
    rc250   = mesh250['zones']['global']['rc']

    heat200 = d200['model_out']['heat_out']
    heat250 = d250['model_out']['heat_out']
    # heat300 = d300['model_out']['heat_out']
    # heat400 = d400['model_out']['heat_out']

    stress200 = d200['model_out']['stress_out']
    stress250 = d250['model_out']['stress_out']
    # stress300 = d300['model_out']['stress_out']
    # stress400 = d400['model_out']['stress_out']
    material_prop   = d200['model_out']['material_prop']
    sy  = material_prop['s_yield']

    s200    = abs( stress200['s_tot'])
    s250    = abs( stress250['s_tot'])
    # s300    = abs( stress300['s_tot'])
    # s400    = abs( stress400['s_tot'])

    st200    = abs( stress200['s_cte'])
    st250    = abs( stress250['s_cte'])

    sp200    = abs( stress200['s_hoop'])
    sp250    = abs( stress250['s_hoop'])

    Tw200   = heat200['Tw']
    Tw250   = heat250['Tw']
    # Tw300   = heat300['Tw']
    # Tw400   = heat400['Tw']

    Tc200   = heat200['T_chan']
    Tc250   = heat250['T_chan']
    # Tc300   = heat300['T_chan']
    # Tc400   = heat400['T_chan']

    qf200   = heat200['qr_flux_iw']
    qf250   = heat250['qr_flux_iw']
    # qf300   = heat300['qr_flux_iw']
    # qf400   = heat400['qr_flux_iw']

    flow200 = d200['model_out']['flow_out']
    flow250 = d250['model_out']['flow_out']
    # flow300 = d300['model_out']['flow_out']
    # flow400 = d400['model_out']['flow_out']

    Tglob200 = heat200['Tall_glob']
    Tglob250 = heat250['Tall_glob']

    Tg200   = flow200['global']['T']
    Tg250   = flow250['global']['T']
    # Tg300   = flow300['global']['T']
    # Tg400   = flow400['global']['T']

    M200   = flow200['global']['M']
    M250   = flow250['global']['M']
    # M300   = flow300['global']['M']
    # M400   = flow400['global']['M']

    u200   = flow200['global']['u']
    u250   = flow250['global']['u']
    # u300   = flow300['global']['u']
    # u400   = flow400['global']['u']

    p200   = flow200['global']['p']
    p250   = flow250['global']['p']
    # p300   = flow300['global']['p']
    # p400   = flow400['global']['p']

    fig  = plt.figure(figsize = (10,8))
    ax1 = fig.add_subplot(111)
    plt.contourf(zc200[:,:nr], rc200[:,:nr], Tw200, levels=30)
    plt.contourf(zc200[:,:nr], -1*rc200[:,:nr], Tw200, levels=30)
    plt.colorbar()
    # plt.axis('equal')
    xl = plt.xlim()
    plt.plot(xl, [0,0], '--')

    fig  = plt.figure(figsize = (10,8))
    # ax2 = fig.add_subplot(212, sharex=ax1)
    # for ii in np.arange(nr):
        # plt.plot(zc200[:,0], Tglob200[:,ii])
    plt.scatter(zc200, rc200, c=Tglob200, s=1)
    plt.scatter(zc200, -1*rc200, c=Tglob200, s=1)
    plt.colorbar()
    plt.axis('equal')
    xl = plt.xlim()
    plt.plot(xl, [0,0], '--')



    fig = plt.figure(figsize=(20,10))
    ax1 = fig.add_subplot(221)
    plt.contourf(zc200[:,:nr], rc200[:,:nr], Tw200, levels=30)
    #
    ax1.axis('equal')
    plt.grid()
    plt.xlabel('$z$ [m]', **plot_opt['axis_font'])
    plt.ylabel('$r$ [m]', **plot_opt['axis_font'])
    plt.colorbar()

    ax2 = fig.add_subplot(222)
    plt.contourf(zc250[:,:nr], rc250[:,:nr], Tw250, levels=30)
    ax2.axis('equal')
    plt.grid()
    plt.xlabel('$z$ [m]', **plot_opt['axis_font'])
    plt.ylabel('$r$ [m]', **plot_opt['axis_font'])
    plt.colorbar()

    ax3 = fig.add_subplot(223, sharex=ax1)
    plt.plot(zc200[:,0], st200, '-r', label='s therm.')
    plt.plot(zc200[:,0], sp200, '-b', label='s press.')
    plt.plot(zc200[:,0], s200, '-k', label='s tot')
    # plt.plot(zc250[:,0], s250, '-r', label=r'With stress cons.')
    # plt.plot(zc300[:,0], s300, '-g', label='340s')
    # plt.plot(zc400[:,0], s400, '-b', label='350s')
    xl = plt.xlim()
    plt.plot(xl, [sy, sy], '--k')
    plt.grid()
    plt.xlabel('$z$ [m]', **plot_opt['axis_font'])
    plt.ylabel('$\sigma$ [Pa] ', **plot_opt['axis_font'])
    plt.legend(prop=plot_opt['cb_font'])
    # fig.tight_layout()

    ax4 = fig.add_subplot(224, sharex=ax2)
    plt.plot(zc250[:,0], st250, '-r', label='s therm.')
    plt.plot(zc250[:,0], sp250, '-b', label='s press.')
    plt.plot(zc250[:,0], s250, '-k', label='s tot')
    # plt.plot(zc250[:,0], s250, '-r', label=r'With stress cons.')
    # plt.plot(zc300[:,0], s300, '-g', label='340s')
    # plt.plot(zc400[:,0], s400, '-b', label='350s')
    xl = plt.xlim()
    plt.plot(xl, [sy, sy], '--k')
    plt.grid()
    plt.xlabel('$z$ [m]', **plot_opt['axis_font'])
    plt.ylabel('$\sigma$ [Pa] ', **plot_opt['axis_font'])
    plt.legend(prop=plot_opt['cb_font'])
    fig.tight_layout()

    fn_save = os.path.join('tex',r'thrust_isp_stress_constr_contourf_compare.png')
    plt.savefig(fn_save, **plot_opt['save'])
    plt.close()





    return

    fig = plt.figure(figsize=(20,12))
    ax1 = fig.add_subplot(231)
    plt.plot(zc200[:,0], Tg200, '-k', label=r'No stress cons.')
    plt.plot(zc250[:,0], Tg250, '-r', label=r'With stress cons.')
    # plt.plot(zc300[:,0], Tg300, '-g', label='340s')
    # plt.plot(zc400[:,0], Tg400, '-b', label='350s')
    plt.grid()
    plt.xlabel('$z$ [m]', **plot_opt['axis_font'])
    plt.ylabel('$T_{gas}$ [K]', **plot_opt['axis_font'])
    plt.legend(prop=plot_opt['cb_font'])

    ax2 = fig.add_subplot(233)
    plt.plot(zc200[:,0], qf200, '-k', label=r'No stress cons.')
    plt.plot(zc250[:,0], qf250, '-r', label=r'With stress cons.')
    # plt.plot(zc300[:,0], qf300, '-g', label='340s')
    # plt.plot(zc400[:,0], qf400, '-b', label='350s')
    plt.grid()
    plt.xlabel('$z$ [m]', **plot_opt['axis_font'])
    plt.ylabel('$q`` [W/m^{2}]$  ', **plot_opt['axis_font'])
    plt.legend(prop=plot_opt['cb_font'])


    ax3 = fig.add_subplot(232)
    plt.plot(zc200[:,0], u200, '-k', label=r'No stress cons.')
    plt.plot(zc250[:,0], u250, '-r', label=r'With stress cons.')
    # plt.plot(zc300[:,0], u300, '-g', label='340s')
    # plt.plot(zc400[:,0], u400, '-b', label='350s')
    plt.grid()
    plt.xlabel('$z$ [m]', **plot_opt['axis_font'])
    plt.ylabel('$u$ [m/s] ', **plot_opt['axis_font'])
    plt.legend(prop=plot_opt['cb_font'])
    sy  = material_prop['s_yield']

    ax4 = fig.add_subplot(234)
    plt.plot(zc200[:,0], Tw200[:,0], '-k', label=r'No stress cons.')
    plt.plot(zc250[:,0], Tw250[:,0], '-r', label=r'With stress cons.')
    # plt.plot(zc300[:,0], Tw300[:,0], '-g', label='340s')
    # plt.plot(zc400[:,0], Tw400[:,0], '-b', label='350s')
    plt.grid()
    plt.xlabel('$z$ [m]', **plot_opt['axis_font'])
    plt.ylabel('$T_{w,i}$ [K]', **plot_opt['axis_font'])
    plt.legend(prop=plot_opt['cb_font'])

    ax5 = fig.add_subplot(235)
    plt.plot(zc200[:,0], Tc200, '-k', label=r'No stress cons.')
    plt.plot(zc250[:,0], Tc250, '-r', label=r'With stress cons.')
    # plt.plot(zc300[:,0], Tc300, '-g', label='340s')
    # plt.plot(zc400[:,0], Tc400, '-b', label='350s')
    plt.grid()
    plt.xlabel('$z$ [m]', **plot_opt['axis_font'])
    plt.ylabel('$T_{c}$ [K]', **plot_opt['axis_font'])
    plt.legend(prop=plot_opt['cb_font'])

    ax6 = fig.add_subplot(236)
    plt.plot(zc200[:,0], s200, '-k', label=r'No stress cons.')
    plt.plot(zc250[:,0], s250, '-r', label=r'With stress cons.')
    # plt.plot(zc300[:,0], s300, '-g', label='340s')
    # plt.plot(zc400[:,0], s400, '-b', label='350s')
    xl = plt.xlim()
    plt.plot(xl, [sy, sy], '--k')
    plt.grid()
    plt.xlabel('$z$ [m]', **plot_opt['axis_font'])
    plt.ylabel('$\sigma$ [Pa] ', **plot_opt['axis_font'])
    plt.legend(prop=plot_opt['cb_font'])
    fig.tight_layout()

    fn_save = os.path.join('tex',r'thrust_isp_stress_constr_compare.png')
    plt.savefig(fn_save, **plot_opt['save'])
    plt.close()



def compare_isp_thrust():
    plot_opt = get_plot_opt()
    read_dir    = 'model_output_norm'

    fn1     = os.path.join(read_dir, 'thrust3.000e+05_mass_higherplim.pickle')
    fn2     = os.path.join(read_dir, 'isp325_thrust3.000e+05_mass.pickle')
    fn3     = os.path.join(read_dir, 'isp340_thrust3.000e+05_mass.pickle')
    fn4     = os.path.join(read_dir, 'isp350_thrust3.000e+05_mass.pickle')

    d200      = pickle.load(open(fn1, 'rb'))
    d250      = pickle.load(open(fn2, 'rb'))
    d300      = pickle.load(open(fn3, 'rb'))
    d400      = pickle.load(open(fn4, 'rb'))

    mesh200 = d200['model_out']['mesh_out']
    mesh250 = d250['model_out']['mesh_out']
    mesh300 = d300['model_out']['mesh_out']
    mesh400 = d400['model_out']['mesh_out']

    zc200   = mesh200['zones']['global']['zc']
    zc250   = mesh250['zones']['global']['zc']
    zc300   = mesh300['zones']['global']['zc']
    zc400   = mesh400['zones']['global']['zc']

    heat200 = d200['model_out']['heat_out']
    heat250 = d250['model_out']['heat_out']
    heat300 = d300['model_out']['heat_out']
    heat400 = d400['model_out']['heat_out']

    stress200 = d200['model_out']['stress_out']
    stress250 = d250['model_out']['stress_out']
    stress300 = d300['model_out']['stress_out']
    stress400 = d400['model_out']['stress_out']
    material_prop   = d200['model_out']['material_prop']
    sy  = material_prop['s_yield']

    s200    = abs( stress200['s_tot'])
    s250    = abs( stress250['s_tot'])
    s300    = abs( stress300['s_tot'])
    s400    = abs( stress400['s_tot'])

    Tw200   = heat200['Tw']
    Tw250   = heat250['Tw']
    Tw300   = heat300['Tw']
    Tw400   = heat400['Tw']

    Tc200   = heat200['T_chan']
    Tc250   = heat250['T_chan']
    Tc300   = heat300['T_chan']
    Tc400   = heat400['T_chan']

    qf200   = heat200['qr_flux_iw']
    qf250   = heat250['qr_flux_iw']
    qf300   = heat300['qr_flux_iw']
    qf400   = heat400['qr_flux_iw']

    flow200 = d200['model_out']['flow_out']
    flow250 = d250['model_out']['flow_out']
    flow300 = d300['model_out']['flow_out']
    flow400 = d400['model_out']['flow_out']

    Tg200   = flow200['global']['T']
    Tg250   = flow250['global']['T']
    Tg300   = flow300['global']['T']
    Tg400   = flow400['global']['T']

    M200   = flow200['global']['M']
    M250   = flow250['global']['M']
    M300   = flow300['global']['M']
    M400   = flow400['global']['M']

    u200   = flow200['global']['u']
    u250   = flow250['global']['u']
    u300   = flow300['global']['u']
    u400   = flow400['global']['u']

    p200   = flow200['global']['p']
    p250   = flow250['global']['p']
    p300   = flow300['global']['p']
    p400   = flow400['global']['p']

    fig = plt.figure(figsize=(20,12))
    ax1 = fig.add_subplot(231)
    plt.plot(zc200[:,0], Tg200, '-k', label='No Isp cons.')
    plt.plot(zc250[:,0], Tg250, '-r', label='325s')
    plt.plot(zc300[:,0], Tg300, '-g', label='340s')
    plt.plot(zc400[:,0], Tg400, '-b', label='350s')
    plt.grid()
    plt.xlabel('$z$ [m]', **plot_opt['axis_font'])
    plt.ylabel('$T_{gas}$ [K]', **plot_opt['axis_font'])
    plt.legend(prop=plot_opt['cb_font'])

    ax2 = fig.add_subplot(233)
    plt.plot(zc200[:,0], qf200, '-k', label='No Isp cons.')
    plt.plot(zc250[:,0], qf250, '-r', label='325s')
    plt.plot(zc300[:,0], qf300, '-g', label='340s')
    plt.plot(zc400[:,0], qf400, '-b', label='350s')
    plt.grid()
    plt.xlabel('$z$ [m]', **plot_opt['axis_font'])
    plt.ylabel('$q`` [W/m^{2}]$  ', **plot_opt['axis_font'])
    plt.legend(prop=plot_opt['cb_font'])


    ax3 = fig.add_subplot(232)
    plt.plot(zc200[:,0], u200, '-k', label='No Isp cons.')
    plt.plot(zc250[:,0], u250, '-r', label='325s')
    plt.plot(zc300[:,0], u300, '-g', label='340s')
    plt.plot(zc400[:,0], u400, '-b', label='350s')
    plt.grid()
    plt.xlabel('$z$ [m]', **plot_opt['axis_font'])
    plt.ylabel('$u$ [m/s] ', **plot_opt['axis_font'])
    plt.legend(prop=plot_opt['cb_font'])


    sy  = material_prop['s_yield']

    ax4 = fig.add_subplot(234)
    plt.plot(zc200[:,0], Tw200[:,0], '-k', label='No Isp cons.')
    plt.plot(zc250[:,0], Tw250[:,0], '-r', label='325s')
    plt.plot(zc300[:,0], Tw300[:,0], '-g', label='340s')
    plt.plot(zc400[:,0], Tw400[:,0], '-b', label='350s')
    plt.grid()
    plt.xlabel('$z$ [m]', **plot_opt['axis_font'])
    plt.ylabel('$T_{w,i}$ [K]', **plot_opt['axis_font'])
    plt.legend(prop=plot_opt['cb_font'])

    ax5 = fig.add_subplot(235)
    plt.plot(zc200[:,0], Tc200, '-k', label='No Isp cons.')
    plt.plot(zc250[:,0], Tc250, '-r', label='325s')
    plt.plot(zc300[:,0], Tc300, '-g', label='340s')
    plt.plot(zc400[:,0], Tc400, '-b', label='350s')
    plt.grid()
    plt.xlabel('$z$ [m]', **plot_opt['axis_font'])
    plt.ylabel('$T_{c}$ [K]', **plot_opt['axis_font'])
    plt.legend(prop=plot_opt['cb_font'])

    ax6 = fig.add_subplot(236)
    plt.plot(zc200[:,0], s200, '-k', label='No Isp cons.')
    plt.plot(zc250[:,0], s250, '-r', label='325s')
    plt.plot(zc300[:,0], s300, '-g', label='340s')
    plt.plot(zc400[:,0], s400, '-b', label='350s')
    xl = plt.xlim()
    plt.plot(xl, [sy, sy], '--k')
    plt.grid()
    plt.xlabel('$z$ [m]', **plot_opt['axis_font'])
    plt.ylabel('$\sigma$ [Pa] ', **plot_opt['axis_font'])
    plt.legend(prop=plot_opt['cb_font'])
    fig.tight_layout()

    fn_save = os.path.join('tex','thrust_isp_constr_compare.png')
    plt.savefig(fn_save, **plot_opt['save'])
    plt.close()


def compare_thrust_constr():
    plt.switch_backend('Qt5Agg')
    plot_opt = get_plot_opt()
    read_dir    = 'model_output_norm'
    fn200       = os.path.join(read_dir, 'thrust2.000e+05_mass_sol1.pickle')
    fn250       = os.path.join(read_dir, 'thrust2.500e+05_mass.pickle')
    fn300       = os.path.join(read_dir, 'thrust3.000e+05_mass.pickle')
    fn400       = os.path.join(read_dir, 'thrust4.000e+05_mass.pickle')

    d200    = pickle.load(open(fn200, 'rb'))
    d250    = pickle.load(open(fn250, 'rb'))
    d300    = pickle.load(open(fn300, 'rb'))
    d400    = pickle.load(open(fn400, 'rb'))

    mesh200 = d200['model_out']['mesh_out']
    mesh250 = d250['model_out']['mesh_out']
    mesh300 = d300['model_out']['mesh_out']
    mesh400 = d400['model_out']['mesh_out']

    zc200   = mesh200['zones']['global']['zc']
    zc250   = mesh250['zones']['global']['zc']
    zc300   = mesh300['zones']['global']['zc']
    zc400   = mesh400['zones']['global']['zc']

    heat200 = d200['model_out']['heat_out']
    heat250 = d250['model_out']['heat_out']
    heat300 = d300['model_out']['heat_out']
    heat400 = d400['model_out']['heat_out']

    stress200 = d200['model_out']['stress_out']
    stress250 = d250['model_out']['stress_out']
    stress300 = d300['model_out']['stress_out']
    stress400 = d400['model_out']['stress_out']
    material_prop   = d200['model_out']['material_prop']
    sy  = material_prop['s_yield']

    s200    = abs( stress200['s_tot'])
    s250    = abs( stress250['s_tot'])
    s300    = abs( stress300['s_tot'])
    s400    = abs( stress400['s_tot'])

    Tw200   = heat200['Tw']
    Tw250   = heat250['Tw']
    Tw300   = heat300['Tw']
    Tw400   = heat400['Tw']

    Tc200   = heat200['T_chan']
    Tc250   = heat250['T_chan']
    Tc300   = heat300['T_chan']
    Tc400   = heat400['T_chan']

    qf200   = heat200['qr_flux_iw']
    qf250   = heat250['qr_flux_iw']
    qf300   = heat300['qr_flux_iw']
    qf400   = heat400['qr_flux_iw']

    flow200 = d200['model_out']['flow_out']
    flow250 = d250['model_out']['flow_out']
    flow300 = d300['model_out']['flow_out']
    flow400 = d400['model_out']['flow_out']

    Tg200   = flow200['global']['T']
    Tg250   = flow250['global']['T']
    Tg300   = flow300['global']['T']
    Tg400   = flow400['global']['T']

    M200   = flow200['global']['M']
    M250   = flow250['global']['M']
    M300   = flow300['global']['M']
    M400   = flow400['global']['M']

    u200   = flow200['global']['u']
    u250   = flow250['global']['u']
    u300   = flow300['global']['u']
    u400   = flow400['global']['u']

    p200   = flow200['global']['p']
    p250   = flow250['global']['p']
    p300   = flow300['global']['p']
    p400   = flow400['global']['p']

    fig = plt.figure(figsize=(20,12))
    ax1 = fig.add_subplot(231)
    plt.plot(zc200[:,0], Tg200, '-k', label='200kN')
    plt.plot(zc250[:,0], Tg250, '-r', label='250kN')
    plt.plot(zc300[:,0], Tg300, '-g', label='300kN')
    plt.plot(zc400[:,0], Tg400, '-b', label='400kN')
    plt.grid()
    plt.xlabel('$z$ [m]', **plot_opt['axis_font'])
    plt.ylabel('$T_{gas}$ [K]', **plot_opt['axis_font'])
    plt.legend(prop=plot_opt['cb_font'])

    ax2 = fig.add_subplot(233)
    plt.plot(zc200[:,0], qf200, '-k', label='200kN')
    plt.plot(zc250[:,0], qf250, '-r', label='250kN')
    plt.plot(zc300[:,0], qf300, '-g', label='300kN')
    plt.plot(zc400[:,0], qf400, '-b', label='400kN')
    plt.grid()
    plt.xlabel('$z$ [m]', **plot_opt['axis_font'])
    plt.ylabel('$q`` [W/m^{2}]$  ', **plot_opt['axis_font'])
    plt.legend(prop=plot_opt['cb_font'])


    ax3 = fig.add_subplot(232)
    plt.plot(zc200[:,0], u200, '-k', label='200kN')
    plt.plot(zc250[:,0], u250, '-r', label='250kN')
    plt.plot(zc300[:,0], u300, '-g', label='300kN')
    plt.plot(zc400[:,0], u400, '-b', label='400kN')
    plt.grid()
    plt.xlabel('$z$ [m]', **plot_opt['axis_font'])
    plt.ylabel('$u$ [m/s] ', **plot_opt['axis_font'])
    plt.legend(prop=plot_opt['cb_font'])
    # fig.tight_layout()

    # fn_save = os.path.join('tex','Tgas_ugas_heatflux_thrust_constr.png')
    # plt.savefig(fn_save, **plot_opt['save'])
    # plt.close()
    # return



    sy  = material_prop['s_yield']


    ax4 = fig.add_subplot(234)
    plt.plot(zc200[:,0], Tw200[:,0], '-k', label='200kN')
    plt.plot(zc250[:,0], Tw250[:,0], '-r', label='250kN')
    plt.plot(zc300[:,0], Tw300[:,0], '-g', label='300kN')
    plt.plot(zc400[:,0], Tw400[:,0], '-b', label='400kN')
    plt.grid()
    plt.xlabel('$z$ [m]', **plot_opt['axis_font'])
    plt.ylabel('$T_{w,i}$ [K]', **plot_opt['axis_font'])
    plt.legend(prop=plot_opt['cb_font'])

    ax5 = fig.add_subplot(235)
    plt.plot(zc200[:,0], Tc200, '-k', label='200kN')
    plt.plot(zc250[:,0], Tc250, '-r', label='250kN')
    plt.plot(zc300[:,0], Tc300, '-g', label='300kN')
    plt.plot(zc400[:,0], Tc400, '-b', label='400kN')
    plt.grid()
    plt.xlabel('$z$ [m]', **plot_opt['axis_font'])
    plt.ylabel('$T_{c}$ [K]', **plot_opt['axis_font'])
    plt.legend(prop=plot_opt['cb_font'])

    ax6 = fig.add_subplot(236)
    plt.plot(zc200[:,0], s200, '-k', label='200kN')
    plt.plot(zc250[:,0], s250, '-r', label='250kN')
    plt.plot(zc300[:,0], s300, '-g', label='300kN')
    plt.plot(zc400[:,0], s400, '-b', label='400kN')
    xl = plt.xlim()
    plt.plot(xl, [sy, sy], '--k')
    plt.grid()
    plt.xlabel('$z$ [m]', **plot_opt['axis_font'])
    plt.ylabel('$\sigma$ [Pa] ', **plot_opt['axis_font'])
    plt.legend(prop=plot_opt['cb_font'])
    fig.tight_layout()

    fn_save = os.path.join('tex','thrust_constr_compare.png')
    plt.savefig(fn_save, **plot_opt['save'])
    plt.close()




    # test = 1

    # #plot the inner wall temperature profiles
    # zc_glob = mesh_out['zones']['global']['zc']
    # rc_glob = mesh_out['zones']['global']['rc']
    # Tw = heat_out['Tw']
    # Tall_glob = heat_out['Tall_glob']



def get_plot_opt():
    '''Get options for plotting, set some plt.rc options
    
    NOTE: use of this function requires installation of Latex on the
    machine under consideration.
    
    Args:
        None
        
    Returns:
        opt (dict) : contains font settings for building plots    
    '''
    
    plt.rc('text', usetex=True)
    axis_font = {'size':'22'}
    title_font = {'size':'24'}
    tick_font = {'size':'22'}
    cb_font = {'size':'20'}
    save = {'bbox_inches':'tight', 'dpi':150}

    opt = {'axis_font':axis_font, 
           'title_font':title_font,
           'tick_font':tick_font,
           'cb_font':cb_font,
           'save':save}
    
    plt.rc('xtick', labelsize=tick_font['size'])
    plt.rc('ytick', labelsize=tick_font['size'])
    plt.rcParams['font.family'] = 'Times New Roman'
    return opt

### if name == __main__
if __name__ == '__main__':
    # compare_thrust_constr()
    # data = load_test()
    # compare_isp_thrust()
    # compare_isp_thrust_stress()
    compare_all_constr()




