#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: jd
contact: jamesduv@umich.edu
affiliation: University of Michigan, Department of Aerospace Eng., CASLAB
"""
import numpy as np
import matplotlib.pyplot as plt

def plot_temp_contourf(heat_out, mesh_out, mesh_param):
    zc_glob = mesh_out['zones']['global']['zc']
    rc_glob = mesh_out['zones']['global']['rc']
    Tw = heat_out['Tw']
    Tall_glob = heat_out['Tall_glob']
    nr = mesh_param['nr']

    # qr_flux = heat_out['qr_flux_iw']
    # fig = plt.figure()
    # plt.plot(zc_glob[:,0], qr_flux, 'ko')
    # plt.title('heat flux vs z')


    fig = plt.figure()
    plt.plot(zc_glob[:,0], Tw[:,0], 'ko')
    plt.title('inner wall temp')

    fig = plt.figure()
    plt.contourf(zc_glob[:,:nr], rc_glob[:,:nr], Tw, levels=30, cmap=plt.cm.jet)
    plt.colorbar()
    plt.title('Inner wall')
    plt.axis('equal')

    fig = plt.figure()
    plt.contourf(zc_glob[:,:], rc_glob[:,:], Tall_glob, levels=30, cmap=plt.cm.jet)
    plt.colorbar()
    plt.title('All')
    plt.axis('equal')


def plot_stress(stress_out, mesh_out, material_prop):

    sy  = material_prop['s_yield']
    plot_opt = get_plot_opt()
    s_hoop  = stress_out['s_hoop']
    s_cte   = stress_out['s_cte']
    s_tot   = stress_out['s_tot']
    zc      = mesh_out['zones']['global']['zc'][:,0]

    fig = plt.figure(figsize=(10,5))
    plt.plot(zc, abs(s_hoop), '-b', label='s hoop')
    plt.plot(zc, abs(s_cte), '-r', label='s cte')
    plt.plot(zc, abs(s_tot), '-k', label='s tot')
    plt.plot([zc[0], zc[-1]], [sy, sy], '--r')
    plt.legend(prop=plot_opt['axis_font'])


def plot_flow_vs_z(flow_out, mesh_out):
    '''Plot various flow variables versus z'''

    plot_opt = get_plot_opt()
    zc = mesh_out['zones']['global']['zc'][:,0]
    flow_glob = flow_out['global']
    pass
    fig = plt.figure(figsize=(20,10))
    ax1 = fig.add_subplot(231)
    ax1.plot(zc, flow_glob['T'], '-k' )
    ax1.set_xlabel('$z$', **plot_opt['axis_font'])
    ax1.set_ylabel('$T(z)$', **plot_opt['axis_font'])
    ax1.set_title('$T(z)$', **plot_opt['axis_font'])
    ax1.grid()

    ax2 = fig.add_subplot(232)
    ax2.plot(zc, flow_glob['p'], '-k' )
    ax2.set_xlabel('$z$', **plot_opt['axis_font'])
    ax2.set_ylabel('$p(z)$', **plot_opt['axis_font'])
    ax2.set_title('$p(z)$', **plot_opt['axis_font'])
    ax2.grid()

    ax3 = fig.add_subplot(233)
    ax3.plot(zc, flow_glob['rho'], '-k' )
    ax3.set_xlabel('$z$', **plot_opt['axis_font'])
    ax3.set_ylabel(r'$\rho(z)$', **plot_opt['axis_font'])
    ax3.set_title(r'$\rho(z)$', **plot_opt['axis_font'])
    ax3.grid()

    ax4 = fig.add_subplot(234)
    ax4.plot(zc, flow_glob['M'], '-k' )
    ax4.set_xlabel('$z$', **plot_opt['axis_font'])
    ax4.set_ylabel('$M(z)$', **plot_opt['axis_font'])
    ax4.set_title('$M(z)$', **plot_opt['axis_font'])
    ax4.grid()

    ax5 = fig.add_subplot(235)
    ax5.plot(zc, flow_glob['u'], '-k' )
    ax5.set_xlabel('$z$', **plot_opt['axis_font'])
    ax5.set_ylabel('$u(z)$', **plot_opt['axis_font'])
    ax5.set_title('$u(z)$', **plot_opt['axis_font'])
    ax5.grid()

    ax6 = fig.add_subplot(236)
    ax6.plot(zc, flow_glob['Re'], '-k' )
    ax6.set_xlabel('$z$', **plot_opt['axis_font'])
    ax6.set_ylabel('$Re(z)$', **plot_opt['axis_font'])
    ax6.set_title('$Re(z)$', **plot_opt['axis_font'])
    ax6.grid()

    fig.tight_layout()

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