# import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def variable_scaling():
    thermochem_vars = {'mdot_fuel':10}


def coor_transform():
    beta = 15 * np.pi / 180
    A = np.array([[1, np.sin(beta)],[0, np.cos(beta)]])
    b = np.array([1,1])
    test = np.linalg.solve(A,b)

    nvec = 10
    rvec = np.linspace(0, 1, nvec)
    zvec = np.linspace(0, 1, nvec)

    R, Z = np.meshgrid(rvec, zvec)
    Rflat = np.reshape(R, (nvec**2), order='F' )
    Zflat = np.reshape(Z, (nvec**2), order='F')

    rpvec = rvec - zvec * np.tan(beta)
    zpvec = zvec / np.cos(beta)

    RP, ZP = np.meshgrid(rpvec, zpvec)
    RPflat = np.reshape(RP, (nvec**2), order='F' )
    ZPflat = np.reshape(ZP, (nvec**2), order='F')

    fig = plt.figure()
    plt.plot(Rflat, Zflat, 'ko')
    plt.plot(RPflat, ZPflat, 'bx')

def build_therm_matrix():
    ns = 5
    nw = 10
    nr = ns + nw + 1
    Tinf = 300 #K





if __name__ == '__main__':
    coor_transform()