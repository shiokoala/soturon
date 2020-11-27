import simulation_settings as ss
import ship as ship
import constants
import math
import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patches as patches
import glob
import os
import driver

gg = constants.gg
rho = constants.rho
beta = constants.beta
dt = constants.dt


def calcFitness(ww,dd):
    pp = ship.Ship(0,0.18)
    dt = constants.dt
    t=0.0
    num =1000
    dl = []
    xl = []
    xvl = []
    xal = []
    zl = []
    zvl= []
    zal =[]
    eng = []
    pl = []
    wil = []
    fl =[]
    fps = 50
    limit = round(round(1./dt)/fps)
    prevForce= 0
    for i in range(num):
        input_array = np.array([pp.accx,
                                pp.accz,
                                pp.velx,
                                pp.velz,
                                pp.angle,
                                prevForce])
        force = dd.output(input_array)
        prevForce = force
        pp.update(ww,t,force)

        dl.append(pp.draught)
        xl.append(pp.posx)
        xvl.append(pp.velx)
        xal.append(pp.accx)
        zl.append(pp.posz)
        zvl.append(pp.velz)
        zal.append(pp.accz)
        wil.append(pp.wave_incline)
        pl.append(pp.angle*(180/math.pi))
        fl.append(force)

        t+=dt

    