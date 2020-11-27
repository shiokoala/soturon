import math
import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import constants


df1=pd.read_csv("Fx.csv",header=-1) #fnq= 41 oq=61  df1[fnq][oq]
df3=pd.read_csv("Fz.csv",header=-1)
df5=pd.read_csv("My.csv",header=-1)
matrixList=[df1,df3,df5]
ddf1=pd.read_csv("Fxd.csv",header=-1)
ddf3=pd.read_csv("Fzd.csv",header=-1)
ddf5=pd.read_csv("Myd.csv",header=-1)
dmatrixList=[ddf1,ddf3,ddf5]


gg = constants.gg
rho = constants.rho
beta = constants.beta
dt = constants.dt

class Wave:
    amp = 0.
    k = 0.
    omega = 0.
    phase = 0.
    def __init__(self, a, o, p):
        self.amp = a
        self.omega = o
        self.k = o*o/gg
        self.phase = p
    
    def get(self,x, t):
        return self.amp*math.sin(self.omega*t + self.k*x + self.phase)

    def calcF(self,x,t,Fn,omega_e,axis):
        fni = round(Fn/(1./40)) #0~40 
        omi = round(omega_e/(1./4)) #0~60

        fni =abs(fni)
        omi =abs(omi)
        if(fni>40): 
            print("Warning! Fn index over 40")
            fni=40
        if(omi>60): 
            print("Warning! Omega_e index over 60")
            omi=60

        j = round((axis-1)/2)

        coeff =  matrixList[j][fni][omi]*100 #N/cm -> N/m
        poffset= dmatrixList[j][fni][omi]
        output = coeff*self.amp*math.sin(self.omega*t + self.k*x + self.phase + poffset)  
        return output #N


class JONSWAP:
    N = 0    
    T_1 = 0
    H_third = 0
    waves = []
    delta_omega = 0
    max_omega = 0
    def __init__(self, num_waves, T_1, H_third):
        self.waves = []
        self.N = num_waves
        self.T_1 = T_1
        self.H_third = H_third
        self.max_omega = 4*math.pi/self.T_1 * 4
        self.delta_omega = self.max_omega/self.N
        for i in range(self.N):
            new_omega = self.delta_omega*(0.5+i)+random.uniform(-self.delta_omega/2,self.delta_omega/2)
            new_wave = Wave(self.getWaveHeight(new_omega),new_omega,random.uniform(0,2*math.pi))
            self.waves.append(new_wave)

    def getWaveHeight(self,omega):
        if(omega!=0):
            sigma = 0.08
            Y = math.exp(-((0.191*omega*self.T_1-1)/(math.sqrt(2)*sigma))**2) 
            S = 155*(self.H_third**2)/(self.T_1**4*omega**5)*math.exp(-944/(self.T_1**4*omega**4))*3.3**Y
            A = math.sqrt(2*S*self.delta_omega)
        else:
            A = 0
            S = 0
        return A
        # return S

    def get(self,x,t):
        waveheight = 0
        for w in self.waves:
            waveheight += w.get(x,t)
        return waveheight


class Water:
    waves = []
    def __init__(self, num, random_wave):
        self.waves = []
        if(random_wave==True):
            for i in range(num):
                newWave = Wave(random.uniform(0.01,0.05),random.uniform(1.2,3),random.uniform(1,10))
                self.waves.append(newWave)
        elif(num==0):
            self.waves.append(Wave(0,0,0))
        else:
            self.waves.append(Wave(0.1,10,0))
            
    def get(self,x,t):
        waveheight = 0
        for w in self.waves:
            waveheight += w.get(x,t)
        return waveheight