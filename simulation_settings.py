import math
import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import constants


# df1=pd.read_csv("Fx.csv",header=-1) #fnq= 41 oq=61  df1[fnq][oq]
# df3=pd.read_csv("Fz.csv",header=-1)
# df5=pd.read_csv("My.csv",header=-1)
# matrixList=[df1,df3,df5]
# ddf1=pd.read_csv("Fxd.csv",header=-1)
# ddf3=pd.read_csv("Fzd.csv",header=-1)
# ddf5=pd.read_csv("Myd.csv",header=-1)
# dmatrixList=[ddf1,ddf3,ddf5]

df = pd.read_csv('ANALYSIS_forces_180.csv')

matrixList = [df['Fx[N/m]'],df['Fz[N/m]'],df['My[Nm/m]'] ]
dmatrixList = [df['Fx[deg]'],df['Fz[deg]'],df['My[deg]'] ]

gg = constants.gg
rho = constants.rho
beta = constants.beta
dt = constants.dt

class Wave:
    amp = 0.
    k = 0.
    omega = 0.
    phase = 0.
    lamda = 0.
    period = 0.
    length = 0.
    timer = 5

    def __init__(self, a, o, p):
        self.amp = a
        self.omega = o
        self.k = self.omega**2/gg
        self.phase = p
        self.period = 2*math.pi/o
        self.length = 2*math.pi/self.k
    
    def get(self,x, t):
        if(t<self.timer):
            return (t/self.timer)*self.amp*math.sin(self.omega*t + self.k*x + self.phase)
        else:
            return self.amp*math.sin(self.omega*t + self.k*x + self.phase)

    def calcF(self,x,t,axis):
        # fni = round(Fn/(1./40)) #0~40 
        # omi = round(omega_e/(1./4)) #0~60

        # fni =abs(fni)
        # omi =abs(omi)
        # if(fni>40): 
        #     # print("Warning! Fn index over 40")
        #     fni=40
        # if(omi>60): 
        #     # print("Warning! Omega_e index over 60")
        #     omi=60

        # j = round((axis-1)/2)

        # coeff =  matrixList[j][fni][omi]*100 #N/cm -> N/m
        # poffset= dmatrixList[j][fni][omi]*math.pi/180 #deg -> rad
        # output = coeff*self.amp*math.sin(self.omega*t + self.k*x + self.phase + poffset)  
        
        pi = 39-math.floor((self.period-0.35)/0.1699)
        if(pi>39): 
            # pi=39
            print(f'period:\t{self.period}')
            return 0
        if(pi<0): 
            # pi = 0
            print(f'period:\t{self.period}')
            return 0

        j = round((axis-1)/2)
        coeff =  matrixList[j][pi] #N/m or Nm/m
        poffset= dmatrixList[j][pi]*math.pi/180 #deg -> rad
        output = self.amp*coeff*math.sin(self.omega*t + self.k*x + self.phase - poffset)  
        
        if(t<self.timer):
            return (t/self.timer)*output
        else:
            return output #N


class JONSWAP:
    N = 0    
    T_1 = 0
    H_third = 0
    waves = []
    delta_omega = 0
    max_omega = 0
    min_omega = 0
    timer = 1
    id=0
    def __init__(self, num_waves, T_1, H_third):
        self.waves = []
        self.N = num_waves
        self.T_1 = T_1
        self.H_third = H_third
        self.max_omega = 8.0*math.pi/self.T_1 
        self.min_omega = 1.0*math.pi/self.T_1 

        self.delta_omega = (self.max_omega-self.min_omega) /self.N

        for i in range(self.N):
            seed = i+self.id
            random.seed(seed)
            new_omega = self.min_omega + self.delta_omega*(0.5+i)+random.uniform(-self.delta_omega/2,self.delta_omega/2)
            random.seed(seed)
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
    def getVel(self,x,t):
        pvelx = 0
        pvelz = 0
        for w in self.waves:
            # z = 
            coeff = w.omega * w.amp * math.exp(w.k*0)
            pvelx += coeff*math.sin(w.omega*t-w.k*x)
            pvelz += coeff*math.cos(w.omega*t-w.k*x)
        return pvelx,pvelz

# class JONSWAP:
#     N = 0    
#     T_1 = 0
#     H_third = 0
#     waves = []
#     delta_period = 0
#     max_period = 0
#     min_period = 0
#     def __init__(self, num_waves, T_1, H_third):
#         self.waves = []
#         self.N = num_waves
#         self.T_1 = T_1
#         self.H_third = H_third
#         self.max_period = 2.0*self.T_1 
#         self.min_period = 0.5*self.T_1 

#         self.delta_period = (self.max_period-self.min_period) /self.N

#         for i in range(self.N):
#             random.seed(i)
#             new_period = self.min_period + self.delta_period*(0.5+i)+random.uniform(-self.delta_period/2,self.delta_period/2)
#             random.seed(i)
#             new_wave = Wave(self.getWaveHeight(new_period),new_period,random.uniform(0,2*math.pi))
#             self.waves.append(new_wave)

#     def getWaveHeight(self,period):
#         if(period!=0):
#             sigma = 0.08
#             omega = 2*math.pi/period
#             Y = math.exp(-((0.191*omega*self.T_1-1)/(math.sqrt(2)*sigma))**2) 
#             S = 155*(self.H_third**2)/(self.T_1**4*omega**5)*math.exp(-944/(self.T_1**4*omega**4))*3.3**Y
#             A = math.sqrt(2*S*self.delta_omega)
#         else:
#             A = 0
#             S = 0
#         return A
#         # return S

#     def get(self,x,t):
#         waveheight = 0
#         for w in self.waves:
#             waveheight += w.get(x,t)
#         return waveheight
#     def getVel(self,x,t):
#         pvelx = 0
#         pvelz = 0
#         for w in self.waves:
#             # z = 
#             coeff = w.omega * w.amp * math.exp(w.k*0)
#             pvelx += coeff*math.sin(w.omega*t-w.k*x)
#             pvelz += coeff*math.cos(w.omega*t-w.k*x)
#         return pvelx,pvelz


class Water:
    waves = []
    def __init__(self, amp, omega):
        self.waves = []
        newWave = Wave(amp,omega,0)
        self.waves.append(newWave)
            
    def get(self,x,t):
        waveheight = 0
        for w in self.waves:
            waveheight += w.get(x,t)
        return waveheight