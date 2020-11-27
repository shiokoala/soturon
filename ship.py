import math
import simulation_settings as ss
import constants

gg = constants.gg
rho = constants.rho
beta = constants.beta
dt = constants.dt
# scaling = 25/3.11
scaling = 1

class Ship:
    posx = 0
    posz = 0
    velx = 0
    velz = 0
    accx = 0
    accz = 0
    
    angle = 0
    anglevel = 0
    angleacc = 0

    wave_incline = 0
    # omega_e = 0

    pl = 3.11*scaling
    ph = 0.34*scaling
    pw = 0.94*scaling
    pv = 0.0878*scaling*scaling*scaling
    mass = pv*100
    cgh = 0.340*scaling

    draught = 0

    A11 = 0.3 * rho * pv
    A13 = 0
    A31 = 0
    A53 = 0.1 * rho * pv * pl
    A35 = 0.05* rho * pv * pl
    A33 = 0.3 * rho * pv
    A51 = 0.01* rho * pv * pl
    A55 = 0.07* rho * pv * pl * pl
    B11 = 2.0* rho * pv * math.sqrt(gg/pl)
    B53 = 0.05* rho * pv * pl * math.sqrt(gg/pl)
    B35 = 0.4* rho * pv * pl * math.sqrt(gg/pl)
    B33 = 2.0* rho * pv * math.sqrt(gg/pl)
    B55 = 0.4* rho * pv * pl * pl * math.sqrt(gg/pl)
    C11 = 0
    C53 = 0
    C35 = 0
    C33 = 11014
    C55 = 6233

    def __init__(self,x,z):
        self.posx = x
        self.posz = z
        self.accx = 0
        self.accz = 0
        self.velx = 0
        self.velz = 0

    def calcAcc(self,ww,t,thr_f):
        tFx=0
        tFz=0
        tMy=0

        sum_draught = 0
        res = 3

        for i in range(res):
          cos = math.cos(self.angle)
          sin = math.sin(self.angle)
          dl = self.pl/2/(res-1)*cos
          dh = self.pl/2/(res-1)*sin
          if(i==0):
            waveheight = ww.get(self.posx,t)
            sum_draught += min(max(waveheight-(self.posz-self.ph),0), self.ph)
          else:
            xpos1 = self.posx+i*dl
            xpos2 = self.posx-i*dl
            zpos1 = self.posz+i*dh
            zpos2 = self.posz-i*dh
            waveheight1 = ww.get(xpos1,t)
            waveheight2 = ww.get(xpos2,t)
            h_tilt = self.ph/cos
            sum_draught += min(max(waveheight1-(zpos1-h_tilt),0), h_tilt)
            sum_draught += min(max(waveheight2-(zpos2-h_tilt),0), h_tilt)

            if(i==res-1): 
              self.wave_incline = math.atan2((waveheight1-waveheight2),self.pl*cos)

        self.draught = sum_draught/(2*res-1)
            
        nx = math.cos(self.angle - math.pi/2)
        ny = math.sin(self.angle - math.pi/2)

        Fn = math.fabs(self.velx)/math.sqrt(9.8*self.pl)

        tFz -= self.mass*gg #Gravity
        
        if(self.draught>0):
            tFx += thr_f * math.cos(self.angle)
            tFz += thr_f * math.sin(self.angle)
            tMy += thr_f / self.mass

            for w in ww.waves:
                omega_e = w.omega + w.omega**2 * self.velx / gg
                # tFx += w.calcF(self.posx,t,Fn,omega_e,axis=1)
                # tFz += w.calcF(self.posx,t,Fn,omega_e,axis=3)
                # tMy += w.calcF(self.posx,t,Fn,omega_e,axis=5)
                tFx += w.calcF(self.posx,t,Fn,omega_e,axis=1)*self.draught/self.ph
                tFz += w.calcF(self.posx,t,Fn,omega_e,axis=3)*self.draught/self.ph
                tMy += w.calcF(self.posx,t,Fn,omega_e,axis=5)*self.draught/self.ph

            z_t = self.posz - 0.185
            heave = (tFz - self.B33*(self.velz + dt/2*self.accz) - self.C33*(z_t+dt*self.velz+(0.5-beta)*dt*dt*self.accz)) / ((self.mass+self.A33) + dt/2*self.B33 + beta*dt*dt*self.C33)
            surge = (tFx - self.B11*(self.velx + dt/2*self.accx)) / ((self.mass+self.A11) + dt/2*self.B11)
            pitch = (tMy - self.B55*(self.anglevel + dt/2*self.angleacc) - self.C55*(self.angle+dt*self.anglevel+(0.5-beta)*dt*dt*self.angleacc)) / ((self.mass+self.A55) + dt/2*self.B55 + beta*dt*dt*self.C55)

            self.accx = surge
            self.accz = heave
            self.angleacc = pitch
        else:
          self.accx = 0
          self.accz = -9.8

    def update(self,ww,t,thr_f):
      self.calcAcc(ww,t,thr_f)
      self.velx += self.accx*dt
      self.velz += self.accz*dt
      self.anglevel += self.angleacc*dt
      self.posx += self.velx*dt
      self.posz += self.velz*dt
      self.angle += self.anglevel*dt
      # print("x=" + "{:.2f}".format(self.posx))
      # print("z=" + "{:.2f}".format(self.posz))
