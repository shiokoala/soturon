import math
import simulation_settings as ss
import constants

gg = constants.gg
rho = constants.rho
beta = constants.beta
dt = constants.dt
# scaling = 25/3.11
scaling = 1

def pwm_curve(pwm_us):
    rev = 0
    if(1475<pwm_us and pwm_us<1525):
        rev = 0
    elif(pwm_us>1500):
        if(pwm_us>1900): pwm_us=1900
        rev = -2588 + 4.11*pwm_us - 2.19*10**-3 * pwm_us**2 + 4.01*10**-7 * pwm_us**3
    elif(pwm_us<1500):
        if(pwm_us<1100): pwm_us=1100
        rev =   775 - 1.69*pwm_us + 1.38*10**-3 * pwm_us**2 - 3.99*10**-7 * pwm_us**3
    return rev


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
    pitchradius = 0.782
    inertialmoment = mass*pitchradius**2
    cgh = 0.340*scaling

    draught = 0

    # A11 = 0.3 * rho * pv * 0.01 
    A11 = mass*0.01 # "Added mass in surge motion is relatively low in comparison to ship mass"
    #https://www.tandfonline.com/doi/full/10.1080/17445302.2019.1615705#:~:text=The%20added%20mass%20is%20referred,considerable%20influence%20on%20ship%20manoeuvrability.
    # A13 = 0
    # A31 = 0
    # A53 = 0.1 * rho * pv * pl
    # A35 = 0.05* rho * pv * pl
    A33 = 0.3 * rho * pv
    # A51 = 0.01* rho * pv * pl
    A55 = 0.07* rho * pv * pl * pl
    B11 = 2.0* rho * pv * math.sqrt(gg/pl)*0.17
    # B53 = 0.05* rho * pv * pl * math.sqrt(gg/pl)
    # B35 = 0.4* rho * pv * pl * math.sqrt(gg/pl)
    B33 = 2.0* rho * pv * math.sqrt(gg/pl)
    B55 = 0.4* rho * pv * pl * pl * math.sqrt(gg/pl)
    # C11 = 0
    # C53 = 0
    # C35 = 0
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
        res = 5
        average_wave_height = 0

        for i in range(res):
            cos = math.cos(self.angle)
            sin = math.sin(self.angle)
            dl = self.pl/2/(res-1)*cos
            dh = self.pl/2/(res-1)*sin
            if(i==0):
                waveheight = ww.get(self.posx,t)
                sum_draught += min(max(waveheight-(self.posz-self.ph),0), self.ph)
                average_wave_height += waveheight
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
                average_wave_height += waveheight1
                average_wave_height += waveheight2
            if(i==res-1): 
                self.wave_incline = math.atan2((waveheight1-waveheight2),self.pl*cos)

        self.draught = sum_draught/(2*res-1)
        average_wave_height = average_wave_height/(2*res-1)
            
        nx = math.cos(self.angle - math.pi/2)
        ny = math.sin(self.angle - math.pi/2)

        Fn = math.fabs(self.velx)/math.sqrt(9.8*self.pl)

        #Thruster Force 
        # hz = pwm_curve(pwm_us)
        # particle_velx = w.getVelx(x,t) 
        # U_rel = (self.velx-particle_velx) 
        # U_in = -U_rel
        # if(hz>0):
        #     thr_f = -0.286*U_in*hz + 0.01375*hz**2
        # else:
        #     thr_f = -0.246*U_in*hz - 0.01185*hz**2

        tFx += thr_f * math.cos(self.angle)
        tFz += thr_f * math.sin(self.angle)
        tMy += thr_f * self.ph

        #Hydrodynamic Force
        for w in ww.waves:
            omega_e = w.omega + w.omega**2 * self.velx / gg
            k= 0.1
            # tFx += w.calcF(self.posx,t,Fn,omega_e,axis=1)*math.exp(-k*(self.posz - 0.185 - average_wave_height))
            # tFz += w.calcF(self.posx,t,Fn,omega_e,axis=3)*math.exp(-k*(self.posz - 0.185 - average_wave_height))
            # tMy += w.calcF(self.posx,t,Fn,omega_e,axis=5)*math.exp(-k*(self.posz - 0.185 - average_wave_height))
            tFx += w.calcF(self.posx,t,Fn,omega_e,axis=1)*self.draught/self.ph
            tFz += w.calcF(self.posx,t,Fn,omega_e,axis=3)*self.draught/self.ph
            tMy += w.calcF(self.posx,t,Fn,omega_e,axis=5)*self.draught/self.ph

        z_t = self.posz - 0.185 - average_wave_height
        heave = (tFz - self.B33*(self.velz + dt/2*self.accz) - 1*self.C33*(z_t+dt*self.velz+(0.5-beta)*dt*dt*self.accz)) / ((self.mass+self.A33) + dt/2*self.B33 + beta*dt*dt*self.C33)
        surge = (tFx - self.B11*(self.velx + dt/2*self.accx)) / ((self.mass+self.A11) + dt/2*self.B11)
        pitch = (tMy - self.B55*(self.anglevel + dt/2*self.angleacc) - self.C55*(self.angle+dt*self.anglevel+(0.5-beta)*dt*dt*self.angleacc)) / ((self.inertialmoment+self.A55) + dt/2*self.B55 + beta*dt*dt*self.C55)


        self.accx = surge
        self.accz = heave
        self.angleacc = pitch


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



class Ship2:
    posx = 0
    posz = 0
    velx = 0
    velz = 0
    accx = 0
    accz = 0

    IMU_accx = 0
    IMU_accz = 0
    IMU_roty = 0
    
    angle       = 0
    anglevel    = 0
    angleacc    = 0

    wave_incline = 0

    pl          = 3.11*scaling
    ph          = 0.34*scaling
    pw          = 0.94*scaling
    pv          = 0.0878*scaling*scaling*scaling
    mass        = 86.80879*scaling*scaling*scaling
    pitchradius = 0.782*scaling
    im          = mass*pitchradius**2
    cgh         = 0.340*scaling
    cgl         = 1.41 *scaling

    power_type      = ''
    prop_type       = ''
    prop_rev        = 0 #hz
    prev_rev        = 0
    rev_err         = 0
    prop_d          = 76/1000 * scaling
    control_delay   = 0
    control_rate    = 0
    control_start   = 0
    time_constant   = 0
    delayed_rate    = 0
    transition_delay= 0
    target_rev      = 0
    Pg = 0
    Ig = 0
    Dg = 0

    draught     = 0

    # A11 = 0.3 * rho * pv * 0.01 
    A11 = mass*0.01 # "Added mass in surge motion is relatively low in comparison to ship mass"
    #https://www.tandfonline.com/doi/full/10.1080/17445302.2019.1615705#:~:text=The%20added%20mass%20is%20referred,considerable%20influence%20on%20ship%20manoeuvrability.
    # A13 = 0
    # A31 = 0
    # A53 = 0.1 * rho * pv * pl
    # A35 = 0.05* rho * pv * pl
    A33 = 0.3 * rho * pv
    # A51 = 0.01* rho * pv * pl
    A55 = 0.07* rho * pv * pl * pl
    B11 = 2.0* rho * pv * math.sqrt(gg/pl)*0.17
    # B53 = 0.05* rho * pv * pl * math.sqrt(gg/pl)
    # B35 = 0.4* rho * pv * pl * math.sqrt(gg/pl)
    B33 = 2.0* rho * pv * math.sqrt(gg/pl)
    B55 = 0.4* rho * pv * pl * pl * math.sqrt(gg/pl)
    # C11 = 0
    # C53 = 0
    # C35 = 0
    C33 = 11014
    C55 = 6233

    def __init__(self,x,z,P,I,D,power_type,prop_type):
        self.posx = x
        self.posz = z
        self.accx = 0
        self.accz = 0
        self.velx = 0
        self.velz = 0

        self.power_type = power_type
        self.prop_type  = prop_type
        self.Pg = P
        self.Ig = I
        self.Dg = D

    def calcAcc(self,ww,t):
        tFx=0
        tFz=0
        tMy=0
        self.IMU_accx = 0
        self.IMU_accz = 0
        self.IMU_roty = 0

        sum_draught = 0
        res = 5
        average_wave_height = 0
        cos = math.cos(self.angle)
        sin = math.sin(self.angle)
        for i in range(res):
            dl = self.pl/2/(res-1)*cos
            dh = self.pl/2/(res-1)*sin
            if(i==0):
                waveheight = ww.get(self.posx,t)
                sum_draught += min(max(waveheight-(self.posz-self.ph),0), self.ph)
                average_wave_height += waveheight
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
                average_wave_height += waveheight1
                average_wave_height += waveheight2
                if(i==res-1): 
                    self.wave_incline = math.atan2((waveheight1-waveheight2),self.pl*cos)

        self.draught = sum_draught/(2*res-1)
        average_wave_height = average_wave_height/(2*res-1)

        Fn = math.fabs(self.velx)/math.sqrt(9.8*self.pl)

        pvelx = 0
        pvelz = 0

        #Hydrodynamic Force
        for w in ww.waves:
            omega_e = w.omega + w.omega**2 * self.velx / gg
            # k= 0.1
            Fx = w.calcF(self.posx,t,Fn,omega_e,axis=1)*self.draught/self.ph
            Fz = w.calcF(self.posx,t,Fn,omega_e,axis=3)*self.draught/self.ph
            My = w.calcF(self.posx,t,Fn,omega_e,axis=5)*self.draught/self.ph
            tFx += Fx
            tFz += Fz
            tMy += My
            self.IMU_accx += (Fz*sin+Fx*cos)/self.mass
            self.IMU_accz += (Fz*cos-Fx*sin)/self.mass
            # particle_velx += w.getVelx(self.posx-self.pl/2,t)
            coeff = w.omega * w.amp * math.exp(w.k*-self.ph/2)
            pvelx += coeff*math.sin(w.omega*t-w.k*self.posx)
            pvelz += coeff*math.cos(w.omega*t-w.k*self.posx)

        #Thruster Force 
        # hz = pwm_curve(pwm_us)
        n = self.prop_rev
        U_rel = (self.velx-pvelx)
        if(n==0):
            thr_f = 0
        elif(n>0):
            # thr_f = -0.286*U_rel*hz + 0.01375*hz**2
            thr_f = rho*self.prop_d**4*n**2*(-0.362*U_rel/(n*self.prop_d) + 0.445)
        else:
            # thr_f = -0.246*U_rel*hz - 0.01185*hz**2
            thr_f = rho*self.prop_d**4*n**2*(-0.290*U_rel/(-n*self.prop_d) - 0.357)
        thr_f *= 4
        tFx += thr_f * cos
        tFz += thr_f * sin
        tMy += thr_f * self.ph/2
        self.IMU_accx += thr_f/self.mass - gg*sin
        self.IMU_accz += -gg*cos

        z_t = self.posz - 0.185 - average_wave_height
        heave = (tFz - self.B33*(self.velz + dt/2*self.accz) - 1*self.C33*(z_t+dt*self.velz+(0.5-beta)*dt*dt*self.accz)) / ((self.mass+self.A33) + dt/2*self.B33 + beta*dt*dt*self.C33)
        surge = (tFx - self.B11*(self.velx + dt/2*self.accx)) / ((self.mass+self.A11) + dt/2*self.B11)
        pitch = (tMy - self.B55*(self.anglevel + dt/2*self.angleacc) - self.C55*(self.angle+dt*self.anglevel+(0.5-beta)*dt*dt*self.angleacc)) / ((self.im+self.A55) + dt/2*self.B55 + beta*dt*dt*self.C55)

        self.accx = surge
        self.accz = heave
        self.angleacc = pitch


    def update(self,ww,t):
        #PID Controller here
        target_vel = 0
        target_acc = 0
        # #update target rev every 0.35s 
        # if(t%(0.35)<dt and self.control_delay < dt):
        #     vel_err = target_vel - self.velx
        #     acc_err = target_acc - self.accx
        #     self.target_rev = self.Pg*vel_err +self.Dg*acc_err
        #     if(self.target_rev>54.8): self.target_rev = 54.8
        #     if(self.target_rev<-54.1): self.target_rev = -54.1
        #     rev_err = self.target_rev - self.prop_rev

        #     #reverse or not
        #     if(self.power_type == 'engine'):
        #         if(self.target_rev*self.prop_rev < 0):
        #             self.control_delay      = 0.35*4
        #             self.transition_delay   = 0.35*0
        #         else: 
        #             self.control_delay      = 0.35*1 
        #             self.transition_delay   = 0.35*1*((abs(rev_err)+0.01)/54.8)
        #     elif(self.power_type == 'motor'):
        #         self.control_delay      = 0.35*0 
        #         self.transition_delay   = 0.35*1*((abs(rev_err)+0.5)/54.8)
        #     # self.delayed_rate = rev_err/self.transition_delay
        #     self.control_rate = rev_err/self.control_delay





        # #change prop rev
        # self.prop_rev += self.control_rate * dt
        # self.prop_rev = 
        # if(self.prop_rev>54.8): self.prop_rev = 54.8
        # if(self.prop_rev<-54.1): self.prop_rev = -54.1
        # #decrement remaining delay
        # if(self.control_delay > dt): 
        #     self.control_delay -= dt
        # else: 
        #     self.control_delay = 0
        #     #keep control rate aligned with delayed rate unless during control delay
        #     # self.control_rate = self.delayed_rate
        #     # print(self.control_rate)

        #update target rev every 0.35s 
        if(t%(0.35)<dt and self.control_delay < dt):
            vel_err = target_vel - self.velx
            acc_err = target_acc - self.accx
            self.target_rev = self.Pg*vel_err +self.Dg*acc_err
            if(self.target_rev>54.8): self.target_rev = 54.8
            if(self.target_rev<-54.1): self.target_rev = -54.1
            self.rev_err = self.target_rev - self.prop_rev

            #reverse or not
            if(self.power_type == 'engine'):
                if(self.target_rev*self.prop_rev < 0):
                    self.control_delay      = 0.35*4
                else: 
                    self.control_delay      = 0.35*1*((abs(self.rev_err)+0.01)/54.8)
            elif(self.power_type == 'motor'):
                self.control_delay      = 0.35*1*((abs(self.rev_err)+0.5)/54.8)
            self.time_constant = self.control_delay
            self.control_rate = self.rev_err/self.control_delay
            self.control_start = t
            self.prev_rev = self.prop_rev

        #change prop rev
        self.prop_rev = self.prev_rev+self.rev_err*(1-math.exp(-(t-self.control_start)/self.time_constant))
        if(self.prop_rev>54.8): self.prop_rev = 54.8
        if(self.prop_rev<-54.1): self.prop_rev = -54.1
        #decrement remaining delay
        if(self.control_delay > dt): 
            self.control_delay -= dt
        else:
            self.control_delay = 0

        self.calcAcc(ww,t)
        self.velx += self.accx*dt
        self.velz += self.accz*dt
        self.anglevel += self.angleacc*dt
        self.posx += self.velx*dt
        self.posz += self.velz*dt
        self.angle += self.anglevel*dt
