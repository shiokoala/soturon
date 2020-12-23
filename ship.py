import math
import simulation_settings as ss
import constants

gg = constants.gg
rho = constants.rho
beta = constants.beta
dt = constants.dt
# scaling = 25/3.11
scaling = 1

def clamp(n, minn, maxn):
    return max(min(maxn, n), minn)

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


class Ship2:
    posx = 0
    posz = 0
    velx = 0
    velz = 0
    accx = 0
    accz = 0

    IMU_en = False
    IMU_accx = 0
    IMU_accz = 0
    IMU_roty = 0
    IMU_velx = 0
    
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
    prop_d          = 76/1000 * scaling
    prop_pitch      = 94/1000 * scaling
    prop_ratio      = prop_pitch/prop_d
    prop_revmax     = 54.8
    prop_revmin     = -54.1
    prop_anglemax   = math.radians(25)
    prop_anglemin   = math.radians(-20)

    prop_rev        = 0 #hz
    prev_rev        = 0
    P_target_val      = 0
    rev_err         = 0
    thr_f           = 0

    control_delay   = 0
    control_start   = 0
    time_constant   = 0
    delayed_rate    = 0
    transition_delay= 0
    Pg = 0
    Ig = 0
    Dg = 0

    prev_angle      = 0
    prop_angle      = 0
    target_angle    = 0
    angle_err       = 0

    P_start_val     = 0
    P_current_val   = 0
    P_target_val    = 0
    P_target_max    = 0
    P_target_min    = 0
    P_err           = 0
    P_vel_err       = 0
    P_acc_err       = 0
    P_sum_err       = 0

    draught     = 0

    # A11 = mass*0.01 # "Added mass in surge motion is relatively low in comparison to ship mass"
    # #https://www.tandfonline.com/doi/full/10.1080/17445302.2019.1615705#:~:text=The%20added%20mass%20is%20referred,considerable%20influence%20on%20ship%20manoeuvrability.
    # A33 = 0.3 * rho * pv
    # A55 = 0.07* rho * pv * pl * pl
    # B11 = 2.0* rho * pv * math.sqrt(gg/pl)*0.17
    # B33 = 2.0* rho * pv * math.sqrt(gg/pl)
    # B55 = 0.4* rho * pv * pl * pl * math.sqrt(gg/pl)
    # C33 = 11014
    # C55 = 6233

    # A11 = 2 # "Added mass in surge motion is relatively low in comparison to ship mass"
    # #https://www.tandfonline.com/doi/full/10.1080/17445302.2019.1615705#:~:text=The%20added%20mass%20is%20referred,considerable%20influence%20on%20ship%20manoeuvrability.
    # A33 = 350
    # A55 = 100
    # B11 = 2
    # B33 = 300
    # B55 = 100
    # C33 = 11014
    # C55 = 6233


    ### T=3
    A11 = 2 # "Added mass in surge motion is relatively low in comparison to ship mass"
    # #https://www.tandfonline.com/doi/full/10.1080/17445302.2019.1615705#:~:text=The%20added%20mass%20is%20referred,considerable%20influence%20on%20ship%20manoeuvrability.
    A33 = 223
    A55 = 137
    B11 = 0.28
    B33 = 386
    B55 = 18.4
    C33 = 11014
    C55 = 6233

    #cpp Kt characteristics
    a0  =  0.0296
    a1  = -0.2758
    a2  =  0.7244
    a3  = -0.1450
    a4  =  0.0501
    a5  =  0.5259

    def __init__(self,x,z,P,I,D,power_type,prop_type, IMU_en):
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
        self.IMU_en = IMU_en

        #constant rotation speed for cpp
        if(prop_type == 'cpp'):
            self.prop_rev = 50
            self.P_target_max = self.prop_anglemax
            self.P_target_min = self.prop_anglemin
        if(prop_type == 'fpp'):
            self.P_target_max = self.prop_revmax
            self.P_target_min = self.prop_revmin

    def calcAcc(self,ww,t):
        tFx=0
        tFz=0
        tMy=0
        self.IMU_accx = 0
        self.IMU_accz = 0
        self.IMU_roty = 0
        pvelx = 0
        pvelz = 0

        cos = math.cos(self.angle)
        sin = math.sin(self.angle)
        Fn = math.fabs(self.velx)/math.sqrt(9.8*self.pl)

        # calculate average wave height
        sum_draught = 0
        res = 5
        average_wave_height = 0
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

        #Calculate Hydrodynamic Force for each wave
        for w in ww.waves:
            omega_e = w.omega + w.omega**2 * self.velx / gg
            # k= 0.1
            # Fx = w.calcF(self.posx,t,Fn,omega_e,axis=1)*self.draught/self.ph
            # Fz = w.calcF(self.posx,t,Fn,omega_e,axis=3)*self.draught/self.ph
            # My = w.calcF(self.posx,t,Fn,omega_e,axis=5)*self.draught/self.ph
            Fx = w.calcF(self.posx,t,Fn,omega_e,axis=1)
            Fz = w.calcF(self.posx,t,Fn,omega_e,axis=3)
            My = w.calcF(self.posx,t,Fn,omega_e,axis=5)
            tFx += Fx
            tFz += Fz
            tMy += My
            self.IMU_accx += (Fz*sin+Fx*cos)/self.mass
            self.IMU_accz += (Fz*cos-Fx*sin)/self.mass
            # particle_velx += w.getVelx(self.posx-self.pl/2,t)
            coeff = w.omega * w.amp * math.exp(w.k*-self.ph/2)
            pvelx += coeff*math.sin(w.omega*t-w.k*self.posx)
            pvelz += coeff*math.cos(w.omega*t-w.k*self.posx)

        #Calculate Thruster Force
        n = self.prop_rev
        U_rel = (self.velx-pvelx)
        J = U_rel / (n * self.prop_d) if n!=0 else 0
        if(self.prop_type == 'fpp'):
            Kt = math.copysign( math.pi/8 * (self.prop_ratio**2 - J**2) , n)
        elif(self.prop_type == 'cpp'):
            Kt = self.a0 + self.a1*J + self.a2*self.prop_angle + self.a3*J**2 + self.a4*J*self.prop_angle + self.a5*self.prop_angle*abs(self.prop_angle)
        thr_f = rho*Kt*self.prop_d**4*n**2
        thr_f *= 4
###############################
        thr_f = 0
###############################
        self.thr_f = thr_f
        tFx += thr_f * cos
        tFz += thr_f * sin
        tMy += thr_f * self.ph/2
        self.IMU_accx += thr_f/self.mass - gg*sin
        self.IMU_accz += -gg*cos


        #Motion Equation (newmark beta)
        z_t = self.posz - 0.185 - average_wave_height
        # z_t = self.posz - 0.19
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

        #update target every 0.35s 
        if(t%(0.35)<dt and self.control_delay < dt):
            # PD error
            if(self.IMU_en==True):
                self.P_vel_err = target_vel - self.velx
                self.P_acc_err = target_acc - self.IMU_accx
            else:
                self.P_vel_err = target_vel - self.velx
                self.P_acc_err = target_acc - self.accx
            self.P_sum_err += self.P_vel_err
            # set target
            self.P_target_val = (self.Pg*self.P_vel_err + self.Ig*self.P_sum_err +self.Dg*self.P_acc_err)/self.P_target_max # P_target_val -> percent
            self.P_target_val = clamp(self.P_target_val, self.P_target_min/self.P_target_max, 1) #clamp to (-0.8~1.0)
            
            # if(self.power_type=='engine' and self.prop_type=='fpp'):
            #     #rev minimum 18%, cut power if target<9%, 
            #     if(abs(self.P_target_val)<5/self.prop_revmax): # ~9%
            #         self.P_target_val = 0
            #     elif(abs(self.P_target_val)<10/self.prop_revmax): #~18%
            #         self.P_target_val = math.copysign(10/self.prop_revmax,self.P_target_val)
            self.P_err = self.P_target_val - self.P_current_val
            self.P_start_val = self.P_current_val

            #delay
            if(self.power_type == 'engine' and self.prop_type=='fpp'):
                if(self.P_target_val*self.P_current_val < 0):
                    self.control_delay      = 0.35*4
                else: 
                    self.control_delay      = 0.35*1#*((abs(self.rev_err)+0.01)/54.8)
            elif(self.power_type == 'motor' or self.prop_type=='cpp'):
                self.control_delay      = 0.35*1#*((abs(self.rev_err)+0.5)/54.8)
            self.time_constant = self.control_delay
            self.control_start = t

        #change rev or angle
        self.P_current_val = self.P_start_val+self.P_err*(1-math.exp(-(t-self.control_start)/self.time_constant))
        # print((t-self.control_start))
        # self.P_current_val = clamp(self.P_current_val, self.P_target_min/self.P_target_max, 1)
        if(self.prop_type=='cpp'): self.prop_angle = self.P_current_val * self.prop_anglemax
        if(self.prop_type=='fpp'): self.prop_rev = self.P_current_val * self.prop_revmax
        #decrement remaining delay
        self.control_delay = self.control_delay-dt if self.control_delay > dt else 0

        self.calcAcc(ww,t)
        self.velx += self.accx*dt
        self.velz += self.accz*dt
        self.anglevel += self.angleacc*dt
        self.posx += self.velx*dt
        self.posz += self.velz*dt
        self.angle += self.anglevel*dt
        self.IMU_velx += self.IMU_accx*dt