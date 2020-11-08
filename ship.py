import math
import random

gg = 9.8
rho = 997

 t = 0.0
 dt = 0.01
 k = 5
beta = 0.25


width = 1000
height = 800

 
 class Particle{
   PVector pos
   PVector vel
   PVector acc
    angle = 0
    anglevel = 0
    angleacc = 0
   
    rho = 997
    gg = 9.8
   
    beta = 0.25
   
    pl = 3.11
    ph = 0.70
    pw = 0.94
    pv = 0.0878
    mass = pv*100
   
    A11 = 0.3 * rho * pv//
    A13 = 0
    A31 = 0
    A53 = 0.1 * rho * pv * pl
    A35 = 0.05* rho * pv * pl
    A33 = 0.3 * rho * pv
    A51 = 0.01* rho * pv * pl
    A55 = 0.07* rho * pv * pl * pl
    B11 = 2.0* rho * pv * sqrt(gg/pl)
    B53 = 0.05* rho * pv * pl * sqrt(gg/pl)
    B35 = 0.4* rho * pv * pl * sqrt(gg/pl)
    B33 = 2.0* rho * pv * sqrt(gg/pl)
    B55 = 0.4* rho * pv * pl * pl * sqrt(gg/pl)
    C11 = 0
    C53 = 0
    C35 = 0
    C33 = 11014
   
   
   
   
   Particle( x,  y){
     pos = new PVector(x,y)
     vel = new PVector(0,0)
     acc = new PVector(0,0)
   }
   void calcAcc(){
 
     PVector totalForce = new PVector(0.,0.)
      totalMoment = 0
         
     totalForce.y += mass*gg //Gravity
     PVector cP = ww.particles[round(pos.x)]
     
      draught = min(max(pos.y-cP.y,0), ph)
     text(draught,800,450)
     
      nx = cos(angle - PI/2)
      ny = sin(angle - PI/2)
     
     if(draught>0){
       h111 =0
       h1111 = 0
       h112 =0
       h113 =0
       for(Wave w:ww.waves){
          ax = w.omega*w.omega*w.amp*exp(w.k*(height/2-cP.y))*cos(w.omega*t + w.k*pos.x)
          ay = -w.omega*w.omega*w.amp*exp(w.k*(height/2-cP.y))*sin(w.omega*t + w.k*pos.x)
          pd = rho*gg*w.amp*exp(w.k*(height/2-cP.y))*sin(w.omega*t + w.k*pos.x)
         totalForce.x += -pw*pl*pd*nx + A11*ax + A13*ay
         totalForce.y += -pw*pl*pd*ny + A31*ax + A33*ay
         totalMoment += A51*ax + A53*ay
         if(pw*pl*pd*ny>h111) h111 = pw*pl*pd*ny
         if(A11*ax>h112) h112 = A11*ax
         if(A13*ay>h113) h113 = A13*ay
         if(exp(w.k*(height/2-cP.y))>h1111) h1111 = exp(w.k*(height/2-cP.y))
       }
        y_t = -(height/2 - pos.y)
       h1 = (totalForce.y - B33*(vel.y + dt/2*acc.y) - C33*(y_t+dt*vel.y+(0.5-beta)*dt*dt*acc.y))
       h11 = totalForce.y
       h12 = - B33*(vel.y + dt/2*acc.y)
       h13 = - C33*(y_t+dt*vel.y+(0.5-beta)*dt*dt*acc.y)
       h2 = ((mass+A33) + dt/2*B33 + beta*dt*dt*C33)
        heave = (totalForce.y - B33*(vel.y + dt/2*acc.y) - C33*(y_t+dt*vel.y+(0.5-beta)*dt*dt*acc.y)) / ((mass+A33) + dt/2*B33 + beta*dt*dt*C33)
        surge = (totalForce.x - B11*(vel.x + dt/2*acc.x)) / ((mass+A11) + dt/2*B11)
       
       text(heave,800,650)
       
       acc.x = surge
       acc.y = heave
     }else{
       acc.x = 0
       acc.y = 9.8
     }
     
 
     pushMatrix()
     translate(pos.x,pos.y)
     line(0,0,acc.x,acc.y)
     popMatrix()
   }
   
   void update(){
     calcAcc()
     if(keyPressed==true){
       if(keyCode == RIGHT){
         acc.x+=1
       }
       if(keyCode == LEFT){
         acc.x-=1
       }
     }
     
     vel = vel.add(acc.copy().mult(dt))
     pos = pos.add(vel.copy().mult(dt))
     pushMatrix()
     translate(pos.x,pos.y)
     line( pl/2, ph/2,-pl/2, ph/2)
     line(-pl/2, ph/2,-pl/2,-ph/2)
     line(-pl/2,-ph/2, pl/2,-ph/2)
     line( pl/2,-ph/2, pl/2, ph/2)
     popMatrix()
   }
 }


##############################################################################

##############################################################################

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



class Water:
    waves = []
    def __init__(self, num):
        for i in range(num):
            newWave = Wave(random.uniform(0.01,0.1),random.uniform(0.01,1),random.uniform(1,10))
            self.waves.append(newWave)
    def get(self,x,t):
        waveheight = 0
        for w in self.waves:
            waveheight += w.get(x,t)
        return waveheight

