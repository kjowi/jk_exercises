#Oscillations of a buoyant object
# Water column: 100m
# surface density: 1025 kg m^-3
# stability frequency squared : N^2=10^-4 s^-2
#
# object with density: 1025.5kg m^-3
# release depth of object: 80 m
#
# Question: predict motion path

# Governing equations
# dw_obj/dt=-g(rho_obj-rho_amb)/rho_amb
# dz/dt=w_obj
import numpy as np
from matplotlib import pyplot as plt
#0. Model constanst
#=====================

delt = 1.
time = np.arange(0,30*60,delt)
rho_surf = 1025.0
N_sq = 10**-4
grav = 9.81
w_obj = np.zeros(len(time))
w_obj[0] = 0.0
rho_obj = 1025.5
pos_obj = np.zeros(len(time))
pos_obj[0] = -80.
#1. Density of the water column
#=================================
def density(depth_obj):
    return rho_surf*(1+N_sq*np.abs(depth_obj)/grav)


#2. Speed of object at given time
#=================================
def vel_obj(wobj,timestep,posobj):
    return wobj-grav*delt*(rho_obj-density(pos_obj[timestep-1]))/rho_obj

#3. Position of object at given time
#====================================

for k in range(1,len(time)):
    w_obj[k] = vel_obj(w_obj[k-1],k,pos_obj[k-1])
    pos_obj[k] = w_obj[k]*delt+pos_obj[k-1]
plt.plot(time/60,pos_obj,color='red')
pos_obj[0] = -50.
for l in range(1,len(time)):
    w_obj[l] = vel_obj(w_obj[l-1],l,pos_obj[l-1])
    pos_obj[l] = w_obj[l]*delt+pos_obj[l-1]
plt.plot(time/60,pos_obj,color='blue')
for m in range(len(time)):
    pos_obj[m] = -1*(1025.5/rho_surf-1)*grav/N_sq

#4. Picture
#=============
plt.plot(time/60,pos_obj,color='green')
plt.axis([0, 30, -100., 0.0])
plt.xlabel('Time (minutes)')
plt.ylabel('Depth (m)')
plt.suptitle('Ex 3: Oscillations of a buoyant object')
plt.show()
