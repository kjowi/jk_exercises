# Exercise 4: The coriolis force in action
#
# Aim: predict the pathway of a non-buoyant fluid parce in
#      a rotating fluid subject to the coriolis force
# Settings:
#     diameter tank: 20 km
#     rotating rate: omega = -0.727*10**-5 s^-1 (clockwise rotation and T= 24h)
     # x0 =0
     # y0= 5km
     # u0=0.5m/s
     # v0=0.5m/s

# Capitals denote static frame of reference, small caps denote rotating frame of reference

# 0. analytical solution
# ==============================================================



# 1. Explicit scheme
# ==============================================================
import numpy as np
from matplotlib import pyplot as plt
import math
r_tank = 20000.0
omega = -0.727*10**-5
f = 2* omega
simulationlength = 20*24*3600

def posvel(delt):
    for i in range(1,len(time)):
        time[i] = time[i-1]+delt
        uvel[i] = f*vvel[i-1]*delt+uvel[i-1]
        vvel[i] = -1*f*uvel[i-1]*delt+vvel[i-1]
        x_pos[i] = uvel[i-1]*delt + x_pos[i-1]
        y_pos[i] = vvel[i-1]*delt + y_pos[i-1]
    return x_pos,y_pos,time,uvel,vvel

deltimes =[1,36,360,3600]
clrs = ['b','g','r','y']
circle = plt.Circle((0, 0), r_tank, color='black', fill=False)
plt.rcParams["figure.figsize"] = [10.0, 10.0]
fig1, ax = plt.subplots()
ax = plt.gca()
ax.cla() # clear things for fresh plot
plt.xlabel('X (m)')
plt.ylabel('Y (m)')

fig2, axs2 = plt.subplots()
plt.rcParams["figure.figsize"] = [7.5, 3.5]
ax2 = plt.gca()
for j in range(len(deltimes)):
    timesteps = int(simulationlength/deltimes[j])
    time = np.zeros(timesteps)
    x_pos = np.zeros(timesteps)
    y_pos = np.zeros(timesteps)
    uvel = np.zeros(timesteps)
    vvel = np.zeros(timesteps)
    x_pos[0] = 0.0
    y_pos[0] = 5000.0
    uvel[0] = 0.05
    vvel[0] = 0.05
    time[0] = 0.0
    xpos,ypos,timemodel,uvelm,vvelm = posvel(deltimes[j])
    # change default range so that new circles will work
    ax.set_xlim((-r_tank,r_tank ))
    ax.set_ylim((-r_tank, r_tank))
    # key data point that we are encircling

    ax.plot(xpos, ypos, color=clrs[j],label = '$\Delta$ t ='+str(deltimes[j])+'s')


    ax.add_patch(circle)
    # plt.show()


    ax2.plot(timemodel/3600, np.sqrt(uvelm**2+vvelm**2), color=clrs[j],label ='$\Delta$ t ='+str(deltimes[j])+'s')
ax.legend()
realspeed = np.full(timesteps,np.sqrt(uvelm[0]**2+vvelm[0]**2))
plt.legend()
plt.ylim(0.07,0.075)
ax2.plot(timemodel/3600, realspeed, color='black',label='analytic result')
plt.xlabel('Time (h)')
plt.ylabel('Parcel speed (m s $^{-1}$)')
plt.show()
