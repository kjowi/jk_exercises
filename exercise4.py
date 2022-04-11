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

# 0. Settings
# ==============================================================

import numpy as np
from matplotlib import pyplot as plt
import math
r_tank = 20000.0
omega = -0.727*10**-5
f = 2* omega
simulationlength = 20*24*3600
modes = ['rotational']#here you can define in which mode you want to run the code,
                    # possibilities are: explicit, semi-implicit,rotational

for m in modes:
    if m == 'explicit':
        # 1. Explicit scheme
        # ==============================================================
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

# 2. Semi-implicit scheme
# ==============================================================
    if m == 'semi-implicit':
        def posvel_si(delt):
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
            for i in range(1,len(time)):
                alfa = delt*f
                beta = 0.25*alfa**2
                time[i] = time[i-1]+delt
                uvel[i] = ((1-beta)*uvel[i-1]+alfa*vvel[i-1])/(1+beta)
                vvel[i] = ((1-beta)*vvel[i-1]-alfa*uvel[i-1])/(1+beta)
                x_pos[i] = 0.5*(uvel[i-1]+uvel[i])*delt + x_pos[i-1]
                y_pos[i] = 0.5*(vvel[i-1]+vvel[i])*delt + y_pos[i-1]
            return x_pos,y_pos,time,uvel,vvel

        deltimes =[1,36,360,3600]
        clrs = ['b','g','r','y']
        circle = plt.Circle((0, 0), r_tank, color='black', fill=False)
        plt.rcParams["figure.figsize"] = [10.0, 10.0]
        fig1, ax = plt.subplots()
        plt.suptitle('Semi-implicit scheme')
        ax = plt.gca()
        ax.cla() # clear things for fresh plot
        plt.xlabel('X (m)')
        plt.ylabel('Y (m)')

        fig2, axs2 = plt.subplots()
        plt.rcParams["figure.figsize"] = [7.5, 3.5]
        plt.suptitle('Semi-implicit scheme')
        ax2 = plt.gca()
        for j in range(len(deltimes)):
            timesteps = int(simulationlength/deltimes[j])

            xpos,ypos,timemodel,uvelm,vvelm = posvel_si(deltimes[j])
            ax.set_xlim((-r_tank,r_tank ))
            ax.set_ylim((-r_tank, r_tank))
            ax.plot(xpos, ypos, color=clrs[j],label = '$\Delta$ t ='+str(deltimes[j])+'s')
            ax.add_patch(circle)
            ax2.plot(timemodel/3600, np.sqrt(uvelm**2+vvelm**2), color=clrs[j],label ='$\Delta$ t ='+str(deltimes[j])+'s')
        ax.legend()
        realspeed = np.full(timesteps,np.sqrt(uvelm[0]**2+vvelm[0]**2))
        plt.legend()
        plt.ylim(0.0705,0.0709)
        ax2.plot(timemodel/3600, realspeed, color='black',label='analytic result')
        plt.xlabel('Time (h)')
        plt.ylabel('Parcel speed (m s $^{-1}$)')
        plt.show()

# 3. Rotational approach, which is basicly an explicit scheme
# ==============================================================
    if m == 'rotational':
        def posvel_rot(delt,angle):
            for i in range(1,len(time)):
                time[i] = time[i-1]+delt
                uvel[i] = np.cos(angle) * uvel[i-1] + np.sin(angle)*vvel[i-1]
                vvel[i] = np.cos(angle) * vvel[i-1] - np.sin(angle)*uvel[i-1]
                x_pos[i] = uvel[i-1]*delt + x_pos[i-1]
                y_pos[i] = vvel[i-1]*delt + y_pos[i-1]
            return x_pos,y_pos,time,uvel,vvel

        deltimes =[1,36,360,3600*5]
        clrs = ['b','g','r','y']
        circle = plt.Circle((0, 0), r_tank, color='black', fill=False)
        plt.rcParams["figure.figsize"] = [10.0, 10.0]

        fig1, ax = plt.subplots()
        plt.suptitle( 'Rotational scheme')
        ax = plt.gca()
        ax.cla() # clear things for fresh plot
        plt.xlabel('X (m)')
        plt.ylabel('Y (m)')

        fig2, axs2 = plt.subplots()
        plt.rcParams["figure.figsize"] = [7.5, 3.5]
        plt.suptitle( 'Rotational scheme')# with $\delta t* f =$ '+str(delt*f))
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
            ang = deltimes[j]*f #np.arcsin(0.5*delt*f)
            xpos,ypos,timemodel,uvelm,vvelm = posvel_rot(deltimes[j],ang)
            # change default range so that new circles will work
            ax.set_xlim((-r_tank,r_tank ))
            ax.set_ylim((-r_tank, r_tank))
            # key data point that we are encircling

            ax.plot(xpos, ypos, color=clrs[j],label = '$\Delta$ t ='+str(deltimes[j])+'s and $\delta t* f =$ '+str(deltimes[j]*f))


            ax.add_patch(circle)
            # plt.show()


            ax2.plot(timemodel/3600, np.sqrt(uvelm**2+vvelm**2), color=clrs[j],label ='$\Delta$ t ='+str(deltimes[j])+'s and $\delta t* f =$ '+str(deltimes[j]*f))
        ax.legend()
        realspeed = np.full(timesteps,np.sqrt(uvelm[0]**2+vvelm[0]**2))
        plt.legend()
        plt.ylim(0.07,0.075)
        ax2.plot(timemodel/3600, realspeed, color='black',label='analytic result')
        plt.xlabel('Time (h)')
        plt.ylabel('Parcel speed (m s $^{-1}$)')
        plt.show()
