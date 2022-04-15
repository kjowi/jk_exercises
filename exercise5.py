# exercise 5: Long waves in a channel
# ===================================================================
# Aim: simulate the progression of shallow-water surface gravity waves
#      in a channel of uniform water depth
# Scaling: no coriolis force: wave periods are short compared to inertial period
#          no non-linear advection of momentum: wave travels faster than Water
#          no frictional effects at bottom or surface
#          waves don't know any variation in the y-direction
# Setup: length channel = 1000 m
#        constant water depth = 10 m
#        dry cells depth= 0 m
#        Forcing senario's: A. region of 110 m centred in the middle of the domain
#                              has a sealevel of waterdepth +1m at the initial time timestep
#                           B. wave paddle at middle of domain with oscilating sea level:
#                                  period_osc = 20s
#                                  ampl_osc = 1m
# discretization techniques: finite difference and shapiro filter
#                            delt = 0.1
import math
import numpy as np
from matplotlib import pyplot as plt
import os
# 1. Model Setup
# ==========================================
domain_x = 1010.0
duration = 120.0
delx = 10.
delz = 10.
delt = 0.1
nx = int(domain_x/delx) + 2
nt = int(duration/delt)
waterdepth = np.zeros((nx,nt))
xpos = np.arange(0,domain_x,delx)
waterdepth[0,:] = 0.0
waterdepth[-1,:] = 0.0
filter_s = 0.01
scenarios = ['A','B']
help_anim = 10
grav = 9.81
if not os.path.exists('ex5'):
    os.makedirs('ex5')
cmd_rm = 'rm -f ex5/*.png'
os.system(cmd_rm)
#2. Initial conditions
# =========================================
# A. Scenario A
wave_el_A = np.zeros(nx)
wave_el_A[45:55] = 1.0

# Scenario B
wave_el_B = np.zeros(nt)
ampl_osc = 2.0
period_osc = 20.
for t_force in range(nt):
    wave_el_B[t_force] = ampl_osc*np.sin(2*math.pi*(delt*t_force/period_osc))
    print(wave_el_B[t_force])

# 3. code
# ===============================================================
# 3.1 functions
# ----------------------------------------------------------------

def shapiro (elevsin, filt):
    elev_shap = np.zeros(nx)
    for i1 in range(2,nx-2):
        elev_shap[i1] = (1-filt)*elevsin[i1]+0.5*filt*(elevsin[i1-1]+elevsin[i1+1])
    elev_shap[nx-2] = (1-0.5*filt)*elevsin[nx-2]+0.5*filt*(elevsin[i1-1])
    elev_shap[1] = (1-0.5*filt)*elevsin[1]+0.5*filt*(elevsin[i1+1])
    elev_shap[0] = elevsin[0]
    elev_shap[nx-1] = elevsin[nx-1]

    return elev_shap
def velocity(uvel,elev):
    velocitynew=np.zeros(len(elev))
    for k in range(1,len(uvel)-1):
        velocitynew[k] = uvel[k] - delt*grav*(elev[k+1]-elev[k])/delx
    velocitynew[0] = 0.0
    velocitynew[-1] = 0.0
    return velocitynew
def elevation(depth,vel,el):
    elevnew=np.zeros(len(vel))
    uvel1 = uvel2 = 0.0
    for m in range(1,len(vel)):
        if vel[m] >= 0:
            uvel1 = vel[m]*depth[m]
        if vel[m] <0:
            uvel1 = vel[m]*depth[m+1]
        if vel[m-1] >= 0:
            uvel2 = vel[m-1]*depth[m-1]
        if vel[m-1] <0:
            uvel2 = vel[m-1]*depth[m]
        grad = delt*(uvel1-uvel2)/delx
        elevnew[m]=el[m]-grad
    return elevnew
# 3.1 Main code
# --------------------------------------------------------------------------------
for s in scenarios:
    wave_uvel = np.zeros((nx,nt))
    wave_elev = np.zeros((nx,nt))

    if s == 'A':
        wave_elev[:,0] = wave_el_A

    for t in range(0,nt-1):
        # scenario B
        if s == 'B':
            wave_elev[51,t] = wave_el_B[t]

        wave_uvel[:,t+1] = velocity(wave_uvel[:,t],wave_elev[:,t])
        wave_elev[:,t+1] = elevation(waterdepth[:,t],wave_uvel[:,t+1],wave_elev[:,t])
        waterdepth[:,t+1]=delz+wave_elev[:,t+1]

        #shapiro has to be applied every timestep
        whave_sh = np.zeros(nx)
        wave_sh = shapiro(wave_elev[:,t+1],filter_s)
        wave_elev[:,t+1] = wave_sh[:]


    # 4. plotting
    # ==========================================================================
    help = help_anim
    fignum = 0

    for tplot in range(nt):
        if help == help_anim:
            fignum = fignum+1
            plt.rcParams["figure.figsize"] = [9, 1.5]
            fig1,ax1 = plt.subplots()

            plt.plot(xpos, wave_elev[1:nx-1,tplot], color='r')
            plt.xlabel('position Channel (m)')
            plt.ylabel('$\eta$ (m)')
            plt.suptitle(str(delt*tplot)+' s')
            if s == 'A':
                plt.ylim(-1,1,0.5)
            if s == 'B':
                plt.ylim(-4,4,0.5)
            plt.xlim(0,1000)
            plt.savefig('ex5/'+"{0:03d}".format(fignum)+s+'elevation.png')
            plt.close()

            fig2,ax2 = plt.subplots()
            plt.rcParams["figure.figsize"] = [9, 1.5]
            ax2.plot(xpos, wave_uvel[1:nx-1,tplot], color='r')
            plt.xlabel('position Channel (m)')
            plt.ylabel('u (m s$^-1$))')
            plt.suptitle(str(delt*tplot)+' s')
            if s == 'A':
                plt.ylim(-1,1,0.5)
            if s == 'B':
                plt.ylim(-4,4,0.5)
            plt.xlim(0,1000)
            plt.savefig('ex5/'+"{0:03d}".format(fignum)+s+'velocity.png')
            help = 1
            plt.close()
        else:
            help = help +1
    w_dir = os.getcwd()
    ex_dir = os.path.join(w_dir,'ex5')
    os.chdir(ex_dir)
    cmd_elanim = 'convert -delay 40 '+ex_dir+'/'+'*'+s+'elevation.png '+ex_dir+'/elevation'+s+'.gif'
    cmd_velanim = 'convert -delay 40 '+ex_dir+'/'+'*'+s+'velocity.png '+ex_dir+'/velocity'+s+'.gif'
    os.system(cmd_elanim)
    os.system(cmd_velanim)
    os.chdir(w_dir)
    cmd_rm = 'rm -f ex5/*.png'
    os.system(cmd_rm)
