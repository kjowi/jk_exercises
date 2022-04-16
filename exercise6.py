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
from matplotlib.patches import Polygon
# 1. Model Setup
# ==========================================
domain_x = 1010.0
duration = 200.0
delx = 10.
delz = 10.
delt = 0.1
nx = int(domain_x/delx) + 2
nt = int(duration/delt)
waterdepth = np.zeros((nx,nt))
xpos = np.arange(0,domain_x,delx)
waterdepth[0,:] = -10.0
waterdepth[-1,:] = -10.0
filter_s = 0.05
scenarios = ['C','D']
help_anim = 10
grav = 9.81
mindepth = np.full(nx,0.1)


if not os.path.exists('ex6'):
    os.makedirs('ex6')
cmd_rm = 'rm -f ex6/*.png'
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
#scenario C
waterdepth_C = np.zeros(nx)
bath_C = np.zeros(nx)
waterdepth_C[0] = 0.0
waterdepth_C[-1] = 0.0
bath_C[0] = 0.0
bath_C[-1] = 0.0
for sc_c1 in range(1,52):
    waterdepth_C[sc_c1] = 10.0-(10.5/500*(sc_c1-1)*10)
for sc_c2 in range(52,nx):
    waterdepth_C[sc_c2] = -11.0+(10.5/500*(sc_c2-1)*10)
for bth_c1 in range(1,52):
    bath_C[bth_c1] = (10.5/500*(bth_c1-1)*10)
for bth_c2 in range(52,103):
    bath_C[bth_c2] = 21.0-(10.5/500*(bth_c2-1)*10)
#scenario D
waterdepth_D = np.zeros(nx)
bath_D = np.zeros(nx)
waterdepth_D[0] = 0.0
waterdepth_D[-1] = 0.0
slope0=10.0
bath_D[0] = slope0
bath_D[-1] = 0.0

for bth_d in range(1,nx):
    bath_D[bth_d] = 10+slope0-((slope0)*((bth_d-1)*10)/(1010.0))
print(bath_D[10])
bath_D[45:55] = 14.0
waterdepth_D[1:nx]=-1*bath_D[1:nx]


# 3. code
# ===============================================================
# 3.1 functions
# ----------------------------------------------------------------

def shapiro (elevsin, filt,waterin):
    elev_shap = np.zeros(nx)
    for i1 in range(1,nx-1):
        if waterin[i1] == 1:
            elev_shap[i1] = (1-0.5*(waterin[i1-1]+waterin[i1+1])*filt)*elevsin[i1]+\
                            0.5*filt*(waterin[i1-1]*elevsin[i1-1]+waterin[i1+1]*elevsin[i1+1])

        else:
            elev_shap[i1]=elevsin[i1]
    elev_shap[0] = (1-0.5*filt)*elevsin[0]+0.5*filt*elevsin[1]
    elev_shap[-1] = (1-0.5*filt)*elevsin[-1]+0.5*filt*elevsin[-2]

    return elev_shap
def velocity(uvel,elev,scenario,dw_vel):
    velocitynew=np.zeros(len(elev))
    for k in range(1,len(uvel)-1):
        grad = grav*(elev[k+1]-elev[k])/delx
        if (dw_vel[k] == 1):
            if ((dw_vel[k+1] == 1) or grad <= 0):
                velocitynew[k] = uvel[k] - delt*grad
        else:
            if ((dw_vel[k+1] == 1) and grad >= 0):
                velocitynew[k] = uvel[k] - delt*grad

    if (scenario != 'D'):
        velocitynew[0] = 0.0
        velocitynew[-1] = 0.0
    else:
        velocitynew[0] = 0.0
        velocitynew[-1] = 0.0
        # velocitynew[0] = velocitynew[1]
        # velocitynew[-1] = velocitynew[-2]

    return velocitynew
def elevation(depth,vel,el,dw_el):
    elevnew=np.zeros(len(vel))
    uvel1 = uvel2 = 0.0
    for m in range(1,len(vel)-1):
        if vel[m] >= 0:
            uvel1 = vel[m]*depth[m]
        if vel[m] <0:
            uvel1 = vel[m]*depth[m+1]
        if vel[m-1] >= 0:
            uvel2 = vel[m-1]*depth[m-1]
        if vel[m-1] <0:
            uvel2 = vel[m-1]*depth[m]
        term2 = delt*(uvel1-uvel2)/delx
        elevnew[m]=el[m]-term2
    return elevnew
# 3.1 Main code
# --------------------------------------------------------------------------------
for s in scenarios:
    wave_uvel = np.zeros((nx,nt))
    wave_elev = np.zeros((nx,nt))

    #landsea, with 1 being wet, 0 being dry
    if s == 'A':
        wave_elev[:,0] = wave_el_A
        wave_elev[:,0] = -1*np.minimum(waterdepth[:,0],np.zeros(nx))
        water = np.full(nx,1.0)
        water[0] = 0.0
        water[-1] = 0.0
    if s == 'B':
        water = np.full(nx,1.0)
        water[0] = 0.0
        water[-1] = 0.0
    if s == 'C':
        wave_elev[:,0] = -1*np.minimum(waterdepth_C[:],np.zeros(nx))
        wave_elev[0:20,0] = 1.0
        waterdepth[:,0] = waterdepth_C+wave_elev[:,0]
        water = np.where(waterdepth_C<=0,0,1)
        water[0] = 0.0
        water[-1] = 0.0

    if s == 'D':
        print(waterdepth_D[10])
        wave_elev[:,0] = -1*np.minimum(waterdepth_D[:],np.zeros(nx))
        wave_elev[0:20,0] = wave_elev[0:20,0]+1.0
        waterdepth[:,0] = waterdepth_D+wave_elev[:,0]

        water = np.where(waterdepth_D<=mindepth,0,1)
        water[0] = 0.0
        water[-1] = 0.0
        print(water[45])

    for t in range(0,nt-1):

        # scenario B
        if s == 'B':
            wave_elev[51,t] = wave_el_B[t]

        wave_uvel[:,t+1] = velocity(wave_uvel[:,t],wave_elev[:,t],s,water)
        wave_elev[:,t+1] = elevation(waterdepth[:,t],wave_uvel[:,t+1],wave_elev[:,t],water)
        #shapiro has to be applied every timestep
        whave_sh = np.zeros(nx)
        wave_sh = shapiro(wave_elev[:,t+1],filter_s,water[:])
        wave_elev[:,t+1] = wave_sh[:]
        if s == 'A' or s == 'B':
            waterdepth[:,t+1]=delz+wave_elev[:,t+1]
        if s == 'C':
            waterdepth[:,t+1]=waterdepth_C[:]+wave_elev[:,t+1]
        if s == 'D':
            waterdepth[:,t+1]=waterdepth_D[:]+wave_elev[:,t+1]
        water[:]=np.where(waterdepth[:,t+1]<mindepth[:],0,1)
        print(water[45])
        waterdepth[:,t+1] = np.where((waterdepth[:,t+1]<mindepth[:]),mindepth[:],waterdepth[:,t+1])





    # 4. plotting
    # ==========================================================================
    help = help_anim
    fignum = 0

    for tplot in range(nt):
        if help == help_anim:
            if ((s == 'A') or (s == 'B')):
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
                plt.savefig('ex6/'+"{0:03d}".format(fignum)+s+'elevation.png')
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
                plt.savefig('ex6/'+"{0:03d}".format(fignum)+s+'velocity.png')
                help = 1
                plt.close()
            if ((s == 'C') or (s == 'D')):
                fignum = fignum+1
                plt.rcParams["figure.figsize"] = [9, 1.5]
                fig1,ax1 = plt.subplots()
                if s == 'C':
                    test = bath_C[1:nx-1]+waterdepth[1:nx-1,tplot]
                else:
                    test = bath_D[1:nx-1]+waterdepth[1:nx-1,tplot]
                plt.plot(xpos, test, color='r')
                if s == 'C':
                    polygon1 = Polygon([(0,0), (500,10.5), (1000,0),],fill=True,alpha=1.0)
                    ax1.add_patch(polygon1)
                if s == 'D':
                    polygon1 = Polygon([(0,0),(0,10+slope0),(430,10.0+slope0-slope0/1010*430),(440,14),(530,14),(540,10.+slope0-slope0/1010*540),(1010,10.),],fill=True,alpha=1.0)
                    ax1.add_patch(polygon1)
                plt.xlabel('position Channel (m)')
                plt.ylabel('$waterdepth$ (m)')
                plt.suptitle(str(delt*tplot)+' s')
                if s == 'C':
                    plt.ylim(9.5,11.5,.0)
                else:
                    plt.ylim(10,10+slope0+1,.0)
                plt.xlim(0,1000)
                plt.savefig('ex6/'+"{0:03d}".format(fignum)+s+'waterdepth.png')
                plt.close()
                help = 0
        help = help +1
    w_dir = os.getcwd()
    ex_dir = os.path.join(w_dir,'ex6')
    os.chdir(ex_dir)
    if ((s=='A') or (s=='B')):
        cmd_elanim = 'convert -delay 40 '+ex_dir+'/'+'*'+s+'elevation.png '+ex_dir+'/elevation'+s+'.gif'
        cmd_velanim = 'convert -delay 40 '+ex_dir+'/'+'*'+s+'velocity.png '+ex_dir+'/velocity'+s+'.gif'
        os.system(cmd_elanim)
        os.system(cmd_velanim)
    if ((s == 'C') or (s=='D')):
        cmd_islanim = 'convert -delay 40 '+ex_dir+'/'+'*'+s+'waterdepth.png '+ex_dir+'/waterdepth'+s+'.gif'
        os.system(cmd_islanim)

    os.chdir(w_dir)
    cmd_rm = 'rm -f ex6/*.png'
    os.system(cmd_rm)
