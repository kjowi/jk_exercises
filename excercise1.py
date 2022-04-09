#excercise 1: wine thief
#===============================================================================
# Decaying problem: dC/dt = -k*C
#                   C(0) = 100%
#                   k = 0.0001 s^-1
#Analytic solution: C(t)=C(0)*exp(-k*t)
#Question: solve the equation with explicit, implicit and semi-implicit scheme
#          Analyse impact of chosen timestep
import numpy as np
from matplotlib import pyplot as plt
import math
C0 = 100.0
k = 0.0001

# delt = 2000000
def analytic(endtime,delt):
    time=np.arange(0,endtime,delt)
    op = time
    for i in range(len(op)):
        op[i]=C0*math.exp(-1*k*time[i])
    # plt.plot(time/3600, c, color='red',label='analytic')
    return op
def explicit(end,delt):
    c = np.zeros(int(end/delt))
    if k < 1/delt:
        op = c
        op[0] = C0
        for i in range(1,len(c)):
            op[i]=(1-k*delt)*op[i-1]
        # plt.plot(time/3600, c, color='blue',label='explicit')
        return op

    else:
        print("The stability condition (k<1/delt) was not fulfilled")
        return c

def implicit(end,delt):
    c = np.zeros(int(end/delt))
    c[0]=C0
    for i in range(1,len(c)):
            c[i]=c[i-1]/(1+k*delt)
    # plt.plot(time/3600, c, color='green',label='implicit')
    return c

def hybrid(end,delt,alfa):
    c = np.zeros(int(end/delt))
    c[0]=C0
    for i in range(1,len(c)):
            c[i]=c[i-1]*(1-(1-alfa)*k*delt)/(1+alfa*k*delt)
    # plt.plot(time/3600, c, color='purple',label='hybrid')
    return c


plt.rcParams["figure.figsize"] = [7.5, 7.5]
fig1, axs = plt.subplots(3)
fig1.suptitle("Impact of timestep on discretized results")

end = 3600*15
delt = 3600
ana_sol = analytic(end,delt)
expl_sol = explicit(end,delt)
impl_sol = implicit(end,delt)
time = np.arange(0,end,delt)
axs[0].plot(time/3600,ana_sol,color='red',label='analytic')
axs[0].plot(time/3600,expl_sol,color='blue',label='explicit')
axs[0].plot(time/3600,impl_sol,color='green',label='implicit')
axs[0].set_title('\u0394t = 1h')
axs[0].legend()


delt = 1800
analytic_36 = analytic(end,delt)
explicit_36 = explicit(end, delt)
implicit_36 = implicit(end, delt)
time36 = np.arange(0,end,delt)
axs[1].plot(time36/3600,analytic_36,color='red',label='analytic')
axs[1].plot(time36/3600,explicit_36,color='blue',label='explicit')
axs[1].plot(time36/3600,implicit_36,color='green',label='implicit')
axs[1].set_title('\u0394t = 0.5h')
axs[1].legend()


delt = 1
analytic_1 = analytic(end,delt)
explicit_1 = explicit(end,delt)
implicit_1 = implicit(end,delt)
time1 = np.arange(0,end,delt)
axs[2].plot(time1/3600,analytic_1,color='red',label='analytic')
axs[2].plot(time1/3600,explicit_1,color='blue',label='explicit')
axs[2].plot(time1/3600,implicit_1,color='green',label='implicit')
axs[2].set_title('\u0394t = 1 s')
axs[2].legend()


# Second figure: impact of alfa on hybrid scheme preformance
plt.rcParams["figure.figsize"] = [15, 7]
fig2, axs2 = plt.subplots(3,2)
fig2.suptitle("Impact of alfa on hybrid schemes")
delt = 3600
hybrid_sol = hybrid(end,delt,0.75)
analytic_sol = analytic(end,delt)
time = np.arange(0,end,delt)
axs2[0, 0].plot(time/3600,analytic_sol,color='red',label='analytic')
axs2[0, 0].plot(time/3600,hybrid_sol,color='blue',label='hybrid')
axs2[0, 0].set_title('\u0394t = 1 h and \u03B1 = 0.75')
axs2[0, 0].legend()

delt = 3600
hybrid_sol = hybrid(end,delt,0.5)
analytic_sol = analytic(end,delt)
time = np.arange(0,end,delt)
axs2[1, 0].plot(time/3600,analytic_sol,color='red',label='analytic')
axs2[1, 0].plot(time/3600,hybrid_sol,color='blue',label='hybrid (semi-implicit)')
axs2[1, 0].set_title('\u0394t = 1 h and \u03B1 = 0.5')
axs2[1, 0].legend()

delt = 3600
hybrid_sol = hybrid(end,delt,0.25)
analytic_sol = analytic(end,delt)
time = np.arange(0,end,delt)
axs2[2, 0].plot(time/3600,analytic_sol,color='red',label='analytic')
axs2[2, 0].plot(time/3600,hybrid_sol,color='blue',label='hybrid')
axs2[2, 0].set_title('\u0394t = 1 h and \u03B1 = 0.25')
axs2[2, 0].legend()

delt = 1800
hybrid_sol = hybrid(end,delt,0.75)
analytic_sol = analytic(end,delt)
time = np.arange(0,end,delt)
axs2[0, 1].plot(time/3600,analytic_sol,color='red',label='analytic')
axs2[0, 1].plot(time/3600,hybrid_sol,color='blue',label='hybrid')
axs2[0, 1].set_title('\u0394t = 0.5h and \u03B1 = 0.75')
axs2[0, 1].legend()

delt = 1800
hybrid_sol = hybrid(end,delt,0.5)
analytic_sol = analytic(end,delt)
time = np.arange(0,end,delt)
axs2[1, 1].plot(time/3600,analytic_sol,color='red',label='analytic')
axs2[1, 1].plot(time/3600,hybrid_sol,color='blue',label='hybrid (semi-implicit)')
axs2[1, 1].set_title('\u0394t = 0.5h and \u03B1 = 0.5')
axs2[1, 1].legend()

delt = 1800
hybrid_sol = hybrid(end,delt,0.25)
analytic_sol = analytic(end,delt)
time = np.arange(0,end,delt)
axs2[2, 1].plot(time/3600,analytic_sol,color='red',label='analytic')
axs2[2, 1].plot(time/3600,hybrid_sol,color='blue',label='hybrid')
axs2[2, 1].set_title('\u0394t = 0.5h and \u03B1 = 0.25')
axs2[2, 1].legend()


for ax in axs.flat:
    ax.set(xlabel='time (h)', ylabel='Concentration (%)')

for ax in axs2.flat:
    ax.set(xlabel='time (h)', ylabel='Concentration (%)')

# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axs.flat:
    ax.label_outer()

for ax in axs2.flat:
    ax.label_outer()
plt.show()
