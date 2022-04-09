#Excercise 2: Wave interference
#The author asks to make a bar animation in scilab, however, I choose to make
#matplotlib files and animate them after

# interference of two waves of the same amplitude A0=1m
# The resultant wave can be described as:
#     A(x,t)=A0[sin[2pi(x/l1-t/T1)]+sin(2pi(x/l2-t/T2))]
#with T1 and T2 the wave periods (s) and l1 and l2 the wave lengths (m)

import numpy as np
import math
from matplotlib import pyplot as plt
import os

plt.rcParams["figure.figsize"] = [7.50, 3.50]
plt.rcParams["figure.autolayout"] = True

def waveinterference(xin,tin,l2_in,T2_in):
    A0 = 1.
    l1 = 100.
    T1 = 60.
    return A0*(np.sin(2*math.pi*(xin/l1-tin/T1))+\
    np.sin(2*math.pi*(xin/l2_in-tin/T2_in)))
def wave1(xin,tin):
    A0 = 1.
    l1 = 100.
    T1 = 60.
    return A0*(np.sin(2*math.pi*(xin/l1-tin/T1)))
def wave2(xin,tin,l2_in,T2_in):
    A0 = 1.
    return A0*(np.sin(2*math.pi*(xin/l2_in-tin/T2_in)))

x = np.arange(0,1000,2)
t = np.arange(0,180,1)

T2 = [50.,60.,50.,-60.,-30.,-30.]
l2 = [100.,90.,90.,100.,50.,95.]

if not os.path.exists('ex2'):
    os.makedirs('ex2')
for j in range(len(T2)):
    if not os.path.exists('ex2/scenario'+str(j+1)):
        os.makedirs('ex2/scenario'+str(j+1))
    for i in range(len(t)):
        plt.axis([0, 1000, -2.0, 2.0])
        plt.plot(x, waveinterference(x,t[i],l2[j],T2[j]), color='purple',label='wave 1 +2')
        plt.plot(x, wave1(x,t[i]), color='blue', label='wave 1')
        plt.plot(x, wave2(x,t[i],l2[j],T2[j]), color='red', label='wave 2')
        plt.legend(loc='upper right')
        plt.suptitle(t[i])
        plt.savefig('ex2/scenario'+str(j+1)+'/'+str(i)+'.png')
        plt.close()
