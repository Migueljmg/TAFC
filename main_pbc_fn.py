from solver_pbc_fn import Problem
from scipy.optimize import curve_fit
import random
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os
from time import time


def plot_energy(period):
    b=p1.getallenergy(period)
    plt.plot(range(len(b)),b)
    plt.grid(True)
    plt.show()

def plot_vel_distr(num_iter,period):
    vlist=p1.getvelocities(num_iter-1-(num_iter-1)%period)
    n,bins,patches=plt.hist(vlist, bins=20)
    bin_centers=bins[:-1]+0.5*(bins[1:]-bins[:-1])

    def func(x, maximum_y, x0, sigma):
        return maximum_y*np.exp(-(x-x0)**2/(2*sigma**2))

    from scipy.optimize import curve_fit
    # Edit the starting values and the function to adapt for your need
    popt, pcov = curve_fit(func, bin_centers, n, p0 = [10, 0, 2])

    plt.plot(bin_centers, func(bin_centers, *popt), 'r-')
    plt.grid(True)
    plt.show()

def animate(num_iter,period):
    traj=p1.getalltrajectories(period)

    sheet_figure = plt.figure()
    ax = plt.axes(xlim=(-delta, N*delta), ylim=(-0.1, 0.1))
    point_set, = ax.plot([traj[0][i] for i in range(N)],[0 for i in range(N)], 'ro')
    plt.axvline(x=-delta/2)
    plt.axvline(x=(N-0.5)*delta)
    pt_circ = plt.Circle((4, 4), 1, color='b', fill=False)
    ax.add_artist(pt_circ)

    def run_animation(i):
        point_set.set_data([traj[i][k] for k in range(N)],[0 for k in range(N)])
        return point_set,

    anim = animation.FuncAnimation(sheet_figure, run_animation, frames=(num_iter-1)/period, interval=10)

    plt.show()

#General way of creating 2D arrays
#a = [[0] * ncol for i in range(nrow)]
os.system("rm ./lixo/*.pkl")
random.seed(42)
t0=time()
N=10
wpdt=0.001
num_iter=10000
delta=0.1 #1/10 of the Debye length
nr_images=0 #number of image charges
period=10
init_pos=[delta*i for i in range(N)]
init_vel=[random.uniform(-0.3,0.3) for i in range(N)]
#init_vel[0]=3
#init_vel[1]=-0.7
#init_vel[2]=0.05
#init_vel[3]=0.05
p1 = Problem(N,init_pos,init_vel,wpdt,delta,num_iter,nr_images,period)

p1.start()
t1=time()

print(t1-t0)
#p1.cenas()

#plot_energy(period)
#plot_vel_distr(num_iter,period)
animate(num_iter,period)


"""
#TRAJECTORY
a=p1.getalltrajectories(period)
print(len(a))
plt.plot(range(num_iter-1), [c[0] for c in a])
plt.plot(range(num_iter-1), [c[1] for c in a])
plt.grid(True)
plt.show()
"""


"""
#ENERGY
b=p1.getallenergy(period)
plt.plot(range(len(b)),b)
plt.grid(True)
plt.show()
"""

"""
#MAXWELLIAN
vlist=p1.getvelocities(num_iter-1-(num_iter-1)%period)
n,bins,patches=plt.hist(vlist, bins=20)
bin_centers=bins[:-1]+0.5*(bins[1:]-bins[:-1])

def func(x, maximum_y, x0, sigma):
    return maximum_y*np.exp(-(x-x0)**2/(2*sigma**2))

from scipy.optimize import curve_fit
# Edit the starting values and the function to adapt for your need
popt, pcov = curve_fit(func, bin_centers, n, p0 = [10, 0, 2])

plt.plot(bin_centers, func(bin_centers, *popt), 'r-')
plt.grid(True)
plt.show()
"""

"""
#ANIMATION

traj=p1.getalltrajectories(period)

#for i in range(num_iter-1):
#    print(i,traj[i][9])


sheet_figure = plt.figure()
ax = plt.axes(xlim=(-delta, N*delta), ylim=(-0.1, 0.1))
point_set, = ax.plot([traj[0][i] for i in range(N)],[0 for i in range(N)], 'ro')
plt.axvline(x=-delta/2)
plt.axvline(x=(N-0.5)*delta)
pt_circ = plt.Circle((4, 4), 1, color='b', fill=False)
ax.add_artist(pt_circ)
#ttl = ax.text(.5, 1.005, '', transform = ax.transAxes)

def run_animation(i):
    point_set.set_data([traj[i][k] for k in range(N)],[0 for k in range(N)])
    #ttl.set_text(str(i))    
    return point_set,

anim = animation.FuncAnimation(sheet_figure, run_animation, frames=(num_iter-1)/period, interval=5)

plt.show()
"""




os.system("rm ./lixo/*.pkl")
