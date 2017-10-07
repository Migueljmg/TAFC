import math
from scipy import integrate
from solver_pbc_fn import Problem
from scipy.optimize import curve_fit
from scipy.integrate import simps
import random
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os
from time import time
import math


def plot_energy(period):
    b=p1.getallenergy(period)
    plt.plot(range(len(b)),b)
    plt.grid(True)
    plt.show()

def plot_vel_distr(num_iter,period):
    vlist=p1.getvelocities(num_iter)#num_iter-1-(num_iter-1)%period)
    n,bins,patches=plt.hist(vlist, bins=20)
    bin_centers=bins[:-1]+0.5*(bins[1:]-bins[:-1])

    vlist2=[vlist[i]**2 for i in range(len(vlist))]
    meansqr=math.sqrt(sum(vlist2)/len(vlist2))

    def func(x, maximum_y, x0, sigma):
        if(0<sigma<100): #bounds on sigma
            return maximum_y*np.exp(-(x-x0)**2/(2*sigma**2))
        return 1e40

    from scipy.optimize import curve_fit
    # Edit the starting values and the function to adapt for your need
    #param_bound=([-100000,-100000,0],[100000,100000,1])
    popt, pcov = curve_fit(func, bin_centers, n,p0 = [10, 0, 0.4])
    
    """
    popt2=[popt[0],0,meansqr]
    plt.plot(bin_centers, func(bin_centers, *popt), 'r-')
    plt.plot(bin_centers, func(bin_centers, *popt2), 'g-')
    plt.grid(True)
    plt.show()
    """


    return popt



def plot_vel_distr_cumul(num_iter,period,min_it):

    vlist=[]
    for i in range(min_it,num_iter):
        if(i%period==0):
            for l in range(len(p1.getvelocities(0))):
                vlist.append(p1.getvelocities(i)[l])

    #vlist=p1.getvelocities(num_iter)
    n,bins,patches=plt.hist(vlist, bins=20)
    bin_centers=bins[:-1]+0.5*(bins[1:]-bins[:-1])


    def func(x, maximum_y, x0, sigma):
        if(0<sigma<100): #bounds on sigma
            return maximum_y*np.exp(-(x-x0)**2/(2*sigma**2))
        return 1e40

    from scipy.optimize import curve_fit
    # Edit the starting values and the function to adapt for your need
    #param_bound=([-100000,-100000,0],[100000,100000,1])
    popt, pcov = curve_fit(func, bin_centers, n,p0 = [10, 0, 0.4])


    #sqrt(<v*v>)
    vlist2=[vlist[i]**2 for i in range(len(vlist))]
    meansqr=math.sqrt(sum(vlist2)/len(vlist2))

    return (popt,meansqr)






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


#INTEGRAL OF v*v*f(v)
def vintegral(k,divs,N):

    vel_vec=[]

    vel_vec.append(p1.getvelocities(k))

    maxk=max(vel_vec[0])
    mink=min(vel_vec[0])

    deltav=(maxk-mink)/divs
    v_distr=[0 for i in range(divs)]


    for i in range(divs):
        for j in range(N):
            if(vel_vec[0][j]>mink+i*deltav and vel_vec[0][j]<mink+(i+1)*deltav):
                v_distr[i]+=1


    vsquared=[] # v*v*f(v)
    for i in range(divs):
        vsquared.append((mink+(i+1./2.)*deltav)*(mink+(i+1./2.)*deltav)*v_distr[i])

    x=[mink+(i+1./2.*deltav) for i in range(divs)]
    y=vsquared

    return simps(y,x)





def timevintegral(k,divs,N,period):

    vel_vec=[]
    vecmax=[]
    vecmin=[]

    ind=0
    for j in range(k+1):
        if(j%period==0):
            vel_vec.append(p1.getvelocities(j))

            vecmax.append(max(vel_vec[ind]))
            vecmin.append(min(vel_vec[ind]))
            
            ind+=1

    #ind now corresponds to the number of time samples

    maxk=max(vecmax)
    mink=min(vecmin)


    deltav=(maxk-mink)/divs
    v_distr=[0 for i in range(divs)]


    for l in range(ind):
        for i in range(divs):
            for j in range(N):
                if(vel_vec[l][j]>mink+i*deltav and vel_vec[l][j]<mink+(i+1)*deltav):
                    v_distr[i]+=1


    vsquared=[] # v*v*f(v)
    for i in range(divs):
        vsquared.append((mink+(i+1./2.)*deltav)*(mink+(i+1./2.)*deltav)*v_distr[i]/(ind)) #Division by ind to get coherent results

    x=[mink+(i+1./2.*deltav) for i in range(divs)]
    y=vsquared

    return simps(y,x)




#General way of creating 2D arrays
#a = [[0] * ncol for i in range(nrow)]
os.system("rm ./lixo/*.pkl")


random.seed(42)
t0=time()
N=1000
wpdt=0.0001
num_iter=100000
delta=0.1 #1/10 of the Debye length
nr_images=0 #number of image charges
period=1000
init_pos=[delta*i for i in range(N)]
init_vel=[random.uniform(-0.5,0.5) for i in range(N)]


#init_vel[0]=3
#init_vel[1]=-0.7
#init_vel[2]=0.05
#init_vel[3]=0.05
p1 = Problem(N,init_pos,init_vel,wpdt,delta,num_iter,nr_images,period)

p1.start()

###MEAN SQUARE VELOCITY IN EQUILIBRIUM

#ele atinge o equilibrio por volta dos 40000

"""
velvec=[]
count=0
for i in range(num_iter-1):
    if(i%period==0 and i>40000):
        velvec.append(p1.getvelocities(i))
        count+=1

#vels=p1.getvelocities(num_iter-1- (num_iter-1)%period)

vsqmedi=[]
for i in range(count):
    vsqmedi.append(sum([vel*vel/N/count for vel in velvec[i]]))

#plt.plot(vsqmedi)
#plt.show()

vsqmed=sum(vsqmedi)

print(vsqmed)
"""

"""
####EQUILIBRIUM VALUE OF SIGMA
mu=[]
sigma=[]

for i in range(num_iter):
    if(i%period==0 and i>40000):
        param=plot_vel_distr(i,period)
        mu.append(param[1])
        sigma.append(param[2])


plt.ylim(0.26,0.6)
plt.plot(sigma)
plt.savefig("sigma5.png")
plt.show()

eq_sigma=[sigma[i] for i in range(len(sigma)-100,len(sigma))]
t_list=[i for i in range(len(sigma)-100,len(sigma))]

def func(x,eq_sigma_val):
    return eq_sigma_val

from scipy.optimize import curve_fit

popt, pcov = curve_fit(func, t_list, eq_sigma, p0 = [0.3])
print(popt[0])
"""


#plt.plot(mu)
#plt.show()

##CUMULATIVE EQUILIBRIUM VALUE OF SIGMA
min_it=40000
param=plot_vel_distr_cumul(num_iter,period,min_it)
print("sigma ", param[0][2])
print("sqrt(<v*v>) ", param[1])


"""
###INTEGRAL OF V*V*F(V)
divs=20 #number of intervals in the distribution function
integral_list=[] #list with the integrals of v*v*f(v)
for k in range(num_iter-1):
    if(k%period==0):
        integral_list.append(vintegral(k,divs,N)) #(k,divs,N)

plt.plot(integral_list)
plt.savefig("test3.png")
plt.show()
"""


#p1.start()
t1=time()

#print(t1-t0)
#p1.cenas()


#plot_energy(period)



#plot_vel_distr(num_iter-1- (num_iter-1)%period,period)



"""
###INTEGRAL OF V*V*F(V) WITH SUM IN TIME
divs=20 #number of intervals in the distribution function
integral_list=[] #list with the integrals of v*v*f(v)
for k in range(num_iter-1):
    if(k%period==0):
        integral_list.append(timevintegral(k,divs,N,period)) #(k,divs,N)

plt.plot(integral_list)
plt.savefig("test4.png")
plt.show()

plot_vel_distr(num_iter-1,period)
"""




#animate(num_iter,period)


##Physical process 

"""
#ANIMATION OF THE DEVIATIONS
random.seed(42)
t0=time()
N=2000
wpdt=0.005
num_iter=3000
delta=0.1 #1/10 of the Debye length
nr_images=0 #number of image charges
period=50
init_pos=[delta*i for i in range(N)]
init_vel=[]

for i in range(N/4):
    init_vel.append(0)

for i in range(N/4,3*N/4):
    init_vel.append(0)#math.sin(2*math.pi*(i-N/4)/N))

for i in range(3*N/4,N):
    init_vel.append(0)

fastpar=50
v_fastpar=2
for i in range(fastpar):
    init_vel[i]=v_fastpar
    init_pos[i]=i*delta/10


#plt.plot(init_vel)
#plt.show()

p1 = Problem(N,init_pos,init_vel,wpdt,delta,num_iter,nr_images,period)

p1.start()

eq_pos=p1.getalleqpos(period)
traj=p1.getalltrajectories(period)

deviation = [[0] * N for i in range(len(eq_pos))]


for k in range((num_iter-1)/period):
    for i in range(N):
        deviation[k][i]=traj[k][i]-eq_pos[k][i]

fig = plt.figure()

# animation function.  This is called sequentially
def animate(i):
    plt.cla()
    #plt.hist(self.pos[i],normed=False,bins=200,range=[0,(self.N-1)*self.delta])
    plt.ylim(-2.5,2.5)
    plt.plot(deviation[i])


# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, frames=len(eq_pos)-1, interval=100)

plt.show()
"""






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
