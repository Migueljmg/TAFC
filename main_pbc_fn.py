#File responsible for collecting initial conditions
#for the problem (number of particles, delta, etc.)
#and for producing the corresponding results

import math #sqrt
from scipy import integrate #perform integrations 
                            #(second momentum of distrib.) 
from solver_pbc_fn import Problem #class with the methods 
                                  #to solve the model
from scipy.optimize import curve_fit #fits
from scipy.integrate import simps #integration with Simpson's method
import random #eases the creation of initial conditions
import numpy as np 
import matplotlib.pyplot as plt #plots
import matplotlib.animation as animation #animations
import os #create directory to store files, and delete such files 
          #when they are no longer needed
from time import time #check runtime
import sys

#Function to plot the total energy of the system
def plot_energy(wpdt,num_iter,period):

    t_list=[]
    for i in range(num_iter-1):
        if(i%period==0):
            t_list.append(i*wpdt)

    b=p1.getallenergy(period)

    energy_dev=(b[len(b)-1]-b[0])/b[0]

    plt.plot(t_list,b)
    plt.xlabel(r'$w_p t$')
    plt.ylabel(r'Energy')
    plt.grid(True)
    plt.savefig("energy.png")
    plt.show()


#Function to return the energy deviation from the initial energy after a time num_iter*wpdt
def energy_deviation(wpdt,num_iter,period):

    b=p1.getallenergy(period)
    energy_dev=(b[len(b)-1]-b[0])/b[0]

    return energy_dev



def plot_vel_distr(num_iter,period,nbins):
    vlist=p1.getvelocities(num_iter-1-(num_iter-1)%period)
    n,bins,patches=plt.hist(vlist, bins=nbins)
    bin_centers=bins[:-1]+0.5*(bins[1:]-bins[:-1])

    plt.xlabel("velocity")
    plt.ylabel("counts")
    plt.savefig("vel_histogram.png")
    print("\n vel_histrogram.png created \n")

    plt.show()





#Function to arrange the velocities of all sheets calculated at some iteration num_iter in a histogram and to fit it to a gaussian distribution. It returns the fit parameters
def fit_vel_distr(num_iter,period):
    vlist=p1.getvelocities(num_iter)
    n,bins=np.histogram(vlist, bins=20)
    bin_centers=bins[:-1]+0.5*(bins[1:]-bins[:-1])

    def func(x, maximum_y, x0, sigma):
        if(0<sigma<100): #bounds on sigma
            return maximum_y*np.exp(-(x-x0)**2/(2*sigma**2))
        return 1e40

    from scipy.optimize import curve_fit
    # Edit the starting values and the function to adapt for your need
    # p0 is the vector that has those values (the second is the average,
    # the third is the sigma and the first is the height at x=x0
    popt, pcov = curve_fit(func, bin_centers, n,p0 = [10, 0, 0.4])

    return popt


#This function is responsible for two computations:
#- Arrange the summed velocities from iteration min_it to num_iter into a histogram to be fitted. It plots the corresponding histogram along with the fit curve and returns the fit parameters.
#- Compute the mean squared velocity based on the summed velocities mentioned above with the purpose of comparing it to the fit parameter sigma
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

    popt, pcov = curve_fit(func, bin_centers, n,p0 = [10, 0, 0.4])


    #sqrt(<v*v>)
    vlist2=[vlist[i]**2 for i in range(len(vlist))]
    meansqr=math.sqrt(sum(vlist2)/len(vlist2))

    error=np.sqrt(n)

    plt.plot(bin_centers, func(bin_centers, *popt), 'r-')
    plt.errorbar(bin_centers,n,yerr=error,fmt='none')
    plt.xlabel("velocity")
    plt.ylabel("counts")
    plt.savefig("cumulated_histogram.png")
    plt.show()

    return (popt,meansqr)



#Function to provide an animation of the sheets, represented by red circles, untill iteration num_iter
def animate(wpdt,num_iter,period):
    traj=p1.getalltrajectories(period)

    sheet_figure = plt.figure()
    ax = plt.axes(xlim=(-delta, N*delta), ylim=(-0.1, 0.1))
    point_set, = ax.plot([traj[0][i] for i in range(N)],[0 for i in range(N)], 'ro')
    plt.axvline(x=-delta/2)
    plt.axvline(x=(N-0.5)*delta)
    plt.xlabel("x")
    pt_circ = plt.Circle((4, 4), 1, color='b', fill=False)
    ax.add_artist(pt_circ)
    ttl = ax.text(.45, 1.005, '', transform = ax.transAxes)

    def run_animation(i):
        point_set.set_data([traj[i][k] for k in range(N)],[0 for k in range(N)])
        ttl.set_text(r"$w_p t$ = " + str(i*wpdt*period))
        return point_set,

    anim = animation.FuncAnimation(sheet_figure, run_animation, frames=(num_iter-1)/period, interval=10)

    plt.show()



os.system("rm ./lixo/*.pkl")


#CHECK THE ARGUMENTS SENT BY THE USER
words=[str(sys.argv[i]) for i in range(len(sys.argv)) if sys.argv[i][0]=='-']
word=""
for case in words:
    word+=case

if(word=="" or "h" in word):
    print("\n HELP \n \n")
    print("-a for animation \n")
    print("-e for energy diagnosis \n")
    print("-t for trajectories' plot \n")
    print("-v for velocity distribution histogram \n")
    print("-p for thermalization physical process analysis \n")



#ANIMATION
if("a" in word):
    #Parameters of the problem
    random.seed(42)
    t0=time() #starting time
    N=10 #Number of sheets
    wpdt=0.001 #time step in units of wp
    num_iter=10000 #number of iterations
    delta=0.1 #distance between sheets in units of the Debye length
    period=10 #iteration interval between sampling points
    init_pos=[delta*i for i in range(N)] #array with the initial positions of the sheets
    init_vel=[random.uniform(-0.2,0.2) for i in range(N)] #array with the initial velocities of the sheets

    p1 = Problem(N,init_pos,init_vel,wpdt,delta,num_iter,period) #initialization of an object of the class Problem responsible for the computation of the sheets' motion

    p1.start() #start the computations

    animate(wpdt,num_iter,period)



#ENERGY CONSERVATION CHECK
if("e" in word):
    #Parameters of the problem
    random.seed(42)
    t0=time() #starting time
    N=20 #Number of sheets
    wpdt=0.001 #time step in units of wp
    num_iter=10000 #number of iterations
    delta=0.1 #distance between sheets in units of the Debye length
    period=10 #iteration interval between sampling points
    init_pos=[delta*i for i in range(N)] #array with the initial positions of the sheets
    init_vel=[random.uniform(-0.5,0.5) for i in range(N)] #array with the initial velocities of the sheets

    p1 = Problem(N,init_pos,init_vel,wpdt,delta,num_iter,period) #initialization of an object of the class Problem responsible for the computation of the sheets' motion

    p1.start() #start the computations

    plot_energy(wpdt,num_iter,period)
    energy_dev=energy_deviation(wpdt,num_iter,period)

    print(" Energy deviation: " +  str(energy_dev))

    print("\n energy.png created \n")



#TRAJECTORIES
if("t" in word):
    #Parameters of the problem
    random.seed(42)
    t0=time() #starting time
    N=10 #Number of sheets
    wpdt=0.001 #time step in units of wp
    num_iter=10000 #number of iterations
    delta=0.1 #distance between sheets in units of the Debye length
    period=1 #iteration interval between sampling points
    init_pos=[delta*i for i in range(N)] #array with the initial positions of the sheets
    init_vel=[random.uniform(-0.08,0.08) for i in range(N)] #array with the initial velocities of the sheets

    p1 = Problem(N,init_pos,init_vel,wpdt,delta,num_iter,period) #initialization of an object of the class Problem responsible for the computation of the sheets' motion

    p1.start() #start the computations

    traj=p1.getalltrajectories(period)

    for i in range(N):
        plt.plot(range(len(traj)), [c[i] for c in traj])

    plt.show()


#HISTOGRAM OF VELOCITIES
if("v" in word):
    #Parameters of the problem
    random.seed(42)
    t0=time() #starting time
    N=1000 #Number of sheets
    wpdt=0.0005 #time step in units of wp
    num_iter=10000 #number of iterations
    delta=0.1 #distance between sheets in units of the Debye length
    period=20 #iteration interval between sampling points
    init_pos=[delta*i for i in range(N)] #array with the initial positions of the sheets
    init_vel=[random.uniform(-0.3,0.3) for i in range(N)] #array with the initial velocities of the sheets

    p1 = Problem(N,init_pos,init_vel,wpdt,delta,num_iter,period) #initialization of an object of the class Problem responsible for the computation of the sheets' motion

    p1.start() #start the computations

    nbins=20
    plot_vel_distr(num_iter-1,period,nbins) #plot histogram of velocities for iteration num_iter-1




#PHYSICAL PROCESS - THERMALIZATION
if("p" in word):
    #Parameters of the problem. Notice that one must should also choose the parameters min_it and starting_time at the beggining of points 2 and 3, respectively (these are explained there)
    random.seed(42)
    N=500 #Number of sheets
    wpdt=0.0005 #time step in units of wp
    num_iter=14000 #number of iterations
    delta=0.1 #distance between sheets in units of the Debye length
    period=250 #iteration interval between sampling points
    init_pos=[delta*i for i in range(N)] #array with the initial positions of the sheets
    init_vel=[random.uniform(-0.2,0.2) for i in range(N)] #array with the initial velocities of the sheets

    p1 = Problem(N,init_pos,init_vel,wpdt,delta,num_iter,period) #initialization of an object of the class Problem responsible for the computation of the sheets' motion

    p1.start() #start the computations

    
    #1 - Check energy conservation
    plot_energy(wpdt,num_iter,period)



    #2 - Comparison between equilibrium values of sigma and sqrt(<v*v>)
    #For these computations the velocity distributions are summed in a cumulative way from min_it to num_iter, "period" after "period". min_it is chosen such that it already corresponds to an equilibrium situation 
    min_it=6000

    if(min_it>=num_iter):
        print("Warning! min_it must be smaller than num_iter. Try again.")
        sys.exit(0)

    param=plot_vel_distr_cumul(num_iter,period,min_it) 
    print("sigma ", param[0][2])
    print("sqrt(<v*v>) ", param[1])

    

    #3 - Plot of the variation of sigma squared with time from starting_time to num_iter, "period" after "period"
    starting_time=1000

    if(starting_time>=num_iter):
        print("Warning! starting_time must be smaller than num_iter. Try again.")
        sys.exit(0)

    sigma=[]
    iter_list=[]

    for i in range(starting_time,num_iter):
        if(i%period==0):
            iter_list.append(i)
            param=fit_vel_distr(i,period)
            sigma.append(param[2]*param[2])


    t_list=[iter_list[i]*wpdt for i in range(len(sigma))]
    plt.xlabel(r'$w_p t$')
    plt.ylabel(r'$\sigma^2$')
    plt.plot(t_list,sigma,'o-')
    plt.savefig("sigma_time.png")

    plt.show()

    print("\n energy.png, cumulated_histogram.png and sigma_time.png created \n")
