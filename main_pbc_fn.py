from solver_pbc_fn import Problem
import random
import numpy as np
import matplotlib.pyplot as plt

#General way of creating 2D arrays
#a = [[0] * ncol for i in range(nrow)]

random.seed(42)

N=10#4
wpdt=0.001
num_iter=5000
delta=0.1 #1/10 of the Debye length
nr_images=0 #I'm sending this but it'll not matter, I have to remove it from the class
init_pos=[delta*i for i in range(N)]
init_vel=[1 for i in range(N)]
#init_vel[0]=-3
#init_vel[1]=-0.5


p1 = Problem(N,init_pos,init_vel,wpdt,delta,num_iter,nr_images)

p1.start()




"""
random.seed(42)

N=50
wpdt=0.05
num_iter=20
delta=0.1 #1/10 of the Debye length
nr_images=0 #number of image charges
init_pos=[delta*i for i in range(N)]
init_vel=[random.uniform(-1,1) for i in range(N)]
init_vel[0]=0.5
init_vel[1]=-0.5
"""


"""
##########Examples of image sheets list handling
#Velocities
n=3
for i in range(n):
    init_vel.append(init_vel[i])

for i in range(n):
    init_vel.insert(0,init_vel[N-1])#don't need to put i as I am inserting elements at the left of the list

print(init_vel)

#Positions
n=5
for i in range(n):
    init_pos.append(init_pos[i]+N*delta)

for i in range(n):
    init_pos.insert(0,init_pos[N-1]-N*delta)#don't need to put i as I am inserting elements at the left of the list

print(init_pos)
"""
