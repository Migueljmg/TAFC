from Unmerged_anim import Problem
import random
import numpy as np
import matplotlib.pyplot as plt

#General way of creating 2D arrays
#a = [[0] * ncol for i in range(nrow)]

N=3
num_iter=1000
delta=0.1 #1/10 of the Debye length
init_pos=[delta*i for i in range(N)]
#init_pos[0]=0
#init_pos[1]=0.01
#init_pos[2]=0.02
#init_vel=[random.uniform(-1,1) for i in range(N)]#velocities such that there are no crossings to test
init_vel=[0 for i in range(N)]
init_vel[0]=-0.06

#init_vel[0]=0.06
#init_vel[1]=-0.06
#init_vel[2]=-0.06


#init_vel=[0.055,-0.055]
wpdt=0.04

p1 = Problem(N,init_pos,init_vel,wpdt,delta,num_iter)

p1.start()

