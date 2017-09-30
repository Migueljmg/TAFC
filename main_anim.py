from Unmerged_anim import Problem
import random
import numpy as np
import matplotlib.pyplot as plt

#General way of creating 2D arrays
#a = [[0] * ncol for i in range(nrow)]

N=1000
wpdt=0.02
num_iter=100
delta=0.1 #1/10 of the Debye length
init_pos=[delta*i for i in range(N)]
#init_pos[0]=0
#init_pos[1]=0
#init_pos[2]=0

init_vel=[random.uniform(-0.2,0.2) for i in range(N)]
#init_vel=[0.1*pow(-1,i)*i/N for i in range(N)]


#init_vel=[0 for i in range(N)]
#init_vel[0]=8
#init_vel[1]=8
#init_vel[2]=8

p1 = Problem(N,init_pos,init_vel,wpdt,delta,num_iter)

p1.start()

