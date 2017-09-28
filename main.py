from Unmerged_mine import Problem
import random
import numpy as np
import matplotlib.pyplot as plt

#General way of creating 2D arrays
#a = [[0] * ncol for i in range(nrow)]

N=10
num_iter=200
delta=0.1 #1/10 of the Debye length
init_pos=[delta*i for i in range(N)]
init_vel=[random.uniform(-0.01,0.01) for i in range(N)]#velocities such that there are no crossings to test
#init_vel=[0.01 for i in range(N)]
wpdt=0.05

print(init_vel) #para verificar que da igual ao que entra na classe (initial_vels). Aquilo do zip como estava antes fazia com que o self.init_vel recebesse posicoes em vez de vels. Troquei e parece tar bem agora. Podes apagar isto dps.

p1 = Problem(N,init_pos,init_vel,wpdt,delta,num_iter)

p1.start()

