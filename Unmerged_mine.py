import math
import matplotlib.pyplot as plt

class Problem:
    def __init__(self, number,initial_positions, initial_velocities,wpdt,delta1,n_iter):
        self.init_pos=sorted(initial_positions) #in units of Debye length
        #self.init_vel=[x for _,x in sorted(zip(initial_velocities,initial_positions))] #in units of thermal velocity #acho que isto estava ao contrario
        self.init_vel=[x for _,x in sorted(zip(initial_positions,initial_velocities))] #in units of thermal velocity
        self.N=number #number of particles
        self.num_iter=n_iter #number of iterations
        self.pos=[[0] * self.N for i in range(self.num_iter)] #final matrix of positions 
        self.vel=[[0] * self.N for i in range(self.num_iter)] #final matrix of velocities
        self.a=wpdt #adimensional measure of time
        self.delta=delta1 #in units of the Debye length
        self.eq_pos=[self.delta*n for n in range(self.N)]

    def start(self):
        "num_iter -> number of iterations; deltat-> iteration time"
        self.pos[0]=self.init_pos
        self.vel[0]=self.init_vel
        
        argument=self.a
        cos=math.cos(argument)
        sin=math.sin(argument)

        print(self.init_vel) #VER main. Isto e so para verificar, podes apagar dps

        for k in range(1,self.num_iter):
            for i in range(self.N):
                self.vel[k][i]=self.vel[k-1][i]*cos-sin*(self.pos[k-1][i]-self.eq_pos[i])
                self.pos[k][i]=self.pos[k-1][i]+self.vel[k-1][i]*sin-(self.pos[k-1][i]-self.eq_pos[i])*(1-cos)

        t_list = [self.a*i for i in range(self.num_iter)]
        vel_list=[self.vel[k][3] for k in range(self.num_iter)]
        pos_list=[self.pos[k][3] for k in range(self.num_iter)]
        plt.plot(t_list,pos_list)#I'm printing particle's 3 positions as an example
        plt.show()
            
                
                    
