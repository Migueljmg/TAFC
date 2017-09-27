import math

class Problem:
    def __init__(self, number,initial_positions, initial_velocities,deltat,plasma_freq,delta1):
        self.init_pos=sorted(initial_positions)
        self.init_vel=[x for _,x in sorted(zip(initial_velocities,initial_positions]
        self.N=number #number of particles
        self.pos=[[]] #final matrix of positions 
        self.vel=[[]] #final matrix of velocities
        self.omegap=plasma_freq
        self.delta=delta1
        self.eq_pos=[self.delta*n for n in range(N)]

    def start(self, num_iter,deltat):
        "num_iter -> number of iterations; deltat-> iteration time"
        self.pos[0]=self.init_pos
        self.vel[0]=self.init_vel
        
        argument=self.omegap*deltat
        cos=math.cos(argument)
        sin=math.sin(argument)

        for k in range(1,num_iter):
            for i in range(N):
                self.vel[k][i]=self.vel[k-1,i]*cos-self.omegap*sin*(self.pos[k-1][i]-self.eq_pos[i])
                self.pos[k][i]=self.pos[k-1][i]+self.vel[k-1,i]*sin-(self.pos[k-1][i]-self.eq_pos[i])*(1-cos)
