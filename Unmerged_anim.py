import math
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
import time

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
        self.energy=[0 for n in range(n_iter)]

    def start(self):
        "num_iter -> number of iterations; deltat-> iteration time"
        self.pos[0]=self.init_pos
        self.vel[0]=self.init_vel
        
        argument=self.a
        a=self.a
        cos=math.cos(argument)
        sin=math.sin(argument)
        pos=self.pos
            
        for k in range(1,self.num_iter):
            for i in range(self.N):
                self.vel[k][i]=self.vel[k-1][i]*cos-sin*(self.pos[k-1][i]-self.eq_pos[i])
                self.pos[k][i]=self.pos[k-1][i]+self.vel[k-1][i]*sin-(self.pos[k-1][i]-self.eq_pos[i])*(1-cos)
                

                
                #Verifying if there were any crossings taking place
                if(i>0):
                    if(self.pos[k][i-1]>self.pos[k][i]):

                        #Computation of tc1
                        tc1 = self.a*(self.pos[k-1][i]-self.pos[k-1][i-1])/(self.pos[k-1][i]-self.pos[k-1][i-1] + self.pos[k][i]-self.pos[k][i-1])
                        # Position at t + tc1
                        pos_tc1_im1 = self.pos[k-1][i-1]+self.vel[k-1][i-1]*math.sin(tc1)-(self.pos[k-1][i-1]-self.eq_pos[i-1])*(1-math.cos(tc1))## i-1 position
                        pos_tc1_i = self.pos[k-1][i]+self.vel[k-1][i]*math.sin(tc1)-(self.pos[k-1][i]-self.eq_pos[i])*(1-math.cos(tc1))## i position

                        #Computation of tc2
                        tc2 = self.a*(self.pos[k-1][i]-self.pos[k-1][i-1])/(self.pos[k-1][i]-self.pos[k-1][i-1] + pos_tc1_i-pos_tc1_im1)

                        # Positions at t + tc2
                        pos_tc2_im1 = self.pos[k-1][i-1]+self.vel[k-1][i-1]*math.sin(tc2)-(self.pos[k-1][i-1]-self.eq_pos[i-1])*(1-math.cos(tc2))## i-1 position
                        pos_tc2_i = self.pos[k-1][i]+self.vel[k-1][i]*math.sin(tc2)-(self.pos[k-1][i]-self.eq_pos[i])*(1-math.cos(tc2))## i position
                        #velocities at t+tc2
                        vel_tc2_im1=self.vel[k-1][i-1]*math.cos(tc2)-math.sin(tc2)*(self.pos[k-1][i-1]-self.eq_pos[i-1])## i-1 velocity
                        vel_tc2_i=self.vel[k-1][i]*math.cos(tc2)-math.sin(tc2)*(self.pos[k-1][i]-self.eq_pos[i])## i velocity

                        # Evolution from tc2 to t+dt
                        #Positions
                        self.pos[k][i-1]=pos_tc2_im1+vel_tc2_im1*math.sin(self.a-tc2)-(pos_tc2_im1-self.eq_pos[i-1])*(1-math.cos(self.a-tc2))+(self.a-tc2)*(self.a-tc2)/2*self.delta
                        self.pos[k][i]=pos_tc2_i+vel_tc2_i*math.sin(self.a-tc2)-(pos_tc2_i-self.eq_pos[i])*(1-math.cos(self.a-tc2))-(self.a-tc2)*(self.a-tc2)/2*self.delta
                        #Velocities
                        self.vel[k][i-1]=vel_tc2_im1*math.cos(self.a-tc2)-math.sin(self.a-tc2)*(pos_tc2_im1-self.eq_pos[i-1])+(self.a-tc2)*self.delta
                        self.vel[k][i]=vel_tc2_i*math.cos(self.a-tc2)-math.sin(self.a-tc2)*(pos_tc2_i-self.eq_pos[i])-(self.a-tc2)*self.delta


            ########### Applying PBC
            if(self.pos[k][self.N-1]>self.pos[k][0]+self.N*self.delta):

                #Translation in the x0 position and equilibrium position to compare with x(N-1). Now it is as x0 is in front of x(N-1)
                self.pos[k][0]=self.pos[k][0]+self.N*self.delta
                self.eq_pos0=self.eq_pos[0]+self.N*self.delta

                N=self.N

                #Computation of tc1
                tc1 = self.a*(self.pos[k-1][0]-self.pos[k-1][N-1])/(self.pos[k-1][0]-self.pos[k-1][N-1] + self.pos[k][0]-self.pos[k][N-1])
                # Position at t + tc1
                pos_tc1_im1 = self.pos[k-1][N-1]+self.vel[k-1][N-1]*math.sin(tc1)-(self.pos[k-1][N-1]-self.eq_pos[N-1])*(1-math.cos(tc1))## i-1 position which now corresponds to N-1
                pos_tc1_i = self.pos[k-1][0]+self.vel[k-1][0]*math.sin(tc1)-(self.pos[k-1][0]-self.eq_pos0)*(1-math.cos(tc1))## i position which now corresponds to x0

                #Computation of tc2
                tc2 = self.a*(self.pos[k-1][0]-self.pos[k-1][N-1])/(self.pos[k-1][0]-self.pos[k-1][N-1] + pos_tc1_i-pos_tc1_im1)

                # Positions at t + tc2
                pos_tc2_im1 = self.pos[k-1][N-1]+self.vel[k-1][N-1]*math.sin(tc2)-(self.pos[k-1][N-1]-self.eq_pos[N-1])*(1-math.cos(tc2))## i-1 position
                pos_tc2_i = self.pos[k-1][0]+self.vel[k-1][0]*math.sin(tc2)-(self.pos[k-1][0]-self.eq_pos0)*(1-math.cos(tc2))## i position
                #velocities at t+tc2
                vel_tc2_im1=self.vel[k-1][N-1]*math.cos(tc2)-math.sin(tc2)*(self.pos[k-1][N-1]-self.eq_pos[N-1])## i-1 velocity
                vel_tc2_i=self.vel[k-1][0]*math.cos(tc2)-math.sin(tc2)*(self.pos[k-1][0]-self.eq_pos0)## i velocity

                # Evolution from tc2 to t+dt
                #Positions
                self.pos[k][N-1]=pos_tc2_im1+vel_tc2_im1*math.sin(self.a-tc2)-(pos_tc2_im1-self.eq_pos[N-1])*(1-math.cos(self.a-tc2))+(self.a-tc2)*(self.a-tc2)/2*self.delta
                self.pos[k][0]=pos_tc2_i+vel_tc2_i*math.sin(self.a-tc2)-(pos_tc2_i-self.eq_pos0)*(1-math.cos(self.a-tc2))-(self.a-tc2)*(self.a-tc2)/2*self.delta
                #Velocities
                self.vel[k][N-1]=vel_tc2_im1*math.cos(self.a-tc2)-math.sin(self.a-tc2)*(pos_tc2_im1-self.eq_pos[N-1])+(self.a-tc2)*self.delta
                self.vel[k][0]=vel_tc2_i*math.cos(self.a-tc2)-math.sin(self.a-tc2)*(pos_tc2_i-self.eq_pos0)-(self.a-tc2)*self.delta

                #And finally, I need to translade N-1 by -N*delta. x0 is already in the right place. When I order after dt N-1 will be the next 0 and 0 will be the next N-1
                self.pos[k][N-1]=self.pos[k][N-1]-self.N*self.delta



            #Order all the particles according to their positions after each time interval dt
            self.vel[k]=[x for _,x in sorted(zip(self.pos[k],self.vel[k]))]
            self.pos[k].sort()



            self.energy[k]=self.energy_from_velocity(self.vel[k])

        #print(self.energy)
        t_list = [self.a*i for i in range(self.num_iter)]
        vel_list=[self.vel[k][0] for k in range(self.num_iter)]
        pos_list=[self.pos[k][0] for k in range(self.num_iter)]
        pos_list2=[self.pos[k][1] for k in range(self.num_iter)]
        pos_list3=[self.pos[k][2] for k in range(self.num_iter)]
        #pos_list4=[self.pos[k][3] for k in range(self.num_iter)]
        #pos_list5=[self.pos[k][4] for k in range(self.num_iter)]
        #pos_list6=[self.pos[k][5] for k in range(self.num_iter)]
        #pos_list7=[self.pos[k][6] for k in range(self.num_iter)]
        #pos_list8=[self.pos[k][7] for k in range(self.num_iter)]
        #pos_list9=[self.pos[k][8] for k in range(self.num_iter)]
        #plt.plot(t_list,pos_list)
        #plt.plot(t_list,pos_list,t_list,pos_list2,t_list,pos_list3)
        #plt.plot(t_list,pos_list,t_list,pos_list2,t_list,pos_list3,t_list,pos_list4,t_list,pos_list5,t_list,pos_list6,t_list,pos_list7,t_list,pos_list8,t_list,pos_list9)

        #plt.hist(self.vel[self.num_iter-1],normed=False,bins=50,range=[-1,1])
        #plt.show()


        """
        ####Animation of the sheets################
        fig = plt.figure()
        ax = plt.axes(xlim=(-0.3, (self.N+2)*self.delta), ylim=(-0.1, 0.1))
        d, = ax.plot([self.pos[0][i] for i in range(self.N)],
                     [0 for i in range(self.N)], 'ro')
        plt.axvline(x=-self.delta/2)
        plt.axvline(x=(self.N-0.5)*self.delta)
        circle = plt.Circle((5, 5), 1, color='b', fill=False)
        ax.add_artist(circle)


        # animation function.  This is called sequentially
        def animate(i):
            d.set_data([self.pos[i][k] for k in range(self.N)],
                       [0 for k in range(self.N)])
            return d,

                # call the animator.  blit=True means only re-draw the parts that have changed.
        anim = animation.FuncAnimation(fig, animate, frames=self.num_iter, interval=30)

        #plt.show()
        """

        ###Animation of the position histogram###
        fig = plt.figure()

        # animation function.  This is called sequentially
        def animate(i):
            plt.cla()
            #plt.hist(self.pos[i],normed=False,bins=200,range=[0,(self.N-1)*self.delta])
            plt.hist(self.vel[i],normed=True,bins=50,range=[-0.1,0.1])

        # call the animator.  blit=True means only re-draw the parts that have changed.
        anim = animation.FuncAnimation(fig, animate, frames=self.num_iter, interval=50)

        plt.show()

    def energy_from_velocity(self,vel):
        return sum([vels*vels for vels in vel])
    def momentum(self,vel):
        return sum (vel)

                        
                        
                        
                        
            

            
                
                    
