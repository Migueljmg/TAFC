import math
import matplotlib.pyplot as plt
import matplotlib.animation as animation

class Problem:
    def __init__(self, number,initial_positions, initial_velocities,wpdt,delta1,n_iter,n_image):

        self.a=wpdt #adimensional measure of time
        self.delta=delta1 #in units of the Debye length
        self.n_i=n_image #number of image charges at each side of the system
        self.init_pos=sorted(initial_positions) #in units of Debye length
        self.init_vel=[x for _,x in sorted(zip(initial_positions,initial_velocities))] #in units of thermal velocity
        self.N=number #number of particles
        self.num_iter=n_iter #number of iterations

        ##Image charges initial positions and velocities
        #Velocity
        for i in range(self.n_i):
            self.init_vel.append(self.init_vel[i])

        for i in range(self.n_i):
            self.init_vel.insert(0,self.init_vel[self.N-1])#don't need to put i as I am inserting elements at the left of the list

        #Positions
        for i in range(self.n_i):
            self.init_pos.append(self.init_pos[i]+self.N*self.delta)

        for i in range(self.n_i):
            self.init_pos.insert(0,self.init_pos[self.N-1]-self.N*self.delta)


        self.pos=[[0] * (self.N+2*self.n_i) for i in range(self.num_iter)] #final matrix of positions 
        self.vel=[[0] * (self.N+2*self.n_i) for i in range(self.num_iter)] #final matrix of velocities
        self.eq_pos=[self.delta*n for n in range(self.N)]
        
        #Equilibrium positions for the image charges
        for i in range(self.n_i):
            self.eq_pos.append(self.eq_pos[i]+self.N*self.delta)

        for i in range(self.n_i):
            self.eq_pos.insert(0,self.eq_pos[self.N-1]-self.N*self.delta)


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

        col_counter=0
            
        for k in range(1,self.num_iter):


            #The evolution from the equations only matters from n to n+N-1. The rest are image sheets and are completely determined 
            for i in range(self.n_i,self.n_i+self.N):
                    
                self.vel[k][i]=self.vel[k-1][i]*cos-sin*(self.pos[k-1][i]-self.eq_pos[i])
                self.pos[k][i]=self.pos[k-1][i]+self.vel[k-1][i]*sin-(self.pos[k-1][i]-self.eq_pos[i])*(1-cos)

            #print(self.pos[k])

            for i in range(self.n_i):
                #Right images
                self.vel[k][i+self.n_i+self.N]=self.vel[k][i+self.n_i]
                self.pos[k][i+self.n_i+self.N]=self.pos[k][i+self.n_i]+self.N*self.delta

            for i in range(self.n_i-1,-1,-1):
                #Left images
                self.vel[k][i]=self.vel[k][i+self.N]
                self.pos[k][i]=self.pos[k][i+self.N]-self.N*self.delta
            #print(self.pos[k])


            #Verifying if there were any crossings taking place
            for i in range(1,self.N+2*self.n_i):

                if(self.pos[k][i-1]>self.pos[k][i]):

                    if(self.n_i-1<i-1<self.n_i+self.N):
                        col_counter=col_counter+1

                   #Computation of tc1
                    tc1 = self.a*(self.pos[k-1][i]-self.pos[k-1][i-1])/(self.pos[k-1][i]-self.pos[k-1][i-1] + self.pos[k][i-1]-self.pos[k][i])
                    # Position at t + tc1
                    pos_tc1_im1 = self.pos[k-1][i-1]+self.vel[k-1][i-1]*math.sin(tc1)-(self.pos[k-1][i-1]-self.eq_pos[i-1])*(1-math.cos(tc1))## i-1 position
                    pos_tc1_i = self.pos[k-1][i]+self.vel[k-1][i]*math.sin(tc1)-(self.pos[k-1][i]-self.eq_pos[i])*(1-math.cos(tc1))## i position


                    #Computation of tc2
                    tc2 = tc1*(self.pos[k-1][i]-self.pos[k-1][i-1])/(self.pos[k-1][i]-self.pos[k-1][i-1] + pos_tc1_im1-pos_tc1_i)
                    #tc2= -math.atan((self.pos[k-1][i]-self.pos[k-1][i-1])/(self.vel[k-1][i]-self.vel[k-1][i-1]))


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



            #Order all the particles according to their positions after each time interval dt
            self.vel[k]=[x for _,x in sorted(zip(self.pos[k],self.vel[k]))]
            self.pos[k].sort()


            self.energy[k]=self.energy_from_velocity(self.particle_vel(k))+self.energy_from_position(self.particle_pos(k), self.particle_eq_pos())

            
        #print(col_counter/self.num_iter)
 


        #energy conservation?
        #print(self.particle_pos(4),self.eq_pos)
        """
        lm=0
        rm=self.num_iter
        t_list=[i for i in range(lm ,rm)]
        pos_list1=[self.pos[i][0] for i in range(lm ,rm)]
        pos_list2=[self.pos[i][1] for i in range(lm ,rm)]
        pos_list3=[self.pos[i][2] for i in range(lm ,rm)]
        pos_list4=[self.pos[i][3] for i in range(lm ,rm)]
        pos_list5=[self.pos[i][4] for i in range(lm ,rm)]
        pos_list6=[self.pos[i][5] for i in range(lm ,rm)]
        pos_list7=[self.pos[i][6] for i in range(lm ,rm)]
        pos_list8=[self.pos[i][7] for i in range(lm ,rm)]
        """
        lm=0
        rm=self.num_iter
        t_list=[i for i in range(lm ,rm)]
        #pos_list=[[self.pos[i][l] for i in range(lm,rm)] for l in range(self.N)]
        #for i in range(self.N):
        #    plt.plot(t_list,pos_list[i])


        #plt.axhline(y=self.delta*(self.N-0.5),xmin=0,xmax=self.num_iter,c="black")
        #plt.axhline(y=-self.delta/2,xmin=0,xmax=self.num_iter,c="black")
        #plt.plot(t_list,pos_list1,t_list,pos_list2,t_list,pos_list3,t_list,pos_list4,t_list,pos_list5,t_list,pos_list6,t_list,pos_list7,t_list,pos_list8,t_list,self.energy)
        


        self.energy[0]=self.energy[1]
        #plt.plot(self.energy)

        plt.xlabel('iters)')
        plt.ylabel('energy')
        plt.title('About as simple as it gets, folks')
        plt.grid(True)
        #plt.savefig("test.png")
        #plt.show()

        #Histogram of the velocity distribution
        vlist=[self.vel[self.num_iter-1][i] for i in range(self.n_i,self.n_i+self.N)]
        plt.hist(vlist,normed=False,bins=20)#,range=[-0.95,0.95])

        plt.show()
        



        ######Animation with sheets####################
        sheet_figure = plt.figure()
        ax = plt.axes(xlim=(-20*self.delta-self.n_i*self.delta-self.delta, (self.N+self.n_i+1+20)*self.delta), ylim=(-0.1, 0.1))
        point_set, = ax.plot([self.pos[0][i] for i in range(self.N+2*self.n_i)],[0 for i in range(self.N+2*self.n_i)], 'ro')
        plt.axvline(x=-self.delta/2)
        plt.axvline(x=(self.N-0.5)*self.delta)
        pt_circ = plt.Circle((4, 4), 1, color='b', fill=False)
        ax.add_artist(pt_circ)
        ttl = ax.text(.5, 1.005, '', transform = ax.transAxes)

        def run_animation(i):
            point_set.set_data([self.pos[i][k] for k in range(self.N+2*self.n_i)],[0 for k in range(self.N+2*self.n_i)])
            ttl.set_text(str(i))    
            return point_set,

        anim = animation.FuncAnimation(sheet_figure, run_animation, frames=self.num_iter, interval=50)


        #plt.show()


    def energy_from_velocity(self,vel):
        return sum([vels*vels for vels in vel])

    def energy_from_position(self,pos,pos_eq):
        diff=[pos[i]-pos_eq[i] for i in range(len(pos))]
        return sum([diffs*diffs for diffs in diff])

    def momentum(self,vel):
        return sum (vel)


    def particle_pos(self,k):
        return [self.pos[k][i] for i in range(self.n_i,self.n_i+self.N)]

    def particle_vel(self,k):
        return [self.vel[k][i] for i in range(self.n_i,self.n_i+self.N)]

    
    def particle_eq_pos(self):
        return [self.eq_pos[i] for i in range(self.n_i, self.n_i+self.N)]
 

                        
                        
                        
                        
            

            
                
                    
