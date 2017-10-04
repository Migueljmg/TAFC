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


            index_list=[]
            vec_tc2=[]
            order_tc2=[]
            i_list=[]
            order_i=[]
            calc_pos=False

            #Verifying if there were any crossings taking place
            for i in range(1,self.N+2*self.n_i):
                if(self.pos[k][i-1]>self.pos[k][i]):
                    index_list.append(i)

            if(k==505):
                print(index_list)

            #If there were, apply corrections
            for i in range(1,self.N+2*self.n_i):
                if(self.pos[k][i-1]>self.pos[k][i]):

                    #Conditions to calculate positions and velocities or save values for tc2
                    if(index_list.index(i)<len(index_list)-1):
                        if(i-index_list[index_list.index(i)+1]==0):
                            calc_pos=False

                        if(i-index_list[index_list.index(i)+1]<1):
                            calc_pos=True
                    else:
                        calc_pos=True #corresponds to the last collision, so we can't save anymore and need to calculate



                    ###Computations done everytime this loop is entered

                    #Computation of tc1
                    tc1 = self.a*(self.pos[k-1][i]-self.pos[k-1][i-1])/(self.pos[k-1][i]-self.pos[k-1][i-1] + self.pos[k][i-1]-self.pos[k][i])
                    # Position at t + tc1
                    pos_tc1_im1 = self.pos[k-1][i-1]+self.vel[k-1][i-1]*math.sin(tc1)-(self.pos[k-1][i-1]-self.eq_pos[i-1])*(1-math.cos(tc1))## i-1 position
                    pos_tc1_i = self.pos[k-1][i]+self.vel[k-1][i]*math.sin(tc1)-(self.pos[k-1][i]-self.eq_pos[i])*(1-math.cos(tc1))## i position

                    #Computation of tc2
                    tc2 = tc1*(self.pos[k-1][i]-self.pos[k-1][i-1])/(self.pos[k-1][i]-self.pos[k-1][i-1] + pos_tc1_im1-pos_tc1_i)
                    vec_tc2.append(tc2)
                    i_list.append(i)


                    if(calc_pos==True):

                        order_i=[x for _,x in sorted(zip(vec_tc2,i_list))]

                        vec_tc2.sort()
                        
                        for j in range(len(order_i)):

                            #now tc2->order_tc2[j] and i->order_i[j]

                            # Positions at t + tc2
                            pos_tc2_im1 = self.pos[k-1][order_i[j]-1]+self.vel[k-1][order_i[j]-1]*math.sin(vec_tc2[j])-(self.pos[k-1][order_i[j]-1]-self.eq_pos[order_i[j]-1])*(1-math.cos(vec_tc2[j]))## i-1 position
                            pos_tc2_i = self.pos[k-1][order_i[j]]+self.vel[k-1][order_i[j]]*math.sin(vec_tc2[j])-(self.pos[k-1][order_i[j]]-self.eq_pos[order_i[j]])*(1-math.cos(vec_tc2[j]))## i position
                            #velocities at t+tc2
                            vel_tc2_im1=self.vel[k-1][order_i[j]-1]*math.cos(vec_tc2[j])-math.sin(vec_tc2[j])*(self.pos[k-1][order_i[j]-1]-self.eq_pos[order_i[j]-1])## i-1 velocity
                            vel_tc2_i=self.vel[k-1][order_i[j]]*math.cos(vec_tc2[j])-math.sin(vec_tc2[j])*(self.pos[k-1][order_i[j]]-self.eq_pos[order_i[j]])## i velocity

                    

                            # Evolution from tc2 to t+dt
                            #Positions
                            self.pos[k][order_i[j]-1]=pos_tc2_im1+vel_tc2_im1*math.sin(self.a-vec_tc2[j])-(pos_tc2_im1-self.eq_pos[order_i[j]-1])*(1-math.cos(self.a-vec_tc2[j]))+(self.a-vec_tc2[j])*(self.a-vec_tc2[j])/2*self.delta
                            self.pos[k][order_i[j]]=pos_tc2_i+vel_tc2_i*math.sin(self.a-vec_tc2[j])-(pos_tc2_i-self.eq_pos[order_i[j]])*(1-math.cos(self.a-vec_tc2[j]))-(self.a-vec_tc2[j])*(self.a-vec_tc2[j])/2*self.delta


                            #Velocities
                            self.vel[k][order_i[j]-1]=vel_tc2_im1*math.cos(self.a-vec_tc2[j])-math.sin(self.a-vec_tc2[j])*(pos_tc2_im1-self.eq_pos[order_i[j]-1])+(self.a-vec_tc2[j])*self.delta
                            self.vel[k][order_i[j]]=vel_tc2_i*math.cos(self.a-vec_tc2[j])-math.sin(self.a-vec_tc2[j])*(pos_tc2_i-self.eq_pos[order_i[j]])-(self.a-vec_tc2[j])*self.delta

                            #print(k,order_i[j],self.vel[k][order_i[j]-1])


                            calc_pos=False#after calculating the positions and velocities ordered according to tc2 search for more cases

                            #Clear the lists
                            del vec_tc2[:] 
                            del i_list[:]


            #print(vec_tc2)

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
        pos_list=[[self.pos[i][l] for i in range(lm,rm)] for l in range(self.N)]
        for i in range(self.N):
            plt.plot(t_list,pos_list[i])


        #plt.axhline(y=self.delta*(self.N-0.5),xmin=0,xmax=self.num_iter,c="black")
        #plt.axhline(y=-self.delta/2,xmin=0,xmax=self.num_iter,c="black")
        #plt.plot(t_list,pos_list1,t_list,pos_list2,t_list,pos_list3,t_list,pos_list4,t_list,pos_list5,t_list,pos_list6,t_list,pos_list7,t_list,pos_list8,t_list,self.energy)
        


        self.energy[0]=self.energy[1]
        plt.plot(self.energy)

        plt.xlabel('iters)')
        plt.ylabel('energy')
        plt.title('About as simple as it gets, folks')
        plt.grid(True)
        #plt.savefig("test.png")
        plt.show()

        #Histogram of the velocity distribution
        #vlist=[self.vel[self.num_iter-1][i] for i in range(self.n_i,self.n_i+self.N)]
        #plt.hist(vlist,normed=False,bins=20)#,range=[-0.95,0.95])

        #plt.show()
        



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

        anim = animation.FuncAnimation(sheet_figure, run_animation, frames=self.num_iter, interval=250)


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


    """
    def particle_pos(self,k):
        pos_vec=[]
        for i in range(2*self.n_i+self.N):
            if(-self.delta/2<self.pos[k][i]<self.delta*(self.N-0.5)):
                pos_vec.append(self.pos[k][i])
        return pos_vec

    def particle_vel(self,k):
        vel_vec=[]
        for i in range(2*self.n_i+self.N):
            if(-self.delta/2<self.pos[k][i]<self.delta*(self.N-0.5)):
                vel_vec.append(self.vel[k][i])
        return vel_vec

    
    def particle_eq_pos(self,k):
        equi_pos=[]
        for i in range(2*self.n_i+self.N):
            if(-self.delta/2<self.pos[k][i]<self.delta*(self.N-0.5)):
                equi_pos.append((i-self.n_i)*self.delta)
        return equi_pos
    """

                        
                        
                        
                        
            

            
                
                    
