import math
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import scipy.stats as stats
import pickle

class Problem:
    def __init__(self, number,initial_positions, initial_velocities,wpdt,delta1,n_iter,n_image,period):

        self.a=wpdt #adimensional measure of time
        self.delta=delta1 #in units of the Debye length
        self.period=period
        self.init_pos=sorted(initial_positions) #in units of Debye length
        self.init_vel=[x for _,x in sorted(zip(initial_positions,initial_velocities))] #in units of thermal velocity
        self.N=number #number of particles
        self.num_iter=n_iter #number of iterations

        self.n_i=0#number of images on each side
        self.vec_n_i_left=[]#Important for calculating the actual equilibrium positions of the particles inside the box for every k
        self.vec_n_i_left.append(0)#we only start calculating the energy for k=1


        #self.pos=[[0] * (self.N) for i in range(self.num_iter)] #final matrix of positions 
        #self.vel=[[0] * (self.N) for i in range(self.num_iter)] #final matrix of velocities
        self.eq_pos=[self.delta*n for n in range(self.N)]
        

        self.energy=[0 for n in range(n_iter)]

    def start(self):
        "num_iter -> number of iterations; deltat-> iteration time"
        self.prepos=self.init_pos
        self.prevel=self.init_vel
        self.pospos=[0 for i in range(len(self.prepos))]
        self.posvel=[0 for i in range(len(self.prepos))]
        
        argument=self.a
        a=self.a
        cos=math.cos(argument)
        sin=math.sin(argument)


        img_count_left=0
        img_count_right=0
        img_count_max=0
            
        for k in range(1,self.num_iter):

            for i in range(len(self.pospos)):
                del self.pospos[0]

            for i in range(len(self.posvel)):
                del self.posvel[0]

            self.pospos = [0 for i in range(self.N)]
            self.posvel = [0 for i in range(self.N)]


            #print(k)
            left_i=0
            right_i=0
            preeq_pos=[self.eq_pos[i] for i in range(len(self.eq_pos))]
            #print(self.prepos)


            #print("\n begin ")
            #print(k-1,self.prepos)

            #if(k>141):

            #Check which particles are inside the box 
            for i in range(len(self.prepos)):
                if(self.prepos[i]>-self.delta/2):
                    left_i=i
                    break
            for i in range(len(self.prepos)):
                if(self.prepos[i]>(self.N-0.5)*self.delta):
                    right_i=i-1
                    break
                if(i==len(self.prepos)-1):#if nothing is on the right of the box
                    right_i=len(self.prepos)-1

            

            #Clear the previous particles outside the box
         
            for i in range(left_i):
                del self.prepos[0]
                del self.prevel[0]
                del self.eq_pos[0]


            #After deleting the previous entries, right_i corresponds to N-1
            for i in range(self.N,len(self.prepos)):
                del self.prepos[self.N]
                del self.prevel[self.N]
                del self.eq_pos[self.N]
            
        
            #print(k,self.eq_pos,left_i,right_i)
            if (k-1)%self.period==0:
                pickle.dump([self.prepos,self.eq_pos], open('./lixo/'+str(k-1)+'pos'+'.pkl', 'wb'))
                pickle.dump(self.prevel, open('./lixo/'+str(k-1)+'vel'+'.pkl', 'wb'))


            
            kinetic=sum([vels*vels for vels in self.prevel])
            diff=[self.prepos[i]-self.eq_pos[i] for i in range(len(self.prepos))]
            potential=sum([diffs*diffs for diffs in diff])
            self.energy[k-1]=kinetic+potential


            img_count_left=0
            img_count_right=0
            img_count_max=0


            #print(self.prevel)


            in1=True
            in2=True
            for i in range(self.N):
                    
                self.posvel[i]=self.prevel[i]*cos-sin*(self.prepos[i]-self.eq_pos[i])
                self.pospos[i]=self.prepos[i]+self.prevel[i]*sin-(self.prepos[i]-self.eq_pos[i])*(1-cos)


                if(self.pospos[i]>(self.N-0.5)*self.delta and in1==True):
                    img_count_left=self.N-i
                    in1=False


                if(self.pospos[self.N-i-1]<-0.5*self.delta and in2==True):
                    img_count_right=self.N-i
                    in2=False
                
              #print(k,self.pospos,self.eq_pos)

            img_count_max=max(img_count_left,img_count_right)


            #Create new images
            self.n_i=img_count_max
            

            #print(self.prepos)

            #print(k,self.pospos,self.eq_pos)

            #Creation of images
            for i in range(self.n_i):
                #Right images
                self.posvel.append(self.posvel[i])
                self.pospos.append(self.pospos[i]+self.N*self.delta)
                self.eq_pos.append(self.eq_pos[i]+self.N*self.delta)
                #We want to guarantee that we have the same number of images for k and k-1. This is needed for the computation of the crossing. So
                self.prevel.append(self.prevel[i])
                self.prepos.append(self.prepos[i]+self.N*self.delta)

            

            for i in range(self.n_i):
                #Left images
                self.posvel.insert(0,self.posvel[self.N-1])
                self.pospos.insert(0,self.pospos[self.N-1]-self.N*self.delta)
                self.prevel.insert(0,self.prevel[self.N-1])
                self.prepos.insert(0,self.prepos[self.N-1]-self.N*self.delta)
                self.eq_pos.insert(0,self.eq_pos[self.N-1]-self.N*self.delta)


            #print(k,self.pospos,self.eq_pos)
            #print(k,self.pospos)


            #Verifying if there were any crossings taking place
            for i in range(1,self.N+2*self.n_i):

                if(self.pospos[i-1]>self.pospos[i]):

                   #Computation of tc1
                    tc1 = self.a*(self.prepos[i]-self.prepos[i-1])/(self.prepos[i]-self.prepos[i-1] + self.pospos[i-1]-self.pospos[i])
                    # Position at t + tc1
                    pos_tc1_im1 = self.prepos[i-1]+self.prevel[i-1]*math.sin(tc1)-(self.prepos[i-1]-self.eq_pos[i-1])*(1-math.cos(tc1))## i-1 position
                    pos_tc1_i = self.prepos[i]+self.prevel[i]*math.sin(tc1)-(self.prepos[i]-self.eq_pos[i])*(1-math.cos(tc1))## i position


                    #Computation of tc2
                    tc2 = tc1*(self.prepos[i]-self.prepos[i-1])/(self.prepos[i]-self.prepos[i-1] + pos_tc1_im1-pos_tc1_i)
                    #tc2= -math.atan((self.prepos[i]-self.prepos[i-1])/(self.prevel[i]-self.prevel[i-1]))


                    # Positions at t + tc2
                    pos_tc2_im1 = self.prepos[i-1]+self.prevel[i-1]*math.sin(tc2)-(self.prepos[i-1]-self.eq_pos[i-1])*(1-math.cos(tc2))## i-1 position
                    pos_tc2_i = self.prepos[i]+self.prevel[i]*math.sin(tc2)-(self.prepos[i]-self.eq_pos[i])*(1-math.cos(tc2))## i position
                    #velocities at t+tc2
                    vel_tc2_im1=self.prevel[i-1]*math.cos(tc2)-math.sin(tc2)*(self.prepos[i-1]-self.eq_pos[i-1])## i-1 velocity
                    vel_tc2_i=self.prevel[i]*math.cos(tc2)-math.sin(tc2)*(self.prepos[i]-self.eq_pos[i])## i velocity



                    # Evolution from tc2 to t+dt
                    #Positions
                    self.pospos[i-1]=pos_tc2_im1+vel_tc2_im1*math.sin(self.a-tc2)-(pos_tc2_im1-self.eq_pos[i-1])*(1-math.cos(self.a-tc2))+(self.a-tc2)*(self.a-tc2)/2*self.delta
                    self.pospos[i]=pos_tc2_i+vel_tc2_i*math.sin(self.a-tc2)-(pos_tc2_i-self.eq_pos[i])*(1-math.cos(self.a-tc2))-(self.a-tc2)*(self.a-tc2)/2*self.delta
                    #Velocities
                    self.posvel[i-1]=vel_tc2_im1*math.cos(self.a-tc2)-math.sin(self.a-tc2)*(pos_tc2_im1-self.eq_pos[i-1])+(self.a-tc2)*self.delta
                    self.posvel[i]=vel_tc2_i*math.cos(self.a-tc2)-math.sin(self.a-tc2)*(pos_tc2_i-self.eq_pos[i])-(self.a-tc2)*self.delta


            #Order all the particles according to their positions after each time interval dt
            self.posvel=[x for _,x in sorted(zip(self.pospos,self.posvel))]
            self.pospos.sort()


            
            #Clear the auxiliary k-1 entries used to calculate crossings
            for i in range(img_count_max):
                del self.prepos[0]
                del self.prevel[0]

            for i in range(self.N,self.N+img_count_max):
                del self.prepos[self.N]
                del self.prevel[self.N]


            self.prevel=[self.posvel[i] for i in range(len(self.posvel))]
            self.prepos=[self.pospos[i] for i in range(len(self.pospos))]



    
    def energy_from_velocity(self,vel):
        return sum([vels*vels for vels in vel])

    def energy_from_position(self,pos,pos_eq):
        diff=[pos[i]-pos_eq[i] for i in range(len(pos))]
        return sum([diffs*diffs for diffs in diff])



    def gettrajectory(self,particle,period):
        lm=0
        rm=self.num_iter-1
        t_list=[i for i in range(lm ,rm) if i%period==0]
        traj=[]
        for i in range(lm,rm):
            if i%period==0:
                pkl_file=open('./lixo/'+str(i)+'pos'+'.pkl','rb')
                data=pickle.load(pkl_file)
                traj.append(data[0][particle])
                pkl_file.close()
        return traj
        
    def getalltrajectories(self,period):
        lm=0
        rm=self.num_iter-1
        traj=[]
        t_list=[i for i in range(lm,rm) if i%period==0]
        for i in range(lm,rm):
            if i%period==0:
                pkl_file=open('./lixo/'+str(i)+'pos'+'.pkl','rb')
                data=pickle.load(pkl_file)
                traj.append(data[0])
                pkl_file.close()
        return traj

    def getallenergy(self,period):
        lm=0
        rm=self.num_iter-1
        energy=[]
        t_list=[i for i in range(lm, rm) if i%period==0]
        for i in range(lm,rm):
            if i%period==0:
                pkl_file_pos=open('./lixo/'+str(i)+'pos'+'.pkl','rb')
                pkl_file_vel=open('./lixo/'+str(i)+'vel'+'.pkl','rb')
                data_pos=pickle.load(pkl_file_pos)
                data_vel=pickle.load(pkl_file_vel)
                energy_vel=self.energy_from_velocity(data_vel)
                energy_pos=self.energy_from_position(data_pos[0],data_pos[1])
                energy.append(energy_vel+energy_pos)
                pkl_file_pos.close()
                pkl_file_vel.close()

        return energy

    def getvelocities(self,k):
        pkl_file=open('./lixo/'+str(k)+'vel.pkl','rb')
        dat=pickle.load(pkl_file)
        data=[dat[0][i] for i in range(len(dat))]
        pkl_file.close()
        return data

 


    def cenas(self):        
        #energy conservation?
        print("oi")
        self.energy[self.num_iter-1]=self.energy[self.num_iter-2]
        plt.plot(self.energy)

        plt.xlabel('iters)')
        plt.ylabel('energy')
        plt.title('About as simple as it gets, folks')
        plt.grid(True)
        #plt.savefig("test.png")
        plt.show()
        
"""
        #Histogram of the velocity distribution
        #vlist=[self.vel[self.num_iter-1][i] for i in range(self.N)]
        #plt.hist(vlist,normed=False,bins=10)#,range=[-0.95,0.95])

        #plt.show()
        



        ######Animation with sheets####################
        sheet_figure = plt.figure()
        ax = plt.axes(xlim=(-self.delta, (self.N)*self.delta), ylim=(-0.1, 0.1))
        point_set, = ax.plot([self.pos[0][i] for i in range(self.N)],[0 for i in range(self.N)], 'ro')
        plt.axvline(x=-self.delta/2)
        plt.axvline(x=(self.N-0.5)*self.delta)
        pt_circ = plt.Circle((4, 4), 1, color='b', fill=False)
        ax.add_artist(pt_circ)
        ttl = ax.text(.5, 1.005, '', transform = ax.transAxes)

        #I'm printing only the particles and not the images
        def run_animation(i):
            point_set.set_data([self.pos[i][k] for k in range(self.N)],[0 for k in range(self.N)])
            ttl.set_text(str(i))    
            return point_set,

        anim = animation.FuncAnimation(sheet_figure, run_animation, frames=self.num_iter, interval=20)


        #plt.show()


    def energy_from_velocity(self,vel):
        return sum([vels*vels for vels in vel])

    def energy_from_position(self,pos,pos_eq):
        diff=[pos[i]-pos_eq[i] for i in range(len(pos))]
        return sum([diffs*diffs for diffs in diff])

    def momentum(self,vel):
        return sum (vel)



    def particle_pos(self,k):
        pos_vec=[]
        num=0
        for i in range(len(self.pospos)):
            if(-self.delta/2<self.pospos[i]<self.delta*(self.N-0.5)):
                pos_vec.append(self.pospos[i])
                num+=1

            #print(k,self.pospos,pos_vec,num)

        return pos_vec

    def particle_vel(self,k):
        vel_vec=[]
        for i in range(len(self.pospos)):
            if(-self.delta/2<self.pospos[i]<self.delta*(self.N-0.5)):
                vel_vec.append(self.posvel[i])

        return vel_vec

    
    def particle_eq_pos(self,k):
        equi_pos=[]
        for i in range(len(self.pospos)):
            if(-self.delta/2<self.pospos[i]<self.delta*(self.N-0.5)):
                #I will have a different number of left images for every k
                equi_pos.append((i-self.vec_n_i_left[k])*self.delta)
        
        return equi_pos

    """
                        
                        
                        
            

            
                
                    
