import math
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import scipy.stats as stats
import pickle

class Problem:
    def __init__(self, number,initial_positions, initial_velocities,wpdt,delta1,n_iter,period):

        self.a=wpdt #adimensional measure of time
        self.delta=delta1 #distance between sheets in units of the Debye length
        self.period=period #iteration interval between saved points
        self.init_pos=sorted(initial_positions) # initial positions of the sheets
        self.init_vel=[x for _,x in sorted(zip(initial_positions,initial_velocities))] #initial velocities of the sheets (in units of vt)
        self.N=number #number of particles
        self.num_iter=n_iter #number of iterations

        self.n_i=0#number of images on each side
        self.vec_n_i_left=[]#Important for calculating the actual equilibrium positions of the particles inside the box for every k
        self.vec_n_i_left.append(0)#we only start calculating the energy for k=1

        self.eq_pos=[self.delta*n for n in range(self.N)]#equilibrium positions of the sheets
        
        self.energy=[0 for n in range(n_iter)]#total energy of the system at each time

    def start(self):
        #The computation of the coordinates at t+dt requires their knowledge at t, dt being the chosen time step. So in each iteration we need variables to describe positions and velocities at t (prepos and prevel) and at t+dt(pospos and posvel) which are initialized below
        self.prepos=self.init_pos 
        self.prevel=self.init_vel
        self.pospos=[0 for i in range(len(self.prepos))]
        self.posvel=[0 for i in range(len(self.prepos))]
        
            
        #Start of the time loop
        for k in range(1,self.num_iter):


            #First of all clear the "pos" coordinates from the previous iteration. The "pos" coordinates from the previous iteration change to the "pre" coordinates from the current iteration at the end of the loop. The "pos" coordinates from this iteration can then be deleted in order to be computed later based on the new "pre" coordinates
            for i in range(len(self.pospos)):
                del self.pospos[0]

            for i in range(len(self.posvel)):
                del self.posvel[0]

            self.pospos = [0 for i in range(self.N)]
            self.posvel = [0 for i in range(self.N)]


            #Now we want to check which sheets are inside the box and erase all the others. Image sheets are created at each iteration (this is better explained below in the code) and here we promote the ones that enter the box to "real" sheets
            left_i=0 #variable used to check the number of sheets on the left of the box
            right_i=0 #variable used to check the number of sheets on the right of the box


            #Check which sheets are inside the box. Knowing that the equilibrium positions of the first and last sheets are, respectively, x=0 and x=(N-1)delta, the box is bounded by [-delta/2,(N-1/2)delta
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

            
            #Clear the sheets on the left of the box
            for i in range(left_i):
                del self.prepos[0]
                del self.prevel[0]
                del self.eq_pos[0]


            #Clear the sheets on the right of the box. After deleting the previous entries, right_i corresponds to N-1
            for i in range(self.N,len(self.prepos)):
                del self.prepos[self.N]
                del self.prevel[self.N]
                del self.eq_pos[self.N]

            
            #Save data at the rate "period" using pickle
            if (k-1)%self.period==0:
                pickle.dump([self.prepos,self.eq_pos], open('./lixo/'+str(k-1)+'pos'+'.pkl', 'wb'))
                pickle.dump(self.prevel, open('./lixo/'+str(k-1)+'vel'+'.pkl', 'wb'))


            #Energy diagnosis at each time interval
            kinetic=sum([vels*vels for vels in self.prevel])
            diff=[self.prepos[i]-self.eq_pos[i] for i in range(len(self.prepos))]
            potential=sum([diffs*diffs for diffs in diff])
            self.energy[k-1]=kinetic+potential


            #Initialization of auxiliar variables to count the needed image sheets to be created
            img_count_left=0
            img_count_right=0
            img_count_max=0

            #Auxiliar variables
            argument=self.a
            cos=math.cos(argument)
            sin=math.sin(argument)


            in1=True
            in2=True
            #Loop to calculate the "pos" from the "pre" coordinates
            for i in range(self.N):
                    
                self.posvel[i]=self.prevel[i]*cos-sin*(self.prepos[i]-self.eq_pos[i])
                self.pospos[i]=self.prepos[i]+self.prevel[i]*sin-(self.prepos[i]-self.eq_pos[i])*(1-cos)


                #Check which is the first sheet at the right the box after the computation (if there are some) and assign the number of left images needed to compensate for that into the variable "img_count_left"
                if(self.pospos[i]>(self.N-0.5)*self.delta and in1==True):
                    img_count_left=self.N-i
                    in1=False

                #Check which is the first sheet at the left the box and assign the number of right images needed to compensate for that into the variable "img_count_right"
                if(self.pospos[self.N-i-1]<-0.5*self.delta and in2==True):
                    img_count_right=self.N-i
                    in2=False
                

            img_count_max=max(img_count_left,img_count_right)#check the maximum between both image counts. This will be the actual number of images creates on each side of the box


            #Assigning the number of images to the corresponding class variable
            self.n_i=img_count_max


            #Creation of images

            #Right images
            for i in range(self.n_i):
                self.posvel.append(self.posvel[i])
                self.pospos.append(self.pospos[i]+self.N*self.delta)
                self.eq_pos.append(self.eq_pos[i]+self.N*self.delta)
                #We want to guarantee that we have the same number of images for k and k-1. This is needed for the computation of the crossings. So
                self.prevel.append(self.prevel[i])
                self.prepos.append(self.prepos[i]+self.N*self.delta)
            

            #Left images
            for i in range(self.n_i):
                self.posvel.insert(0,self.posvel[self.N-1])
                self.pospos.insert(0,self.pospos[self.N-1]-self.N*self.delta)
                self.prevel.insert(0,self.prevel[self.N-1])
                self.prepos.insert(0,self.prepos[self.N-1]-self.N*self.delta)
                self.eq_pos.insert(0,self.eq_pos[self.N-1]-self.N*self.delta)



            #Verifying if there were any crossings taking place
            for i in range(1,self.N+2*self.n_i):

                if(self.pospos[i-1]>self.pospos[i]):#condition for a crossing to occur

                    #This algorithm is performed based on Dawson's article from 1970 and so the correspondig notation is used (e.g: tc1, tc2,...)

                    #Computation of tc1
                    tc1 = self.a*(self.prepos[i]-self.prepos[i-1])/(self.prepos[i]-self.prepos[i-1] + self.pospos[i-1]-self.pospos[i])
                    # Positions at t + tc1
                    pos_tc1_im1 = self.prepos[i-1]+self.prevel[i-1]*math.sin(tc1)-(self.prepos[i-1]-self.eq_pos[i-1])*(1-math.cos(tc1))## i-1 position
                    pos_tc1_i = self.prepos[i]+self.prevel[i]*math.sin(tc1)-(self.prepos[i]-self.eq_pos[i])*(1-math.cos(tc1))## i position


                    #Computation of tc2
                    tc2 = tc1*(self.prepos[i]-self.prepos[i-1])/(self.prepos[i]-self.prepos[i-1] + pos_tc1_im1-pos_tc1_i)


                    # Positions at t + tc2
                    pos_tc2_im1 = self.prepos[i-1]+self.prevel[i-1]*math.sin(tc2)-(self.prepos[i-1]-self.eq_pos[i-1])*(1-math.cos(tc2))## i-1 position
                    pos_tc2_i = self.prepos[i]+self.prevel[i]*math.sin(tc2)-(self.prepos[i]-self.eq_pos[i])*(1-math.cos(tc2))## i position
                    #velocities at t+tc2
                    vel_tc2_im1=self.prevel[i-1]*math.cos(tc2)-math.sin(tc2)*(self.prepos[i-1]-self.eq_pos[i-1])## i-1 velocity
                    vel_tc2_i=self.prevel[i]*math.cos(tc2)-math.sin(tc2)*(self.prepos[i]-self.eq_pos[i])## i velocity


                    # Evolution from tc2 to t+dt. Here the constant acceleration wp*wp*delta (or delta in our units) is added or subtracted
                    #Positions
                    self.pospos[i-1]=pos_tc2_im1+vel_tc2_im1*math.sin(self.a-tc2)-(pos_tc2_im1-self.eq_pos[i-1])*(1-math.cos(self.a-tc2))+(self.a-tc2)*(self.a-tc2)/2*self.delta
                    self.pospos[i]=pos_tc2_i+vel_tc2_i*math.sin(self.a-tc2)-(pos_tc2_i-self.eq_pos[i])*(1-math.cos(self.a-tc2))-(self.a-tc2)*(self.a-tc2)/2*self.delta
                    #Velocities
                    self.posvel[i-1]=vel_tc2_im1*math.cos(self.a-tc2)-math.sin(self.a-tc2)*(pos_tc2_im1-self.eq_pos[i-1])+(self.a-tc2)*self.delta
                    self.posvel[i]=vel_tc2_i*math.cos(self.a-tc2)-math.sin(self.a-tc2)*(pos_tc2_i-self.eq_pos[i])-(self.a-tc2)*self.delta


            #Finally order all particles according to their positions 
            self.posvel=[x for _,x in sorted(zip(self.pospos,self.posvel))]#order the velocities based on the position ordering
            self.pospos.sort()


            
            #Clear the auxiliary k-1 entries used to calculate crossings
            for i in range(img_count_max):
                del self.prepos[0]
                del self.prevel[0]

            for i in range(self.N,self.N+img_count_max):
                del self.prepos[self.N]
                del self.prevel[self.N]


            #Prepare the "pre" coordinates for the next iteration which are the "pos" from the current one
            self.prevel=[self.posvel[i] for i in range(len(self.posvel))]
            self.prepos=[self.pospos[i] for i in range(len(self.pospos))]



    #Function that returns the total energy of the system at time iterations separated by "period"
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


    def getalleqpos(self,period):
        lm=0
        rm=self.num_iter-1
        traj=[]
        t_list=[i for i in range(lm,rm) if i%period==0]
        for i in range(lm,rm):
            if i%period==0:
                pkl_file=open('./lixo/'+str(i)+'pos'+'.pkl','rb')
                data=pickle.load(pkl_file)
                traj.append(data[1])
                pkl_file.close()
        return traj



    def getvelocities(self,k):
        pkl_file=open('./lixo/'+str(k)+'vel.pkl','rb')
        dat=pickle.load(pkl_file)
        data=[dat[i] for i in range(len(dat))]
        pkl_file.close()
        return data
                        
                        
            

            
                
                    
