

N=10 # nr of sheets 

x=[] #array with the positions of sheets at each time
x.append([])
x[0] = [0 for i in range(N)] #0 aqui corresponde a t=0

vx=[] #array with the velocities of sheets at each time
vx.append([])
vx[0] = [0 for i in range(N)] #0 aqui corresponde a t=0

xeq=[] #equilibrium positions
xeq=[0 for i in range(N)]

#Parameters
a=0.05 #wpe*dt
#dt
#L #com L e N tiro delta 


#Initializing positions and velocities
pos=0
for i in range(N):
    print(i)
    x[0][i]=pos
    pos+=0.5
    print(x[0][i])

