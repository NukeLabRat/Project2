import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA
plt.rcParams['axes.linewidth'] = 2




g = 40
u = np.zeros((g+1, g))
un=np.zeros((g,g))
vn=np.zeros((g,g))
v = np.zeros((g, g+1))
p = np.zeros((g+1, g+1))
dx = float(1.0/(g-1))
dy = float(1.0/(g-1))

Re = 1.0
#Time step selection
fac=(1.0/(0.25*Re))
fac2=(0.25*Re*dx**2)
if fac>fac2:
	dt=fac2
elif fac<=fac2:
	dt=fac


#dt = 0.01

p[:,:]=1.0
w = 1.5
C=np.copy(p)


F=np.zeros((g+1,g))
G=np.zeros((g,g+1))
error=1.0

t=0
Tf=1000.0
err=1
while (t<=Tf) and (err>10**-5):
	
    #print t
    #Boundary Values
    u[:,0]=0
    u[:,-1]=0
    v[:,0]=-v[:,1]
    v[:,-1]=-v[:,-2]
    v[0,:]=0
    v[-1,:]=0
    u[-1,:]=2-u[-2,:]
    u[0,:]=-u[1,:]
    speed=np.sqrt(un**2+vn**2)
    
    
    
    
    
    #print u
   
    for j in range(1,g):
        for i in range(1,g-1):   
        	F[j,i]=u[j,i]+((dt/(Re*dx**2))*(u[j,i+1]-2*u[j,i]+u[j,i-1]))\
                +(((dt/(Re*dy**2)))*(u[j+1,i]-2*u[j,i]+u[j-1,i]))\
                -(dt/dx)*(((u[j,i]+u[j,i+1])/2)**2-((u[j,i-1]+u[j,i])/2)**2)\
                -(dt/dy)*(((v[j,i]+v[j,i+1])/2)*((u[j,i]+u[j+1,i])/2)-(((v[j-1,i]+v[j,i+1])/2)*((u[j-1,i]+u[j,i])/2)))
    #print u[-2,:]
    F[-1,:]=2-F[-2,:]
    F[0,:]=-F[-1,:]
    F[:,0]=0.0
    F[:,-1]=0.0
    #print F   

    for j in range(1,g-1):
        for i in range(1,g):
        	G[j,i]=v[j,i]+((dt/(Re*dx**2))*(v[j,i+1]-2*v[j,i]+v[j,i-1]))\
                +(((dt/(Re*dy**2)))*(v[j+1,i]-2*v[j,i]+v[j-1,i]))\
                -(dt/dy)*(((v[j,i]+v[j+1,i])/2)**2-((v[j-1,i]+v[j,i])/2)**2)\
                -(dt/dx)*((((u[j,i]+u[j+1,i])/2)*((v[j,i]+v[j,i+1])/2))-(((u[j,i-1]+u[j+1,i-1])/2)*((v[j,i-1]+v[j,i])/2)))
    G[0,:]=0.0
    G[-1,:]=0.0
    G[:,0]=-G[:,-1]
    G[:,-1]=-G[:,-2]
    #C=[]
    temp=1
 
  
    p2 = np.copy(p)
    #print type(F)
    #print type(G)
    error = 1
    while error > 10**-5:
        p[:,0] = p[:,1]
        p[:,-1] = p[:,-2]
        p[0,:] = p[1,:]
        p[-1,:] = p[-2,:]  
        p1 = np.copy(p)            
        for i in range(1,g):
            for j in range(1,g):
                C[j,i]=((F[j,i]-F[j,i-1])/(dx*dt))+((G[j,i]-G[j-1,i])/(dy*dt))
                p2[j,i] = (1-w)*p[j,i]+w*(((dy**2)*(p[j,i-1]+p[j,i+1])+(dx**2)*(p[j-1,i]+p[j+1,i])-(dx**2)*(dy**2)*C[j,i])/(2*(dx**2+dy**2)))

                p[j,i] = p2[j,i]
                #print const
        error = LA.norm(p1-p)
    for j in range(1,g):
        for i in range(1,g-1):
            u[j,i]=F[j,i]-(dt/dx)*(p[j,i+1]-p[j,i])

    for j in range(1,g-1):
        for i in range(1,g):
            v[j,i]=G[j,i]-(dt/dx)*(p[j+1,i]-p[j,i])   
    t=t+dt
   


    for j in range(g):
    	for i in range(g):
    		un[j,i]=(u[j,i]+u[j+1,i])/2
    for j in range(g):
    	for i in range(g):
    		vn[j,i]=((v[j,i]+v[j,i+1]))/2

    err=LA.norm(np.sqrt(un**2+vn**2)-speed)
    print("Error",err,"Time", t)






X1=np.linspace(0,1,np.size(un[:,0]))
Y1=np.linspace(0,1,np.size(un[:,0]))

Y,X=np.meshgrid(Y1,X1)
speed=np.sqrt(un*un+vn*vn)

plt.streamplot(X1,Y1,un,vn,density=1,  color=un, linewidth=2, cmap=plt.cm.autumn)
plt.xlim(0,1)
plt.ylim(0,1)
plt.savefig('9.png',format='png',dpi=300)
plt.show()
plt.pcolor(un)
plt.xlim(0,g)
plt.ylim(0,g)
plt.savefig('10.png',format='png',dpi=300)
plt.colorbar()
plt.show()
plt.pcolor(vn)
plt.xlim(0,g)
plt.ylim(0,g)
plt.savefig('11.png',format='png',dpi=300)
plt.show()
np.savetxt('U2.txt',un)
np.savetxt('V2.txt',vn)
