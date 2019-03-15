import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def distance(x1,x2,y1,y2,z1,z2):
	return np.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)

def Ecinetique(m,vx,vy,vz):
	return 0.5*m*(vx**2+vy**2+vz**2)

def Epotentielle(m1,m2,d):
	return -0.5*G*m1*m2/d   #Multiplication par 0.5 parceque je compte deux fois l'energie: pour un couple {i,j} de particule, je compte Eij et Eji, donc il faut diviser par deux


nombrePlan=4
Tmasse=[1.988544e30,3.302e23,4.8685e24, 5.97219e24]

AU = 149.6e9
G = 6.67428e-11*( (24. * 60. * 60.) ** 2)/AU**3

temps=1000
dt=0.1
N=int(temps/dt)

T=np.zeros(N)

TPosx=np.zeros((N,nombrePlan))
TPosy=np.zeros((N,nombrePlan))
TPosz=np.zeros((N,nombrePlan))
TVitx=np.zeros((N,nombrePlan))
TVity=np.zeros((N,nombrePlan))
TVitz=np.zeros((N,nombrePlan))
Epot=np.zeros(N-1)
Ecin=np.zeros(N-1)


TPosx[0,0]=0
TPosy[0,0]=0
TPosz[0,0]=0
TVitx[0,0]=0
TVity[0,0]=0
TVitz[0,0]=0

TPosx[0,1]=3.373130583178054e-01
TPosy[0,1]=-2.035999304350679e-01
TPosz[0,1]=-4.758290141863730e-02
TVitx[0,1]=9.017706724241308e-03
TVity[0,1]=2.537646569534078e-02
TVitz[0,1]=1.246132243839972e-03


TPosx[0,2]=5.496573190048547e-01
TPosy[0,2]=-4.760713087406738e-01
TPosz[0,2]=-3.824501332858326e-02
TVitx[0,2]=1.311279754946672e-02
TVity[0,2]=1.521238322677600e-02
TVitz[0,2]=-5.482326391309875e-04



TPosx[0,3]=-1.711651664398659e-01
TPosy[0,3]=9.682993411474208e-01
TPosz[0,3]= -2.976388969151488e-05
TVitx[0,3]=-1.721626633528917e-02
TVity[0,3]=-3.064115847218701e-03
TVitz[0,3]=6.740893165190773e-07




ttab=np.linspace(0,dt,N-1)

for i in range(1,N):
    for planete in range(nombrePlan):
        ax=0
        ay=0
        az=0
        for planeteAutre in range(nombrePlan):
            if planeteAutre!=planete:
                d=distance(TPosx[i-1,planete],TPosx[i-1,planeteAutre],TPosy[i-1,planete],TPosy[i-1,planeteAutre],TPosz[i-1,planete],TPosz[i-1,planeteAutre])
                a=-G*Tmasse[planeteAutre]/d**2
                ax+=a*(TPosx[i-1,planete]-TPosx[i-1,planeteAutre])/d               
                ay+=a*(TPosy[i-1,planete]-TPosy[i-1,planeteAutre])/d
                az+=a*(TPosz[i-1,planete]-TPosz[i-1,planeteAutre])/d
                Epot[i-1]=Epot[i-1]+Epotentielle(Tmasse[planete],Tmasse[planeteAutre],d)

        TVitx[i,planete]=TVitx[i-1,planete]+dt*ax
        TVity[i,planete]=TVity[i-1,planete]+dt*ay	
        TVitz[i,planete]=TVitz[i-1,planete]+dt*az
        TPosx[i,planete]=TPosx[i-1,planete]+dt*TVitx[i,planete]
        TPosy[i,planete]=TPosy[i-1,planete]+dt*TVity[i,planete]
        TPosz[i,planete]=TPosz[i-1,planete]+dt*TVitz[i,planete]
        
        Ecin[i-1]=Ecin[i-1]+Ecinetique(Tmasse[planete],TVitx[i,planete],TVity[i,planete],TVitz[i,planete])

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(TPosx[:,0],TPosy[:,0],TPosz[:,0],label="Soleil",color="orange")
ax.plot(TPosx[:,1],TPosy[:,1],TPosz[:,1],label="Mercure",color="red")
ax.plot(TPosx[:,2],TPosy[:,2],TPosz[:,2],label="Venus",color="grey")
ax.plot(TPosx[:,3],TPosy[:,3],TPosz[:,3],label="Terre",color="blue")
plt.legend()

ax.set_xlim3d([-1.0, 1.0])

ax.set_ylim3d([-1.0, 1.0])

ax.set_zlim3d([-1.0, 1.0])


plt.show()



