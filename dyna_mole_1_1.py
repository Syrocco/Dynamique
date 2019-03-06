import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def distance(x1,x2,y1,y2,z1,z2):
	return np.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)

def Ecinetique(m,vx,vy,vz):
	return 0.5*m*(vx**2+vy**2+vz**2)

def Epotentielle(m1,m2,d):
	G = 6.67428*10**-11
	return -0.5*G*m1*m2/d   #Multiplication par 0.5 parceque je compte deux fois l'energie: pour un couple {i,j} de particule, je compte Eij et Eji, donc il faut diviser par deux





nombrePlan=5
Tmasse=[10**20,10**5,10**5,10**5,10**5]

G = 6.67428e-11

temps=10
dt=0.001
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

TPosx[0,1]=0
TPosy[0,1]=-1350
TPosz[0,1]=0
TVitx[0,1]=150
TVity[0,1]=0
TVitz[0,1]=1500


TPosx[0,2]=0
TPosy[0,2]=1500
TPosz[0,2]=0
TVitx[0,2]=1500
TVity[0,2]=0
TVitz[0,2]=0


TPosx[0,3]=0
TPosy[0,3]=0
TPosz[0,3]=1500
TVitx[0,3]=-1500
TVity[0,3]=0
TVitz[0,3]=500



TPosx[0,4]=500
TPosy[0,4]=150
TPosz[0,4]=550
TVitx[0,4]=-250
TVity[0,4]=2800
TVitz[0,4]=500



"""TPosx[0,5]=500
TPosy[0,5]=-500
TPosz[0,5]=-500
TVitx[0,5]=500
TVity[0,5]=0
TVitz[0,5]=-1000"""



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
Etot=Epot+Ecin
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(TPosx[:,0],TPosy[:,0],TPosz[:,0],label="0")
ax.plot(TPosx[:,1],TPosy[:,1],TPosz[:,1],label="1")
ax.plot(TPosx[:,2],TPosy[:,2],TPosy[:,2],label="2")
ax.plot(TPosx[:,3],TPosy[:,3],TPosz[:,3],label="3")
ax.plot(TPosx[:,4],TPosy[:,4],TPosz[:,4],label="4")
"""ax.plot(TPosx[:,5],TPosy[:,5],TPosy[:,5],label="5")"""
plt.legend()
plt.figure()
plt.plot(ttab,Epot)
plt.plot(ttab,Ecin)
plt.plot(ttab,Etot)
plt.show()
