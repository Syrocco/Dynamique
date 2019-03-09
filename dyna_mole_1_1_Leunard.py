import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def distance(x1,x2,y1,y2,z1,z2):
	return np.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)







nombrePlan=2


G = 6.67428e-11

temps=100
dt=0.01
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


TPosx[0,0]=-2
TPosy[0,0]=2
TPosz[0,0]=0
TVitx[0,0]=0.01
TVity[0,0]=0
TVitz[0,0]=0

TPosx[0,1]=-2
TPosy[0,1]=-2
TPosz[0,1]=0
TVitx[0,1]=0.01
TVity[0,1]=0
TVitz[0,1]=0



for i in range(1,N):
    for planete in range(nombrePlan):
        ax=0
        ay=0
        az=0
        for planeteAutre in range(nombrePlan):
            if planeteAutre!=planete:
                d=distance(TPosx[i-1,planete],TPosx[i-1,planeteAutre],TPosy[i-1,planete],TPosy[i-1,planeteAutre],TPosz[i-1,planete],TPosz[i-1,planeteAutre])
                a=((1/d**12)-(1/d**6))
                ax+=a*(TPosx[i-1,planete]-TPosx[i-1,planeteAutre])/d               
                ay+=a*(TPosy[i-1,planete]-TPosy[i-1,planeteAutre])/d
                az+=a*(TPosz[i-1,planete]-TPosz[i-1,planeteAutre])/d


        TVitx[i,planete]=TVitx[i-1,planete]+dt*ax
        TVity[i,planete]=TVity[i-1,planete]+dt*ay	
        TVitz[i,planete]=TVitz[i-1,planete]+dt*az
        TPosx[i,planete]=TPosx[i-1,planete]+dt*TVitx[i,planete]
        TPosy[i,planete]=TPosy[i-1,planete]+dt*TVity[i,planete]
        TPosz[i,planete]=TPosz[i-1,planete]+dt*TVitz[i,planete]

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(TPosx[:,0],TPosy[:,0],TPosz[:,0],label="0")
ax.plot(TPosx[:,1],TPosy[:,1],TPosz[:,1],label="1")
"""ax.plot(TPosx[:,5],TPosy[:,5],TPosy[:,5],label="5")"""
plt.legend()
plt.show()
