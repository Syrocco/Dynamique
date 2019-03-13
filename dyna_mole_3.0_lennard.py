import numpy as np
import matplotlib.pyplot as plt
import numpy.random as rd
from mpl_toolkits.mplot3d import Axes3D
import time




def norme(Vec3D):
    return np.sqrt(Vec3D[0]**2+Vec3D[1]**2+Vec3D[2]**2)


def GenVitesse(T):
    vit_xyz=np.zeros(3)
    for i in range(3):
        vit_xyz[i]=rd.random()
    NormEtVit=norme(vit_xyz)/rd.normal(T,0)   #constante permettant de normer vit_xyz pour ensuite lui donner une norme prise au hasard dans une gaussienne

    for i in range(3):
        signe=rd.random()
        if signe<0.5:
            vit_xyz[i]=vit_xyz[i]/NormEtVit
        else:
            vit_xyz[i]=-vit_xyz[i]/NormEtVit
    
    return np.array(vit_xyz)


def GenPosition(nombrePlanete,rayon,methode):
    if methode=="Cube":
        Ecart=(2*rayon)/((nombrePlanete)**(1/3))
        T=[]
        for i in np.arange(-rayon,rayon,Ecart):
            for j in np.arange(-rayon+Ecart/2,rayon,Ecart):
                for z in np.arange(-rayon,rayon,Ecart):
                    T.append([i,j,z])
        if len(T)<nombrePlanete:
            print("L'attribution est mauvaise")
            return "L'attribution est mauvaise"
        return np.array(T[:nombrePlanete])   #Je ne suis pas arrivé à faire un programme renvoyant tout le temps une matrice avec le nombre de planete exacte, du coup, j'augmente un peu les
                                   #limites du cube et je ne prend que les "nombrePlanetes" premières valeurs de la distribution.
                                  
def AttributionInitiale():
	particule=GenPosition(nombrePlan, 1, "Cube")
	for i in range(nombrePlan):
  	  TPosx[0,i]=particule[i,0]
  	  TPosy[0,i]=particule[i,1]
  	  TPosz[0,i]=particule[i,2]

	for i in range(nombrePlan):
		A=GenVitesse(0.2)
		TVitx[0,i]=A[0]
		TVity[0,i]=A[1]
		TVitz[0,i]=A[2]


def CalculAcceleration():
	d=distance(TPosx[i-1,planete],TPosx[i-1,planeteAutre],TPosy[i-1,planete],TPosy[i-1,planeteAutre],TPosz[i-1,planete],TPosz[i-1,planeteAutre])
	a=((12/d**13)-(6/d**7))/1
	return a*(TPosx[i-1,planete]-TPosx[i-1,planeteAutre])/d, a*(TPosy[i-1,planete]-TPosy[i-1,planeteAutre])/d,a*(TPosz[i-1,planete]-TPosz[i-1,planeteAutre])/d, d 

def CalculVitesse():
    TVitx[i,planete]=TVitx[i-1,planete]+dt*ax
    TVity[i,planete]=TVity[i-1,planete]+dt*ay	
    TVitz[i,planete]=TVitz[i-1,planete]+dt*az
    TPosx[i,planete]=TPosx[i-1,planete]+dt*TVitx[i,planete]
    TPosy[i,planete]=TPosy[i-1,planete]+dt*TVity[i,planete]
    TPosz[i,planete]=TPosz[i-1,planete]+dt*TVitz[i,planete]                

def distance(x1,x2,y1,y2,z1,z2):
	return np.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)

def Ecinetique(vx,vy,vz):
	return 0.5*1*(vx**2+vy**2+vz**2)

def Epotentielle(d):
	return 0.5*((1/d**12)-(1/d**6))   #Multiplication par 0.5 parceque je compte deux fois l'energie: pour un couple {i,j} de particule, je compte Eij et Eji, donc il faut diviser par deux



rd.seed(4)


nombrePlan=3

temps=200
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
ttab=np.arange(dt,temps,dt)#ATTENTION AU TEMPS MODULAIRE


AttributionInitiale()
  




print("Début des calculs")
for i in range(1,N):
    for planete in range(nombrePlan):
        ax=0
        ay=0
        az=0
        for planeteAutre in range(nombrePlan):
            if planeteAutre!=planete:      
                a=CalculAcceleration()
                ax+=a[0]
                ay+=a[1]
                az+=a[2]
                Epot[i-1]=Epot[i-1]+Epotentielle(a[3])
        CalculVitesse()
        Ecin[i-1]=Ecin[i-1]+Ecinetique(TVitx[i,planete],TVity[i,planete],TVitz[i,planete])        
        
       
        
        
     
fig = plt.figure()
ax = fig.gca(projection='3d')
for i in range(nombrePlan):
    ax.plot(TPosx[:,i],TPosy[:,i],TPosz[:,i],label=str(i))
plt.legend()
plt.figure()
Etot=Epot+Ecin
plt.plot(ttab,Epot,label="Epot")
plt.plot(ttab,Ecin,label="Ecin")
plt.plot(ttab,Etot,label="Etot")
plt.legend()
plt.show()
