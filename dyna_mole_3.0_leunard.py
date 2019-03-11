import numpy as np
import matplotlib.pyplot as plt
import numpy.random as rd
from mpl_toolkits.mplot3d import Axes3D




def norme(Vec3D):
    return np.sqrt(Vec3D[0]**2+Vec3D[1]**2+Vec3D[2]**2)


def AttributionVitesse(T):
    vit_xyz=np.zeros(3)
    for i in range(3):
        vit_xyz[i]=rd.random()
    NormEtVit=norme(vit_xyz)/rd.normal(T,0.005)   #constante permettant de normer vit_xyz pour ensuite lui donner une norme prise au hasard dans une gaussienne

    for i in range(3):
        signe=rd.random()
        if signe<0.5:
            vit_xyz[i]=vit_xyz[i]/NormEtVit
        else:
            vit_xyz[i]=-vit_xyz[i]/NormEtVit
    
    return np.array(vit_xyz)


def AttributionPosition(nombrePlanete,rayon,methode):
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
                                  


def distance(x1,x2,y1,y2,z1,z2):
	return np.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)

def Ecinetique(m,vx,vy,vz):
	return 0.5*m*(vx**2+vy**2+vz**2)

def Epotentielle(m1,m2,d):
	return -0.5*G*m1*m2/d   #Multiplication par 0.5 parceque je compte deux fois l'energie: pour un couple {i,j} de particule, je compte Eij et Eji, donc il faut diviser par deux






nombrePlan=20

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


particule=AttributionPosition(nombrePlan, 2, "Cube")
for i in range(nombrePlan):
    TPosx[0,i]=particule[i,0]
    TPosy[0,i]=particule[i,1]
    TPosz[0,i]=particule[i,2]

for i in range(nombrePlan):
    A=AttributionVitesse(0.1)
    TVitx[0,i]=A[0]
    TVity[0,i]=A[1]
    TVitz[0,i]=A[2]
    





ttab=np.arange(dt,temps,dt)#ATTENTION AU TEMPS MODULAIRE
print("Début des calculs")
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
for i in range(nombrePlan):
    ax.plot(TPosx[:,i],TPosy[:,i],TPosz[:,i],label=str(i))
plt.legend()
plt.show()
