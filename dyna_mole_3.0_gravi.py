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
    NormEtVit=norme(vit_xyz)/rd.normal(T,500)   #constante permettant de normer vit_xyz pour ensuite lui donner une norme prise au hasard dans une gaussienne

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
                                  
                                  
def PosVit_t1(px,py,pz,vx,vy,vz):
    for planete in range(nombrePlan):
        ax=0
        ay=0
        az=0
        for planeteAutre in range(nombrePlan):
            if planeteAutre!=planete:
                d=distance(px[0,planete],px[0,planeteAutre],py[0,planete],py[0,planeteAutre],pz[0,planete],pz[0,planeteAutre])
                a=-G*Tmasse[planeteAutre]/d**2
                ax+=a*(px[0,planete]-px[0,planeteAutre])/d               
                ay+=a*(py[0,planete]-py[0,planeteAutre])/d
                az+=a*(pz[0,planete]-pz[0,planeteAutre])/d


        vx[1,planete]=vx[0,planete]+dt*ax
        vy[1,planete]=vy[0,planete]+dt*ay	
        vz[1,planete]=vz[0,planete]+dt*az
        px[1,planete]=px[0,planete]+dt*vx[1,planete]
        py[1,planete]=py[0,planete]+dt*vy[1,planete]
        pz[1,planete]=pz[0,planete]+dt*vz[1,planete]  
                                   


def distance(x1,x2,y1,y2,z1,z2):
	return np.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)

def Ecinetique(m,vx,vy,vz):
	return 0.5*m*(vx**2+vy**2+vz**2)

def Epotentielle(m1,m2,d):
	return -0.5*G*m1*m2/d   #Multiplication par 0.5 parceque je compte deux fois l'energie: pour un couple {i,j} de particule, je compte Eij et Eji, donc il faut diviser par deux






#Tmasse=[10**20,10**20,10**20,10**20,10**20]
Tmasse=[10**20,10**5,10**5,10**5,10**5]    
nombrePlan=len(Tmasse)

G = 6.67428e-11

temps=3
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



TPosx[0,0]=5
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

"""
particule=AttributionPosition(nombrePlan, 10000, "Cube")
for i in range(nombrePlan):
    TPosx[0,i]=particule[i,0]
    TPosy[0,i]=particule[i,1]
    TPosz[0,i]=particule[i,2]
 
    

for i in range(nombrePlan):
    A=AttributionVitesse(15)
    TVitx[0,i]=A[0]
    TVity[0,i]=A[1]
    TVitz[0,i]=A[2]"""





ttab=np.arange(dt,temps,dt)#ATTENTION AU TEMPS MODULAIRE

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
axes = fig.gca(projection='3d')
for i in range(len(Tmasse)):
    axes.plot(TPosx[:,i],TPosy[:,i],TPosz[:,i],label=str(i))
plt.legend()
plt.figure()
plt.plot(ttab,Epot,label='epot')
plt.plot(ttab,Ecin,label='ecin')
plt.plot(ttab,Etot,label='etot')
plt.legend()
plt.show()


