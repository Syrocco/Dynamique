import numpy as np
import matplotlib.pyplot as plt
import numpy.random as rd
from mpl_toolkits.mplot3d import Axes3D




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
                                  
def AttributionInitiale(rayon,Vitesse,):
	particule=GenPosition(nombrePlan, rayon, "Cube")
	for i in range(nombrePlan):
  	  TPosx[0,i]=particule[i,0]
  	  TPosy[0,i]=particule[i,1]
  	  TPosz[0,i]=particule[i,2]

	for i in range(nombrePlan):
		A=GenVitesse(Vitesse)
		TVitx[0,i]=A[0]
		TVity[0,i]=A[1]
		TVitz[0,i]=A[2]


def CalculAcceleration():
	d=distance(TPosx[i,planete],TPosx[i,planeteAutre],TPosy[i,planete],TPosy[i,planeteAutre],TPosz[i,planete],TPosz[i,planeteAutre])
	a=((12/d**13)-(6/d**7))/1
	return a*(TPosx[i,planete]-TPosx[i,planeteAutre])/d, a*(TPosy[i,planete]-TPosy[i,planeteAutre])/d,a*(TPosz[i,planete]-TPosz[i,planeteAutre])/d, d 

def CalculVitesseEtPosition():
    TVitx[i,planete]=TVitx[i-1,planete]+dt*ax
    TVity[i,planete]=TVity[i-1,planete]+dt*ay	
    TVitz[i,planete]=TVitz[i-1,planete]+dt*az
    TPosx[i,planete]=TPosx[i-1,planete]+dt*TVitx[i,planete]
    TPosy[i,planete]=TPosy[i-1,planete]+dt*TVity[i,planete]
    TPosz[i,planete]=TPosz[i-1,planete]+dt*TVitz[i,planete] 

def CalculPosition():
    TPosx[i+1,planete]=2*TPosx[i,planete]-TPosx[i-1,planete]+ax*dt**2
    TPosy[i+1,planete]=2*TPosy[i,planete]-TPosy[i-1,planete]+ay*dt**2
    TPosz[i+1,planete]=2*TPosz[i,planete]-TPosz[i-1,planete]+az*dt**2
    TVitx[i,planete]=(TPosx[i+1,planete]-TPosx[i-1,planete])/(2*dt)
    TVity[i,planete]=(TPosy[i+1,planete]-TPosy[i-1,planete])/(2*dt)
    TVitz[i,planete]=(TPosz[i+1,planete]-TPosz[i-1,planete])/(2*dt)
    
                
    
    
def PosVit_t1():
    for planete in range(nombrePlan):
        ax=0
        ay=0
        az=0
        for planeteAutre in range(nombrePlan):
            if planeteAutre!=planete:
                d=distance(TPosx[0,planete],TPosx[0,planeteAutre],TPosy[0,planete],TPosy[0,planeteAutre],TPosz[0,planete],TPosz[0,planeteAutre])
                a=((12/d**13)-(6/d**7))
                ax+=a*(TPosx[0,planete]-TPosx[0,planeteAutre])/d               
                ay+=a*(TPosy[0,planete]-TPosy[0,planeteAutre])/d
                az+=a*(TPosz[0,planete]-TPosz[0,planeteAutre])/d


        TVitx[1,planete]=TVitx[0,planete]+dt*ax
        TVity[1,planete]=TVity[0,planete]+dt*ay
        TVitz[1,planete]=TVitz[0,planete]+dt*ax
        TPosx[1,planete]=TPosx[0,planete]+dt*TVitx[1,planete]
        TPosy[1,planete]=TPosy[0,planete]+dt*TVity[1,planete]
        TPosz[1,planete]=TPosz[0,planete]+dt*TVitz[1,planete]  
                

def distance(x1,x2,y1,y2,z1,z2):
	return np.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)

def Ecinetique(vx,vy,vz):
	return 0.5*1*(vx**2+vy**2+vz**2)

def Epotentielle(d):
	return 0.5*((1/d**12)-(1/d**6))   #Multiplication par 0.5 parceque je compte deux fois l'energie: pour un couple {i,j} de particules, je compte Eij et Eji, donc il faut diviser par deux

##############################################################################################################################################################################################

rd.seed(4)


nombrePlan=10

temps=25
dt=0.001
N=int(temps/dt)

T=np.zeros(N)

TPosx=np.zeros((N+1,nombrePlan))
TPosy=np.zeros((N+1,nombrePlan))
TPosz=np.zeros((N+1,nombrePlan))
TVitx=np.zeros((N+1,nombrePlan))
TVity=np.zeros((N+1,nombrePlan))
TVitz=np.zeros((N+1,nombrePlan))



Epot=np.zeros(N)
Ecin=np.zeros(N)
ttab=np.arange(2*dt,temps+2*dt,dt)#ATTENTION AU TEMPS MODULAIRE


AttributionInitiale(2,0.1)

  
PosVit_t1()


print("Début des calculs")
for i in range(1,N):
    #print(i,"/",N)
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
                Epot[i]=Epot[i]+Epotentielle(a[3])
        CalculPosition()
        Ecin[i-1]=Ecin[i-1]+Ecinetique(TVitx[i-2,planete],TVity[i-2,planete],TVitz[i-2,planete])  
        
       
        
        
Etot=Epot[3:-1]+Ecin[3:-1]   
fig = plt.figure()
axes = fig.gca(projection='3d')
for i in range(nombrePlan):
    axes.plot(TPosx[:,i],TPosy[:,i],TPosz[:,i],label=str(i))
plt.figure()
plt.plot(ttab[2:-2],Epot[3:-1],label="Epot")
plt.plot(ttab[2:-2],Ecin[3:-1],label="Ecin")
plt.plot(ttab[2:-2],Etot,label="Etot")
plt.legend()
plt.show()
