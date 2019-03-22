import numpy as np
import matplotlib.pyplot as plt
import numpy.random as rd
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection




################################################################
###-----------------Définition des fonctions-----------------###    
################################################################ 

def plot_cube(cube_definition):  #La fonction volée
    cube_definition_array = [
        np.array(list(item))
        for item in cube_definition
    ]

    points = []
    points += cube_definition_array
    vectors = [
        cube_definition_array[1] - cube_definition_array[0],
        cube_definition_array[2] - cube_definition_array[0],
        cube_definition_array[3] - cube_definition_array[0]
    ]

    points += [cube_definition_array[0] + vectors[0] + vectors[1]]
    points += [cube_definition_array[0] + vectors[0] + vectors[2]]
    points += [cube_definition_array[0] + vectors[1] + vectors[2]]
    points += [cube_definition_array[0] + vectors[0] + vectors[1] + vectors[2]]

    points = np.array(points)

    edges = [
        [points[0], points[3], points[5], points[1]],
        [points[1], points[5], points[7], points[4]],
        [points[4], points[2], points[6], points[7]],
        [points[2], points[6], points[3], points[0]],
        [points[0], points[2], points[4], points[1]],
        [points[3], points[6], points[7], points[5]]
    ]


    faces = Poly3DCollection(edges, linewidths=1, edgecolors='k')
    faces.set_facecolor((0,0,1,0.1))

    axes.add_collection3d(faces)




def norme(Vec3D):
    return np.sqrt(normeCarree(Vec3D))

def normeCarree(Vec3D):
    return Vec3D[0]**2+Vec3D[1]**2+Vec3D[2]**2

#Avec la vitesse/température T et l'écart type E, la fonction renvoie UNE vitesse prise dans une gaussienne
def GenVitesse(T,E):
    vit_xyz=np.zeros(3)
    for i in range(3):
        vit_xyz[i]=rd.random()
    NormEtVit=norme(vit_xyz)/rd.normal(T,E)   #constante permettant de normer vit_xyz pour ensuite lui donner une norme prise au hasard dans une gaussienne

    for i in range(3):
        signe=rd.random()
        if signe<0.5:
            vit_xyz[i]=vit_xyz[i]/NormEtVit
        else:
            vit_xyz[i]=-vit_xyz[i]/NormEtVit
    
    return np.array(vit_xyz)

#La fonction renvoie la distribution des corps selon la demi-longueur envoyée
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
        return np.array(T[:nombrePlanete])   #Je ne suis pas arrivé à faire un programme renvoyant tout le temps une matrice avec le nombre exacte de planetes, du coup, j'augmente un peu les
                                             #limites du cube et je ne prend que les "nombrePlanetes" premières valeurs de la distribution.


#La fonction attribue les positions et vitesses initiales aux corps                                  
def AttributionInitiale(rayon,Vitesse,Ecart):
    particule=GenPosition(nombrePlan, rayon, "Cube")
    for i in range(nombrePlan):
        TPosx[0,i]=particule[i,0]
        TPosy[0,i]=particule[i,1]
        TPosz[0,i]=particule[i,2]

    for i in range(nombrePlan):
        A=GenVitesse(Vitesse,Ecart)
        TVitx[0,i]=A[0]
        TVity[0,i]=A[1]
        TVitz[0,i]=A[2]
        


#Calcul de l'acceleration
def CalculAcceleration():
	d=distance(TPosx[i-1,planete],TPosx[i-1,planeteAutre],TPosy[i-1,planete],TPosy[i-1,planeteAutre],TPosz[i-1,planete],TPosz[i-1,planeteAutre])
	a=((12/d**13)-(6/d**7))/1
	return a*(TPosx[i-1,planete]-TPosx[i-1,planeteAutre])/d, a*(TPosy[i-1,planete]-TPosy[i-1,planeteAutre])/d,a*(TPosz[i-1,planete]-TPosz[i-1,planeteAutre])/d, d 

#Calcul de la position et de la vitesse avec la méthode d'euler semi-explicite
def CalculVitesseEtPosition():
    TVitx[i,planete]=TVitx[i-1,planete]+dt*ax
    TVity[i,planete]=TVity[i-1,planete]+dt*ay	
    TVitz[i,planete]=TVitz[i-1,planete]+dt*az
    TPosx[i,planete]=TPosx[i-1,planete]+dt*TVitx[i,planete]
    TPosy[i,planete]=TPosy[i-1,planete]+dt*TVity[i,planete]
    TPosz[i,planete]=TPosz[i-1,planete]+dt*TVitz[i,planete] 
    
    
#Teste si la particule à l'étape i se trouve dans la boite ou non                
def DansBoite(r):
        
    #r la demi longueur du cube
    if TPosx[i,planete]<-r:
        return "x-"
    if TPosx[i,planete]>r:
        return "x+"
    if TPosy[i,planete]<-r:
        return "y-"
    if TPosy[i,planete]>r:
        return "y+"
    if TPosz[i,planete]<-r:
        return "z-"
    if TPosz[i,planete]>r:
        return "z+"
    return("no")

#Modifie, selon le resultat de la fonction DansBoite(), la vitesse et la position de la particule à l'étape i pour simuler une collision, enregistre aussi la quantité de mouvement transmise aux parois
def modif(info):
    if info=="no":
        return()
    if info=='x-':
        TPosx[i,planete]=TPosx[i,planete] + 2*(-TailleBoite-TPosx[i,planete])
        Moment[i]=Moment[i]+abs(2*TVitx[i,planete])
        TVitx[i,planete]=-TVitx[i,planete]
        return()
    if info=='x+':
        TPosx[i,planete]=TPosx[i,planete] + 2*(TailleBoite-TPosx[i,planete])
        Moment[i]=Moment[i]+abs(2*TVitx[i,planete])
        TVitx[i,planete]=-TVitx[i,planete]
        return()
    if info=='y-':
        TPosy[i,planete]=TPosy[i,planete] + 2*(-TailleBoite-TPosy[i,planete])
        Moment[i]=Moment[i]+abs(2*TVity[i,planete])
        TVity[i,planete]=-TVity[i,planete]
        return()
    if info=='y+':
        TPosy[i,planete]=TPosy[i,planete] + 2*(TailleBoite-TPosy[i,planete])
        Moment[i]=Moment[i]+abs(2*TVity[i,planete])
        TVity[i,planete]=-TVity[i,planete]
        return()
    if info=='z-':
        TPosz[i,planete]=TPosz[i,planete] + 2*(-TailleBoite-TPosz[i,planete])
        Moment[i]=Moment[i]+abs(2*TVitz[i,planete])
        TVitz[i,planete]=-TVitz[i,planete]
        return()
    if info=='z+':
        TPosz[i,planete]=TPosz[i,planete] + 2*(TailleBoite-TPosz[i,planete])
        Moment[i]=Moment[i]+abs(2*TVitz[i,planete])
        TVitz[i,planete]=-TVitz[i,planete]
        

def distance(x1,x2,y1,y2,z1,z2):
	return np.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)

def Ecinetique(vx,vy,vz):
	return 0.5*1*(vx**2+vy**2+vz**2)

def Epotentielle(d):
	return 0.5*((1/d**12)-(1/d**6))   #Multiplication par 0.5 parceque je compte deux fois l'energie: pour un couple {i,j} de particule, je compte Eij et Eji, donc il faut diviser par deux



#####################################################################
###-----------------Initialisation des paramètres-----------------###    
#####################################################################    


#Nombre de corps
nombrePlan=30

#Durée de la simulation
temps=30

#Intervalle de temps (Il vaut mieux garder un multiple de 10 sinon, le ttab peut avoir des problemes de dimensionnement (N+1 colonnes plutot que N))
dt=0.01

#Nombre de simulation(s)
N=int(temps/dt)


#Définition des tableaux contenant les différentes données
ttab=np.arange(dt,temps,dt)

#Position et vitesse
TPosx=np.zeros((N,nombrePlan))
TPosy=np.zeros((N,nombrePlan))
TPosz=np.zeros((N,nombrePlan))
TVitx=np.zeros((N,nombrePlan))
TVity=np.zeros((N,nombrePlan))
TVitz=np.zeros((N,nombrePlan))

#Energie et quantité de mouvement
Epot=np.zeros(N)
Ecin=np.zeros(N)
Moment=np.zeros(N)  #Quantité de mouvement transmise aux parois, en valeur absolue puisqu'elle ne sert qu'à trouver une pression

#Demi longueur du cube dans lequel on place les corps et leurs vitesses + ecart type de la gaussienne en t=0
TailleInitiale=1.7
VitesseInitiale=3
EcartType=0.00

#Taille de la boite dans laquelle se passe les collisions (à garder STRICTEMENT inférieur à: TailleInitiale)
TailleBoite=4

cube_definition = [
    (-TailleBoite,-TailleBoite,-TailleBoite), (TailleBoite,-TailleBoite,-TailleBoite), (-TailleBoite,TailleBoite,-TailleBoite), (-TailleBoite,-TailleBoite,TailleBoite)
]

#Generation des conditions initiales
AttributionInitiale(TailleInitiale,VitesseInitiale,EcartType)


#temperature= 0.159402512101

############################################################
###-----------------Programme Principale-----------------###    
############################################################  


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
                Epot[i-1]=Epot[i-1]+Epotentielle(a[3])
        CalculVitesseEtPosition()
        modif(DansBoite(TailleBoite))    #A mettre en commentaire pour desactiver les collisions
        Ecin[i]=Ecin[i]+Ecinetique(TVitx[i,planete],TVity[i,planete],TVitz[i,planete])        
        
       
###############################################
###-----------------Calculs-----------------###    
###############################################


        
#Calcul de l'énergie totale     
Etot=Epot[1:-1]+Ecin[1:-1]


#Calcul de la pression exercée sur les parois 
aireBoite=6*(2*TailleBoite)**2
pas=100
nombreDivision=int(N/pas)

ttab2=np.linspace(0,dt*N,nombreDivision)
TMoment=np.zeros(nombreDivision)
for i in range(nombreDivision):
    TMoment[i]=sum(Moment[i*pas:(i*pas)+pas])
Tpress=TMoment/(pas*dt*aireBoite)

#Calcul temperature
kb=1

eneCinMoyenne=0
for i in range(int(N/2),N):
    eneCinMoyenne+=2*Ecin[i]/N
temperature=eneCinMoyenne/(1.5*kb*nombrePlan)


pression=sum(Tpress[:int(len(Tpress)/2)])*2/len(Tpress)   #Pression exercée en moyenne sur les parois pendant le temps (dt*N)/2

print("pression=", pression,"temperature=",temperature )  

        
###########################################################
###-----------------Affichage Graphique-----------------###    
###########################################################   
     
fig = plt.figure()
axes = fig.gca(projection='3d')
for i in range(nombrePlan):
    axes.plot(TPosx[:,i],TPosy[:,i],TPosz[:,i])
    
plot_cube(cube_definition)
    
plt.legend()
plt.figure()
plt.plot(ttab[:-1],Epot[1:-1],label="Epot")
plt.plot(ttab[:-1],Ecin[1:-1],label="Ecin")
plt.plot(ttab[:-1],Etot,label="Etot")

plt.figure()
plt.plot(ttab2,Tpress,label='Pression')
plt.figure()
plt.plot(ttab[:-1],Moment[1:-1],label="Quantité de mouvement")

plt.legend()
plt.show()
