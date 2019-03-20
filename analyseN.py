import matplotlib.pyplot as plt
import numpy as np


N=10
l=6
pression=[0.04891216252280707,0.13031318849004142,0.011151842608221524]
energie=[1.1567332343198415,2.8215423771161885, 0.4391645209927484]
#b=  1.5338358181213643
#a=  4.629290983057639


N=20
l=6
pression=[0.19033162670297407,0.40765062979185845,0.0353772162114431]
energie=[2.0966208487182825,4.24775858051845, 0.7288682005676108]
#b=  1.318750906391287
#a=  4.380430676745182

N=15
l=8
pression=[0.13866616702728904,0.07111193989578121,0.018454452538211045]
energie=[4.652862524537208,2.399793347896256,0.7088505667250177]
#b=  1.299364662157295
#a=  3.166670559277799

N=25
l=8
pression=[0.010988083393956209,0.19406474091886247,0.5255060874673629]
energie=[0.560359125718088,3.788333118384508,10.379221159019384]
#b=  1.2997597758572788
#a=  5.216339506252895


v=l**3

coeff=np.polyfit(pression, energie, 1)
print(coeff)
b=v/N-coeff[0]
a=coeff[1]/(N/v-b*(N/v)**2)

print('b= ',b)
print('a= ',a)
    

plt.plot(pression,energie,"o")
plt.show()