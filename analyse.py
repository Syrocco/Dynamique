import matplotlib.pyplot as plt
import numpy as np


N=20
l=6
pression=[0.112398949165,0.0399243526013,0.0121793541182,0.0656867985354 ]
energie=[1.08089947941,0.43678112101, 0.205587024662,0.680256852826]
#b=  2.01905368757
#a=  1.26936350515


N=10
l=6
pression=[0.0323864143119,0.0130739135327 ,0.0630879942723  ]
energie=[0.710231060915, 0.299286552782,1.29844195507]
b=  1.69749898878
a=  1.15317636996

"""
N=30
l=8
pression=[0.0558863356351,0.102735683851  ,0.158534874872  ]
energie=[0.903691612285, 1.64688191676,2.5688057558]
"""

v=l**3

coeff=np.polyfit(pression, energie, 1)
print(coeff)
b=v/N-coeff[0]
a=coeff[1]/(N/v-b*(N/v)**2)

print('b= ',b)
print('a= ',a)
    

plt.plot(pression,energie,"o")
plt.show()
