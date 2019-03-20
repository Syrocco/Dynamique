import matplotlib.pyplot as plt
import numpy as np
"""
pression= 0.529343220585 temperature= 1.04680483866
pression= 0.517005334184 temperature= 1.03307557194
"""
"""
N=20
l=6
pression=[0.6319088998643482,0.10156223721809352,0.25480608720225734,0.9767347951538655]
energie=[6.334141264750083,0.9702756246836863,2.633886564515918,9.573286505731403]
#b=  1.0044675978229574
#a=  0.7836634476555281
"""
"""
N=20
l=4
pression=[0.0485405046442334,0.3960252845357385,0.18545220464826132,0.7494009238298456]
energie=[0.21233952444206908,1.0418656344196735,0.5241644729997249,1.7542557484372823]
#b=  0.9941112410177757
#a=  0.5684051729365055
"""


"""
N=10
l=6
pression=[0.012178885604389829, 0.1309613237823359,0.5597852738550967,0.33117999031307194]
energie=[0.2300719605734232,2.8714060060827513,11.899855013691347,6.6070777184672504]
#b=  0.5474359189450055
#a=  -0.8984399606753093
"""

"""
N=10
l=4
pression=[2.054909704308492, 0.20144410999369294,0.5073851104147705,0.04516078889076431]
energie=[11.739466805324565,1.2205338205844152,2.8809016552129223,0.31652760176729317]
#b=  0.7145262766265112
#a=  0.3377106003705477
"""
"""
N=100
l=10
pression=[0.0773544805432702, 0.29787518765179355,0.8925834606972957]
energie=[ 0.7551581569616538,2.68184020865967,8.04799757730767]
#b=  1.0380960679176905
#a=  0.457404583129348
"""
"""
N=50
l=10
pression=[0.4413456842810395, 0.14993922655454617,0.2767044756740931]
energie=[8.379892367841835,2.91452661895624,5.2593156567441]
#b=  1.235198987793165
#a=  1.8906783667816462
"""

"""
N=5
l=3
pression=[0.09292977134282138,0.15947874423655817,0.25996489821706076]
energie=[0.479567299661353,0.8124898818268574, 1.2200355390596562]
#b=  0.9973152669400873
#a=  0.5657875907806607
"""


N=10
l=10
pression=[0.012020345964736652,0.02919821351670362,0.05299914988979195,0.00011989238049744969]
energie=[1.204142305029632,2.8670973197640133, 5.226671526862875,0.0867634417706489]
#b=  2.6750114251314585
#a=  5.221387223775845


v=l**3

coeff=np.polyfit(pression, energie, 1)
print(coeff)
b=v/N-coeff[0]
a=coeff[1]/(N/v-b*(N/v)**2)

print('b= ',b)
print('a= ',a)
    

plt.plot(pression,energie,"o")
plt.show()
