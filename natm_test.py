import numpy as np
from mmsim import natm
from mmsim import natmslalib

warr=np.arange(3500, 13000, 500)

t=273.15+12
p=101325*0.77

#narr=natm(warr, t, p)
narr=natmslalib(warr, t, p)

for i in range(len(narr)):
    print(narr[i])

print("T ; P = ", t-273.15, p/101325)