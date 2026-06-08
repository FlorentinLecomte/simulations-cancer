import numpy as np
import matplotlib.pyplot as plt

data=np.fromfile('cadherin.bin', dtype=np.float64).reshape(-1,5)

t=data[:,0]
Eb=data[:,1]
C=data[:,2]
Ec=data[:,3]
b=data[:,4]

plt.figure(figsize=(8,6))
plt.plot(t, Eb, '-', color='black', label=r'[E/$\beta$]')
plt.plot(t, C, '-.', color='black', label='[C]')
plt.plot(t, Ec, '--', color='black', label='[Ec]')
plt.plot(t, b, ':', color='black', label=r'[$\beta$]')

plt.xlabel('Time')
plt.ylabel('Concentrations')
plt.legend(loc='upper right')

plt.savefig(f'figures/EDO.png')
plt.close()