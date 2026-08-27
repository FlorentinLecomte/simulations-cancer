import numpy as np
import matplotlib.pyplot as plt

data=np.fromfile('data.bin', dtype=np.float64).reshape(-1,4)

times=np.unique(data[:,0])

for i, t in enumerate(times):
    time=data[data[:,0]==t]
    x=time[:,1]
    y=time[:,2]
    b=time[:,3]

    fig,ax=plt.subplots(figsize=(5,5))
    ax.scatter(x, y, c=b, cmap='grey', marker='h', s=30, edgecolors='black', linewidths=0.4, vmin=0, vmax=1)
    ax.set_title(f't={t:.2f}')
    ax.axis('off')

    plt.savefig(f'figures/cancer_{t:.1f}.png')
    plt.close(fig)
