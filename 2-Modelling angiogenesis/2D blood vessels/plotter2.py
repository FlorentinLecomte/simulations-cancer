import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib.patches as mpatches

data = np.fromfile('continuous_data.bin', dtype=np.float64)
df = data.reshape(-1, 5, 402, 402)

couleurs_types = ['lightgray', 'darkgray', 'dimgray', 'black', (0.0, 0.35, 0.2, 0.4)]
cmap_personnalisee = ListedColormap(couleurs_types)
cmap_vaisseaux = ListedColormap([(0, 0, 0, 0), 'crimson'])

cancer_cells = np.loadtxt('discrete_data.txt')
cancer_cells[:, 0].astype(np.int64)
t_steps = np.unique(cancer_cells[:, 0]).astype(np.int64)

for i, t_step in enumerate(t_steps):
    mask = (cancer_cells[:, 0] == t_step)
    x = cancer_cells[mask, 1]
    y = cancer_cells[mask, 2]
    types = cancer_cells[mask, 3]

    fig, ax = plt.subplots(figsize=(6, 6))

    v_matrix = df[i, 4, :, :]
    ax.imshow(v_matrix, extent=[0, 1, 0, 1], origin='lower', cmap=cmap_vaisseaux, vmin=0, vmax=1, zorder=1)

    patch_vaisseaux = mpatches.Patch(color='crimson', label='Blood vessels')
    handles = [patch_vaisseaux]
    labels_modifies = ['Blood vessels']

    if len(types) > 0:
        ordre_dessin = np.argsort(types != 5)

        x_trie = x[ordre_dessin]
        y_trie = y[ordre_dessin]
        types_trie = types[ordre_dessin]

        scatter = ax.scatter(x_trie/400, y_trie/400, s=1, c=types_trie, cmap=cmap_personnalisee, vmin=1, vmax=5)
        scat_handles, scat_labels = scatter.legend_elements()
        handles=list(scat_handles)
        labels_modifies = []

        for label in scat_labels:
            if '5' in label: labels_modifies.append('Dead cell')
            elif '1' in label: labels_modifies.append('Type 1 cell')
            elif '2' in label: labels_modifies.append('Type 2 cell')
            elif '3' in label: labels_modifies.append('Type 3 cell')
            elif '4' in label: labels_modifies.append('Type 4 cell')
            else: labels_modifies.append(label)

        handles.insert(0, patch_vaisseaux)
        labels_modifies.insert(0, 'Blood vessels')

    ax.legend(handles=handles, labels=labels_modifies, title="Types of cells", loc="upper right")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_title(f"Cancer cells - t={t_step/2000:3f}")
    ax.set_xlabel("x")
    ax.set_ylabel("y")

    filename = f"figures/frame_{i:03d}.png"
    plt.savefig(filename)
    plt.close(fig)


saving_step=2000

for t_idx in range(df.shape[0]):
    grid = df[t_idx, 2, 1:-1, 1:-1]

    fig, ax = plt.subplots()
    im = ax.imshow(grid, cmap='Greys', origin='lower', extent=[0, 1, 0, 1])
    ax.set_title(f"oxygen - t={t_idx*saving_step/2000:3f}")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    fig.colorbar(im, ax=ax, label="Concentration d'oxygène")

    filename = f"figures/frame2_{t_idx:03d}.png"
    plt.savefig(filename)
    plt.close(fig)

import imageio.v2 as imageio
import glob

images = sorted(glob.glob("figures/frame_*.png"))
frames = [imageio.imread(f) for f in images]

imageio.mimsave('animations/Cancer.mp4', frames, fps=5)
images = sorted(glob.glob("figures/frame2_*.png"))
frames = [imageio.imread(f) for f in images]

imageio.mimsave('animations/Oxygen.mp4', frames, fps=5)

for t_idx in range(df.shape[0]):
    grid = df[t_idx, 0, 1:-1, 1:-1]

    fig, ax = plt.subplots()
    im = ax.imshow(grid, cmap='Greys', origin='lower', extent=[0, 1, 0, 1])
    ax.set_title(f"ECM - t={t_idx*saving_step/2000:3f}")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    fig.colorbar(im, ax=ax, label="ECM concentration")

    filename = f"figures/frame3_{t_idx:03d}.png"
    plt.savefig(filename)
    plt.close(fig)

images = sorted(glob.glob("figures/frame3_*.png"))
frames = [imageio.imread(f) for f in images]

imageio.mimsave('animations/ECM.mp4', frames, fps=5)

for t_idx in range(df.shape[0]):
    grid = df[t_idx, 1, 1:-1, 1:-1]

    fig, ax = plt.subplots()
    im = ax.imshow(grid, cmap='Greys', origin='lower', extent=[0, 1, 0, 1])
    ax.set_title(f"MDE - t={t_idx*saving_step/2000:3f}")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    fig.colorbar(im, ax=ax, label="MDE concentration")

    filename = f"figures/frame4_{t_idx:03d}.png"
    plt.savefig(filename)
    plt.close(fig)

images = sorted(glob.glob("figures/frame4_*.png"))
frames = [imageio.imread(f) for f in images]

imageio.mimsave('animations/MDE.mp4', frames, fps=5)

for t_idx in range(df.shape[0]):
    grid = df[t_idx, 3, 1:-1, 1:-1]

    fig, ax = plt.subplots()
    im = ax.imshow(grid, cmap='Greys', origin='lower', extent=[0, 1, 0, 1])
    ax.set_title(f"TAF - t={t_idx*saving_step/2000:3f}")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    fig.colorbar(im, ax=ax, label="TAF concentration")

    filename = f"figures/frame5_{t_idx:03d}.png"
    plt.savefig(filename)
    plt.close(fig)

images = sorted(glob.glob("figures/frame5_*.png"))
frames = [imageio.imread(f) for f in images]

imageio.mimsave('animations/TAF.mp4', frames, fps=5)