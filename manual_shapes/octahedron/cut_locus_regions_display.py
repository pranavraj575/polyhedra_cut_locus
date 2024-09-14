from matplotlib import pyplot as plt
import numpy as np
from src.utils import rowtation
import os

plt.rcParams["figure.figsize"] = (2*2*np.sqrt(3), 2*3)
plt.rcParams['font.size'] = 18

vertices = np.array(
    (rowtation(np.pi/6).flatten(),
     rowtation(5*np.pi/6).flatten(),
     rowtation(-np.pi/2).flatten(),)
)*2

print(vertices.shape)
idxs = [i%3 for i in range(4)]
plt.plot(vertices[idxs, 0], vertices[idxs, 1], '--',
         # label='border',
         color='black')

alpha = .5
plt.fill_between([0, np.sqrt(3)/2, np.sqrt(3)], [0, -1/2, 1], [0, 1/2, 1], label='$\\tau_0$', alpha=alpha)
plt.fill_between([0, np.sqrt(3)], 1, [0, 1], label='$\\tau_1$', alpha=alpha)
plt.fill_between([0, -np.sqrt(3)], 1, [0, 1], label='$\\tau_2$', alpha=alpha)
plt.fill_between([0, -np.sqrt(3)/2, -np.sqrt(3)], [0, -1/2, 1], [0, 1/2, 1], label='$\\tau_3$', alpha=alpha)
plt.fill_between([-np.sqrt(3)/2, 0], [-1/2, -2], [-1/2, 0], label='$\\tau_4$', alpha=alpha)
plt.fill_between([np.sqrt(3)/2, 0], [-1/2, -2], [-1/2, 0], label='$\\tau_5$', alpha=alpha)

plt.plot([0, 0], [-2, 1], '--',
         color='black')
plt.plot([-np.sqrt(3)/2, np.sqrt(3)], [-1/2, 1], '--',
         color='black')
plt.plot([np.sqrt(3)/2, -np.sqrt(3)], [-1/2, 1], '--',
         color='black')

plt.legend()
plt.xticks([])
plt.yticks([])

plt.savefig(os.path.join('..', '..', 'images', 'octahedron', 'cut_locus_isometric_regions.png'))
# plt.show()
