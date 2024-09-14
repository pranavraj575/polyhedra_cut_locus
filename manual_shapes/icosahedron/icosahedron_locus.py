import os
import numpy as np

from src.utils import coltation
from src.shape_creation import Icosahedron

"""
ICOSAHEDRON IS AT LEAST 5!!
"""
if __name__ == "__main__":

    cube = Icosahedron()

    fn = 1
    """
    two dimensions of freedom for (p,q)
    """
    p = np.array([[0.], [0.]])
    p = np.array([[0.], [0.2]])
    p = .1*coltation(-np.pi/6)
    p = .3*coltation(-np.pi/6)
    p = np.array([[0.], [0.2]])
    p = np.array([[0.], [0.2]]) + .125*coltation(-np.pi/6)
    p = np.array([[0.], [0.2]]) - .125*coltation(-np.pi/6)

    """
    now for each of these, we have two dimensions of freedom for q (vertex)

    thus, we can embed these as 4-hypercubes, and we need at least 5 sets
    """

    cube.add_point_to_face(p, fn, {'color': 'purple', 's': 40})

    folder = os.path.join('../../images', 'icosahedron')
    if not os.path.exists(folder):
        os.makedirs(folder)
    name = 'p_' + str(tuple(p.flatten())) + '_face_' + str(fn) + '.png'
    cube.plot_faces(save_image=os.path.join(folder, name), show=False, figsize=(15, 12), voronoi=(p, fn, None))
