import os
import numpy as np

from src.shape_creation import Tetrahedron

if __name__ == "__main__":

    cube = Tetrahedron()

    fn = 3
    p = np.array([[-.4], [-.6]])

    cube.add_point_to_face(p, fn, {'color': 'purple', 's': 40})

    folder = os.path.join('../../images', 'tetrahedron')
    if not os.path.exists(folder):
        os.makedirs(folder)
    name = 'p_' + str(tuple(p.flatten())) + '_face_' + str(fn) + '.png'
    cube.plot_faces(save_image=os.path.join(folder, name), show=False, figsize=(12, 8), voronoi=(p, fn, None))
