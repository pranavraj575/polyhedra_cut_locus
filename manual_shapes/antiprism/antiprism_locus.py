from src.shapes import *
from src.shape_creation import Antiprism
import os

for n in range(3, 10):
    print(n)
    c = 5
    r = np.sqrt(3)/np.tan(np.pi/n)
    R = r/np.cos(np.pi/n)
    fn = 0
    folder = os.path.join('../../images', 'antiprisms', str(n))
    if not os.path.exists(folder):
        os.makedirs(folder)

    for mult in (np.arange(c))*.5/(c - 1):
        p1 = coltation(-np.pi/2 - np.pi/n)*R*mult
        cube = Antiprism(n)
        cube.add_point_to_face(p1, fn, {'color': 'purple', 's': 40})
        name = 'p_' + str(tuple(p1.flatten())) + '_face_' + str(fn) + '.png'
        cube.plot_faces(save_image=os.path.join(folder, name), show=False, figsize=((n)*2.5, 8), legend=lambda i, j: i == 1 or i == 2,
                        voronoi=(p1, fn, 5))

        for mult2 in (np.arange(c))*.25/(c - 1):
            cube = Antiprism(n)
            p2 = coltation(-np.pi/2 - np.pi/n)*R*mult + coltation(-np.pi/2)*mult2*R

            cube.add_point_to_face(p2, fn, {'color': 'purple', 's': 40})
            name = 'p_' + str(tuple(p2.flatten())) + '_face_' + str(fn) + '.png'
            cube.plot_faces(save_image=os.path.join(folder, name), show=False, figsize=((n)*2.5, 8), legend=lambda i, j: i == 1 or i == 2,
                            voronoi=(p2, fn, 5))
