from shapes import *
from shape_creation import Antiprism
from display_utils import *
import os

# if __name__ == "__main__":
for n in range(3, 20):
    # n=7
    # for mult in (np.arange(10)+1)/10:
    print(n)
    cube = Antiprism(n)

    r = np.sqrt(3)/np.tan(np.pi/n)
    R = r/np.cos(np.pi/n)

    fn = 0

    p = np.array([[0.], [0.]])
    # p = coltation(-np.pi/2 - np.pi/n)*r*.1
    # p = np.array([[0.], [-r*mult]])
    # p = coltation(-np.pi/2 - 2*np.pi/n)*r*mult
    # p = coltation(-np.pi/2 - np.pi/n)*.05*r**2+coltation(-np.pi/2)*.05*np.sqrt(r)

    cube.add_point_to_face((p, {'color':'black', 's':40}), fn)

    folder = os.path.join('images', 'antiprisms', str(n))
    if not os.path.exists(folder):
        os.makedirs(folder)
    name = 'p_' + str(tuple(p.flatten())) + '_face_' + str(fn) + '.png'
    cube.plot_faces(save_image=os.path.join(folder, name), show=False, figsize=((n)*2.5, 8), legend=lambda i, j:i == 1 or i == 2,
                    voronoi=(p, fn, 5))  # DIAMETER ONLY TRUE FOR TOP OR BOTTOM FACES
