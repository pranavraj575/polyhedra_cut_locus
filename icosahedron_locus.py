from shapes import *
from shape_creation import Icosahedron
from display_utils import *
import os

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
    # p = np.array([[0.], [0.2]])
    # p=.1*coltation(-np.pi/6)
    # p=.3*coltation(-np.pi/6)
    # p = np.array([[0.], [0.2]])
    # p = np.array([[0.], [0.2]])+.125*coltation(-np.pi/6)
    # p = np.array([[0.], [0.2]])-.125*coltation(-np.pi/6)

    """
    now for each of these, we have two dimensions of freedom for q (vertex)

    thus, we can embed these as 4-hypercubes, and we need at least 5 sets
    """

    cube.add_point_to_face((p, {'color':'black', 's':40}), fn)

    radii = [1.4 + np.sqrt(i)/1.4 for i in range(110)]
    i = 0
    color_dic = {r:color for (r, color) in zip(radii, rb_gradient(len(radii)))}
    for r, color in zip(radii, rb_gradient(len(radii))):
        i += 1
        A = Arc(p, 0, np.pi*2, r)
        plot = False
        # cube.add_arc_end_to_face(A, fn, arc_info={"color":color, 'plot':plot})
    # cube.add_all_cut_locus_points(point_info={'color':'black','s':2},conditional_point_info=lambda r:{'color':color_dic[r]})

    folder = os.path.join('images', 'icosahedron')
    if not os.path.exists(folder):
        os.makedirs(folder)
    name = 'p_' + str(tuple(p.flatten())) + '_face_' + str(fn) + '.png'
    cube.plot_faces(save_image=os.path.join(folder, name), show=False, figsize=(15, 12), voronoi=(p, fn, None))
