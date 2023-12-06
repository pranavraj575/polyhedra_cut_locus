from shapes import *
from shape_creation import Octahedron
from display_utils import *
import os
"""
OCTAHEDRON IS AT LEAST 5!!
"""
if __name__ == "__main__":



    cube = Octahedron()

    fn = 1

    """
    two dimensions of freedom for (p,q)
    """
    p = np.array([[0.], [0.]])
    p = np.array([[0.], [0.1]])
    p=.1*coltation(-np.pi/6)
    p = np.array([[0.], [0.2]])
    p = np.array([[0.], [0.2]])+.075*coltation(-np.pi/6)
    p = np.array([[0.], [0.2]])-.075*coltation(-np.pi/6)

    """
    now for each of these, we have two dimensions of freedom for q (vertex)
    
    thus, we can embed these as 4-hypercubes, and we need at least 5 sets
    """

    cube.add_point_to_face((p, {'color':'black','s':40}), fn)

    radii = [1.4 + np.sqrt(i)/3 for i in range(140)]
    i=0
    color_dic={r:color for (r, color) in zip(radii, rb_gradient(len(radii)))}
    for r, color in zip(radii, rb_gradient(len(radii))):
        i+=1
        A = Arc(p, 0, np.pi*2, r)
        plot=False
        cube.add_arc_end_to_face(A, fn, arc_info={"color":color, 'plot':plot})
    cube.add_all_cut_locus_points(point_info={'color':'black','s':2},conditional_point_info=lambda r:{'color':color_dic[r]})

    folder=os.path.join('images','octahedron')
    if not os.path.exists(folder):
        os.makedirs(folder)
    name='p_'+str(tuple(p.flatten()))+'_face_'+str(fn)+'.png'
    cube.plot_faces(save_image=os.path.join(folder,name),show=False,figsize=(16,8))
    quit()
    c = 100

    radii = [1 + i/2 for i in range(6)]
    for r, color in zip(radii, rb_gradient(len(radii))):
        cc = int(r*c)
        for theta in [np.pi*2/cc*i for i in range(cc)]:
            v = r*np.array([[np.cos(theta)], [np.sin(theta)]])
            cube.add_path_end_to_face(p, v, fn, {'color':color, 's':4, 'plot':True})

    cube.plot_faces()

    for theta in (-np.pi/6,):
        rng = 3
        vbar = np.array([[np.cos(theta)], [np.sin(theta)]])
        for r in [i*rng/c for i in range(c)]:
            v = r*vbar
            cube.add_path_end_to_face(p, v, fn, {'color':'red', 's':4, 'plot':True})
