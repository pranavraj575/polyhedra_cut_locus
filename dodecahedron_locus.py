from shapes import *
from shape_creation import Dodecahedron
from display_utils import *
import os
"""
DODECAHEDRON IS AT LEAST 5!!
"""
if __name__ == "__main__":


    cube = Dodecahedron()

    fn = 0
    """
    two dimensions of freedom for (p,q)
    """
    p = np.array([[0.], [0.]])
    p = np.array([[0.], [0.2]])
    p = .2*coltation(np.pi/2-np.pi*2/5)
    #p = .2*coltation(np.pi/2-np.pi*2/5)+.075*coltation(np.pi/2)
    #p = .2*coltation(np.pi/2-np.pi*2/5)+.075*coltation(np.pi/2)


    """
    now for each of these, we have two dimensions of freedom for q (vertex)

    thus, we can embed these as 4-hypercubes, and we need at least 5 sets
    """

    cube.add_point_to_face((p, {'color':'black','s':40}), fn)

    radii = [1.4 + np.sqrt(i)/2 for i in range(80)]
    i=0
    color_dic={r:color for (r, color) in zip(radii, rb_gradient(len(radii)))}
    for r, color in zip(radii, rb_gradient(len(radii))):
        i+=1
        A = Arc(p, 0, np.pi*2, r)
        plot=False
        cube.add_arc_end_to_face(A, fn, arc_info={"color":color, 'plot':plot})
    cube.add_all_cut_locus_points(point_info={'color':'black','s':2},conditional_point_info=lambda r:{'color':color_dic[r]})

    folder=os.path.join('images','dodecahedron')
    if not os.path.exists(folder):
        os.makedirs(folder)
    name='p_'+str(tuple(p.flatten()))+'_face_'+str(fn)+'.png'
    cube.plot_faces(save_image=os.path.join(folder,name),show=False,figsize=(15,9))
