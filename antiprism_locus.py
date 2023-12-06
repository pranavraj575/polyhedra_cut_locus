from shapes import *
from shape_creation import Antiprism
from display_utils import *
import os

#if __name__ == "__main__":
for n in range(4,17):
    print(n)
    #n=8
    cube = Antiprism(n)

    r = np.sqrt(3)/np.tan(np.pi/n)
    R=r/np.cos(np.pi/n)

    fn = 0

    p = np.array([[0.], [0.]])


    cube.add_point_to_face((p, {'color':'black','s':40}), fn)

    c=100

    radii = [r + (R+3)*i/c for i in range(c)]
    i=0
    color_dic={r:color for (r, color) in zip(radii, rb_gradient(len(radii)))}
    for r, color in zip(radii, rb_gradient(len(radii))):
        i+=1
        A = Arc(p, 0, np.pi*2, r)
        plot=False
        cube.add_arc_end_to_face(A, fn, arc_info={"color":color, 'plot':plot})
    cube.add_all_cut_locus_points(point_info={'color':'black','s':2},conditional_point_info=lambda r:{'color':color_dic[r]})

    folder=os.path.join('images','antiprisms',str(n))
    if not os.path.exists(folder):
        os.makedirs(folder)
    name='p_'+str(tuple(p.flatten()))+'_face_'+str(fn)+'.png'
    cube.plot_faces(save_image=os.path.join(folder,name),show=False,figsize=((n)*4,8),legend=lambda i,j:i==1 or i==2)