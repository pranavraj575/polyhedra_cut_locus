import numpy as np

from shapes import *

class Cube(Shape):
    # makes cube where faces 0,1,2,3 are the 'sides' in order (1 is to the right of 0)
    # 4 is the top, aligned with 0 on the right, 1 on the top
    # 5 is the bottom, aligned with 0 on the right, 1 on the bottom
    def __init__(self):
        super().__init__()
        for i in range(6):
            self.add_face()

        for i in range(4):
            curr: Face = self.faces[i]
            neigh = self.faces[(i + 1)%4]
            I = np.identity(2)
            s = np.array([[1], [0]])
            curr.add_boundary_paired(neigh, np.array([[1, 0]]), 1, -s, I, -s)

        top = self.faces[4]
        bot = self.faces[5]

        for i in range(4):
            curr = self.faces[i]
            angle = i*np.pi/2
            # "angle" of boundary of top face.
            # 0th face is the right one, so angle is 0
            s = np.array([[0], [-1]])
            theta = np.pi + (angle - np.pi/2)  # since we move top (pi/2) to the angle
            # adding pi to the angle because portal logic
            # 'exiting' through top, then rotating to the left leads to 'exiting' through left,
            # which is same as 'entering' on right
            T = rotation_T(theta)
            si = np.array([[np.cos(angle)], [np.sin(angle)]])
            # add_boundary_pair(curr, top, np.array([[0, 1]]), 1, s, T, si)
            curr.add_boundary_paired(top, np.array([[0, 1]]), 1, s, T, si)

        for i in range(4):
            curr = self.faces[i]
            angle = -i*np.pi/2
            # "angle" of boundary of bottom face.
            # 0th face is the right one, so angle is 0
            # inverted since we are going clockwise instead of counter
            s = np.array([[0], [1]])
            theta = np.pi + (angle - (-np.pi/2))  # since we move bottom (-pi/2) to the angle
            # adding pi to the angle because portal logic
            # 'exiting' through top, then rotating to the left leads to 'exiting' through left,
            # which is same as 'entering' on right
            T = rotation_T(theta)
            si = np.array([[np.cos(angle)], [np.sin(angle)]])
            # add_boundary_pair(curr, bot, np.array([[0, -1]]), 1, s, T, si)
            curr.add_boundary_paired(bot, np.array([[0, -1]]), 1, s, T, si)

    def faces_to_plot_n_m(self):
        def face_map(i, j):
            if i == 1:
                return self.faces[(j + 2)%4]
            if i == 0 and j == 1:
                return self.faces[4]
            if i == 2 and j == 1:
                return self.faces[5]
            return None

        return face_map, 3, 4


class Tetrahedron(Shape):
    # makes tetrahedron where face 0 borders face 1 on the right, 2 on the left and 3 on the bottom
    # 3 is upside down triangle for display purposes
    # each face has a circumcenter of radius 1
    def __init__(self):
        super().__init__()

        for i in range(4):
            self.add_face()

        for i in range(3):
            curr = self.faces[i]
            neigh = self.faces[(i + 1)%3]
            curr.add_boundary_paired(neigh, rowtation(np.pi/6), 1, -coltation(np.pi/6), rotation_T(-np.pi/3), coltation(np.pi*5/6))

        front, right, left, bottom = [self.faces[i] for i in range(4)]

        front.add_boundary_paired(bottom, np.array([[0, -1]]), 1, np.array([[0], [1]]), np.identity(2), np.array([[0], [1]]))

        right.add_boundary_paired(bottom, np.array([[0, -1]]), 1, np.array([[0], [1]]), rotation_T(-np.pi*2/3), coltation(-np.pi/6))

        left.add_boundary_paired(bottom, np.array([[0, -1]]), 1, np.array([[0], [1]]), rotation_T(np.pi*2/3), coltation(-np.pi*5/6))

    def faces_to_plot_n_m(self):
        def face_map(i, j):
            if i == 1 and j == 1:
                return self.faces[3]
            if i == 0:
                return self.faces[(j + 2)%3]
            return None

        return face_map, 2, 3


class Octahedron(Shape):
    # makes octahedron where the top half are faces 0,1,2,3 (1 is to the right of 0)
    # bottom half is 4,5,6,7 (4 is below 0, 5 below 1)
    # each face has a circumcenter of radius 1
    # we will make 0,1,2,3 normal triangles and 4,5,6,7 upside down triangles for visualization purposes
    def __init__(self):
        super().__init__()

        for i in range(8):
            self.add_face()

        for i in range(4):
            curr = self.faces[i]
            neigh = self.faces[(i + 1)%4]
            curr.add_boundary_paired(neigh, rowtation(np.pi/6), 1, -coltation(np.pi/6), rotation_T(-np.pi/3), coltation(np.pi*5/6))

        for i in range(4):
            curr = self.faces[4 + i]
            neigh = self.faces[4 + (i + 1)%4]
            curr.add_boundary_paired(neigh, rowtation(-np.pi/6), 1, -coltation(-np.pi/6), rotation_T(np.pi/3), coltation(-np.pi*5/6))

        for i in range(4):
            top = self.faces[i]
            bot = self.faces[4 + i]
            top.add_boundary_paired(bot, np.array([[0, -1]]), 1, np.array([[0], [1]]), np.identity(2), np.array([[0], [1]]))

    def faces_to_plot_n_m(self):
        def face_map(i, j):
            return self.faces[i*4 + j]

        return face_map, 2, 4


class Icosahedron(Shape):
    # makes icosahedron
    # top 5 are faces 0-4 (1 to the right of 0)
    # bottom 5 are faces 15-19 (16 to the right of 15, and all upside down)
    # middle 10 are faces 5-14 (11 to the right of 10, and alternating upside down)
    # (0,1,2,3,4) is above (5,7,9,11,13). (also 5,7,9,11,13 are all down)
    # (6,8,10,12,14) is above (15,16,17,18,19)
    def __init__(self):
        super().__init__()

        for i in range(20):
            self.add_face()

        for i in range(5):
            curr = self.faces[i]
            neigh = self.faces[(i + 1)%5]
            curr.add_boundary_paired(neigh, rowtation(np.pi/6), 1, -coltation(np.pi/6), rotation_T(-np.pi/3), coltation(np.pi*5/6))

        for i in range(5):
            curr = self.faces[i + 15]
            neigh = self.faces[15 + (i + 1)%5]
            curr.add_boundary_paired(neigh, rowtation(-np.pi/6), 1, -coltation(-np.pi/6), rotation_T(np.pi/3), coltation(-np.pi*5/6))

        for i in range(5):
            down = self.faces[5 + 2*i]
            up = self.faces[6 + 2*i]
            next_down = self.faces[5 + 2*((i + 1)%5)]

            down.add_boundary_paired(up, rowtation(-np.pi/6), 1, -coltation(-np.pi/6), np.identity(2), coltation(np.pi*5/6))
            up.add_boundary_paired(next_down, rowtation(np.pi/6), 1, -coltation(np.pi/6), np.identity(2), coltation(-np.pi*5/6))

        for i in range(5):
            top = self.faces[6 + 2*i]
            bottom = self.faces[15 + i]

            top.add_boundary_paired(bottom, rowtation(-np.pi/2), 1, coltation(np.pi/2), np.identity(2), coltation(np.pi/2))

        for i in range(5):
            top = self.faces[i]
            bottom = self.faces[5 + 2*i]

            top.add_boundary_paired(bottom, rowtation(-np.pi/2), 1, coltation(np.pi/2), np.identity(2), coltation(np.pi/2))

    def faces_to_plot_n_m(self):
        def face_map(i, j):
            if i == 0:
                return self.faces[j]
            if i == 3:
                return self.faces[j + 15]
            if i == 1:
                return self.faces[5 + 2*j]
            if i == 2:
                return self.faces[6 + 2*j]
            return None
        return face_map, 4, 5

class Dodecahedron(Shape):
    # makes dodecahedron
    # 0-2 on the top (1 to the right of 0)
    # 10-12 on bottom

    def __init__(self):
        super().__init__()
        for i in range(12):
            self.add_face()
        tau=2*np.pi
        for i in range(3):
            curr = self.faces[i]
            neigh = self.faces[(i + 1)%3]
            curr.add_boundary_paired(neigh, rowtation(tau/4-tau/5), 1, -coltation(tau/4-tau/5), rotation_T(tau/2-tau/5), coltation(tau/2+tau/5))

        for i in range(3):
            curr = self.faces[10+i]
            neigh = self.faces[10+(i + 1)%3]
        raise NotImplementedError

class NTorus(Shape):
    # makes n-torus
    def __init__(self,n):
        super().__init__()
        self.add_face()
        face=self.faces[0]
        face:Face
        for i in range(n):
            ei=np.zeros((n,1))
            ei[i,0]=1
            face.add_boundary_paired(face,ei.T,1,-ei,np.identity(n),-ei)


    def faces_to_plot_n_m(self):
        def face_map(i, j):
            return self.faces[0]
        return face_map, 1,1


class Torus(NTorus):
    # makes torus
    def __init__(self):
        super().__init__(2)

if __name__=="__main__":
    from display_utils import *
    cube = Tetrahedron()

    face:Face=cube.faces[0]
    fn = 0
    p = np.array([[0.5], [0.8]])

    cube.add_point_to_face((p, {'color':'red','s':20}), fn)

    radii = [1 + i/20 for i in range(90)]
    for r, color in zip(radii, rb_gradient(len(radii))):
        A = Arc(p, 0, np.pi*2, r)
        cube.add_arc_end_to_face(A, fn, arc_info={"color":color, 'plot':False})
    cube.add_all_cut_locus_points(point_info={'color':'black','s':10})

    cube.plot_faces()