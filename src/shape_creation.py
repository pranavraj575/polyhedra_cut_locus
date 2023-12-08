import numpy as np
from src.shapes import *


class Cube(Shape):
    def __init__(self):
        """
        makes cube where faces 0,1,2,3 are the 'sides' in order (1 is to the right of 0)
        4 is the top, aligned with 0 on the right, 1 on the top
        5 is the bottom, aligned with 0 on the right, 1 on the bottom
        """

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
    def __init__(self):
        """
        makes tetrahedron where face 0 borders face 1 on the right, 2 on the left and 3 on the bottom
        3 is upside down triangle for display purposes
        each face has a circumcenter of radius 1
        """
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
    def __init__(self):
        """
        makes octahedron where the top half are faces 0,1,2,3 (1 is to the right of 0)
        bottom half is 4,5,6,7 (4 is below 0, 5 below 1)
        each face has a circumcenter of radius 1
        we will make 0,1,2,3 normal triangles and 4,5,6,7 upside down triangles for visualization purposes
        """
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
    def __init__(self):
        """
        makes icosahedron
        top 5 are faces 0-4 (1 to the right of 0)
        bottom 5 are faces 15-19 (16 to the right of 15, and all upside down)
        middle 10 are faces 5-14 (11 to the right of 10, and alternating upside down)
        (0,1,2,3,4) is above (5,7,9,11,13). (also 5,7,9,11,13 are all down)
        (6,8,10,12,14) is above (15,16,17,18,19)
        """
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
    def __init__(self):
        """
        makes dodecahedron
        0 on the top
        1-5 surrounding 0 (upside down, 2 to the right of 1)
        11 on the bottom (upside down)
        6-10 surrounding 11 (7 to right of 6)
        middle faces zig zag
         1 2 3 4  5
        6 7 8 9 10
        """
        super().__init__()
        for i in range(12):
            self.add_face()
        tau = 2*np.pi
        top = self.faces[0]

        for i in range(5):
            curr = self.faces[1 + i]
            theta = -tau/4 + i*tau/5
            top.add_boundary_paired(curr, rowtation(theta), 1, -coltation(theta), rotation_T(-i*tau/5), coltation(tau/4))
        for i in range(5):
            curr = self.faces[1 + i]
            neigh = self.faces[1 + (i + 1)%5]
            curr.add_boundary_paired(neigh, rowtation(tau/4 - tau/5), 1, -coltation(tau/4 - tau/5), rotation_T(tau/2 + tau*2/5), coltation(tau/4 + tau/5))

        bottom = self.faces[11]

        for i in range(5):
            curr = self.faces[6 + i]
            theta = tau/4 - i*tau/5
            bottom.add_boundary_paired(curr, rowtation(theta), 1, -coltation(theta), rotation_T(i*tau/5), coltation(-tau/4))
        for i in range(5):
            curr = self.faces[6 + i]
            neigh = self.faces[6 + (i + 1)%5]
            curr.add_boundary_paired(neigh, rowtation(-tau/4 + tau/5), 1, -coltation(-tau/4 + tau/5), rotation_T(tau/2 - tau*2/5), coltation(-tau/4 - tau/5))

        for i in range(5):
            floor = self.faces[6 + i]
            ceil = self.faces[1 + i]
            next_floor = self.faces[6 + (i + 1)%5]
            floor.add_boundary_paired(ceil, rowtation(-tau/4 + 2*tau/5), 1, -coltation(-tau/4 + 2*tau/5), np.identity(2), coltation(tau/4 + tau*2/5))
            ceil.add_boundary_paired(next_floor, rowtation(tau/4 - 2*tau/5), 1, -coltation(tau/4 - 2*tau/5), np.identity(2), coltation(-tau/4 - 2*tau/5))

    def faces_to_plot_n_m(self):
        def face_map(i, j):
            if i == 0 and j == 2:
                return self.faces[0]
            if i == 1:
                return self.faces[j + 1]
            if i == 2:
                return self.faces[j + 6]
            if i == 3 and j == 0:
                return self.faces[11]
            return None

        return face_map, 4, 5


class Antiprism(Shape):
    def __init__(self, n):
        """
        makes uniform antiprism with n-gon (2n+2 faces)
        0 on the top
        1 on the bottom
        2 to n+1 are upside down triangles surrounding 0 (3 to right of 2)
        n+2 to 2n+1 are triangles surrounding 1 (n+3 to right of n+2)
        middle faces zig zag
        (n+1)      2       3 ...
             (n+2)   (n+3)   ...
        2 is at bottom of 0
        n+2 is at top of 1

        :param n: n-gon for root of antiprism
        """
        super().__init__()
        assert n >= 2
        if n == 2:
            print("WARNING: 2 antiprism will display/act weird because of the line face")
        self.n_gon = n
        for i in range(2*n + 2):
            self.add_face()
        dtheta = 2*np.pi/n
        top = self.faces[0]
        bottom = self.faces[1]
        r = np.sqrt(3)/np.tan(np.pi/n)

        for i in range(n):
            curr = self.faces[i + 2]
            theta = -np.pi/2 + i*dtheta
            top.add_boundary_paired(curr, rowtation(theta), r, -r*coltation(theta), rotation_T(-i*dtheta), coltation(np.pi/2))

        for i in range(n):
            curr = self.faces[n + i + 2]
            theta = np.pi/2 - i*dtheta
            bottom.add_boundary_paired(curr, rowtation(theta), r, -r*coltation(theta), rotation_T(i*dtheta), coltation(-np.pi/2))

        for i in range(n):
            floor = self.faces[i + 2 + n]
            ceil = self.faces[i + 2]
            next_floor = self.faces[(i + 1)%n + 2 + n]
            floor.add_boundary_paired(ceil, rowtation(np.pi/6), 1, -coltation(np.pi/6), np.identity(2), coltation(-np.pi*5/6))
            ceil.add_boundary_paired(next_floor, rowtation(-np.pi/6), 1, -coltation(-np.pi/6), np.identity(2), coltation(np.pi*5/6))

    def faces_to_plot_n_m(self):
        center = self.n_gon//2

        def face_map(i, j):
            if i == 0 and j == center:
                return self.faces[0]
            if i == 3 and j == center - ((self.n_gon + 1)%2):
                return self.faces[1]
            if i == 1:
                return self.faces[2 + (j + center)%self.n_gon]
            if i == 2:
                return self.faces[self.n_gon + 2 + (j + center + 1)%self.n_gon]
            return None

        return face_map, 4, self.n_gon


class NTorus(Shape):
    def __init__(self, n):
        """
        makes n-torus, single face

        :param n: dimension of torus
        """
        super().__init__()
        self.add_face()
        face = self.faces[0]
        face: Face
        for i in range(n):
            ei = np.zeros((n, 1))
            ei[i, 0] = 1
            face.add_boundary_paired(face, ei.T, 1, -ei, np.identity(n), -ei)


class LargeNTorus(Shape):
    def __init__(self, n):
        """
        makes n-torus, 2^n faces, each dimension is 2 faces long
        """
        super().__init__()

        def name(i):
            """
            binary name of face i
            """
            name = bin(i)[2:]
            while len(name) < n:
                name = '0' + name
            return name[::-1]

        for i in range(int(2**n)):
            self.add_face(Face(name=name(i)))
        for i in range(int(2**n)):
            nm = name(i)
            for k in range(n):
                neigh_nm = nm[:k] + str((int(nm[k]) + 1)%2) + nm[k + 1:]
                ek = np.zeros((n, 1))
                ek[k][0] = 1

                self.faces[nm].add_boundary_paired(self.faces[neigh_nm], ek.T, 1, -ek, np.identity(n), -ek)


class Large2Torus(LargeNTorus):
    def __init__(self):
        """
        makes 2-torus, 4 faces for plotting reasons

        0 1
        2 3
        """
        super().__init__(2)

    def faces_to_plot_n_m(self):
        def face_map(i, j):
            return self.faces[str(j) + str(i)]

        return face_map, 2, 2


if __name__ == "__main__":
    from display_utils import *

    cube = Large2Torus()
    p = np.array([[0.0], [0.0]])
    # cube.add_point_to_face(p, '00', {'color':'black', 's':20})

    cube.interactive_vornoi_plot(event_key='motion_notify_event')
    quit()
    top = cube.faces[4]
    bottom = cube.faces[5]
    top: Face
    for path in top.face_paths_to(bottom.name):
        for (_, f) in path:
            print(f.name, end=', ')
        print()
    cube.plot_faces()

    face: Face = cube.faces[0]
    fn = 0
    p = np.array([[-1.0], [0.0]])

    cube.add_point_to_face(p, fn, {'color':'red', 's':20})

    radii = [1 + i/11 for i in range(115)]
    for r, color in zip(radii, rb_gradient(len(radii))):
        A = Arc(p, 0, np.pi*2, r)
        cube.add_arc_end_to_face(A, fn, arc_info={"color":color, 'plot':False})
    cube.add_all_cut_locus_points(point_info={'color':'black', 's':2})

    cube.plot_faces(show=True, legend=lambda i, j:i == 1 or i == 2)
