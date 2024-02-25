from src.shapes import *


class Prism(Shape):
    def __init__(self, n, tolerance=.001):
        """
        makes n-gon prism where faces 0 to (n-1) are the 'sides' in order (1 is to the right of 0)
        n is the top, aligned with 0 on the bottom
        (n+1) is the bottom, aligned with 0 on the top
        """

        super().__init__(tolerance=tolerance)
        if n == 2:
            raise Exception("Use Mirror(4) instead of Prism(2)")
        if n < 2 or type(n) is not int:
            raise Exception("n=" + str(n) + " is not a valid Prism")
        self.n = n
        for i in range(n + 2):
            self.add_face()
        I = np.identity(2)
        e1 = np.array([[1], [0]])
        e2 = np.array([[0], [1]])
        for i in range(n):
            curr: Face = self.faces[i]
            neigh = self.faces[(i + 1)%n]
            curr.add_boundary_paired(neigh, e1.T, 1, -e1, I, -e1)

        top = self.faces[n]
        bot = self.faces[n + 1]
        r = 1/np.tan(np.pi/n)

        for i in range(n):
            curr = self.faces[i]
            theta = -np.pi/2 + 2*np.pi*i/n
            # "angle" of boundary of top face.
            top.add_boundary_paired(curr, rowtation(theta), r, -r*coltation(theta), rotation_T(-2*np.pi*i/n), e2)
        for i in range(n):
            curr = self.faces[i]
            theta = np.pi/2 - 2*np.pi*i/n
            # "angle" of boundary of bottom face.
            bot.add_boundary_paired(curr, rowtation(theta), r, -r*coltation(theta), rotation_T(2*np.pi*i/n), -e2)

    def faces_to_plot_n_m(self):
        def face_map(i, j):
            shift = self.n//2
            if i == 1:
                return self.faces[(j - shift)%self.n]
            if i == 0 and j == shift:
                return self.faces[self.n]
            if i == 2 and j == shift:
                return self.faces[self.n + 1]
            return None

        return face_map, 3, self.n


class Pyramid(Shape):
    def __init__(self, n, tolerance=.001):
        """
        makes pyramid with an n-gon base (3<=n<=6)
        face n is the bottom
        face 0 to n-1 are triangles, with face 1 to the right of face 0

        the top of face n borders the bottom of face 0
        """
        super().__init__(tolerance=tolerance)
        if n == 6:
            raise Exception("Use Mirror(6) instead of Pyramid(6)")
        if n <= 2 or n > 6 or type(n) is not int:
            raise Exception("n=" + str(n) + " is not a valid Pyramid")
        self.n = n
        for i in range(n + 1):
            self.add_face()
        bottom = self.faces[n]

        for i in range(n):
            curr = self.faces[i]
            neigh = self.faces[(i + 1)%n]
            curr.add_boundary_paired(neigh, rowtation(np.pi/6), 1, -coltation(np.pi/6), rotation_T(-np.pi/3), coltation(np.pi*5/6))
        r = np.sqrt(3)/np.tan(np.pi/n)
        for i in range(n):
            curr = self.faces[i]
            theta = np.pi/2 - 2*i*np.pi/n
            curr.add_boundary_paired(bottom, np.array([[0, -1]]), 1, np.array([[0], [1]]), rotation_T(-i*2*np.pi/n), r*coltation(theta))

    def _tetrahedron_faces_to_plot_n_m(self):
        def face_map(i, j):
            if i == 1 and j == 1:
                return self.faces[3]
            if i == 0:
                return self.faces[(j + 2)%3]
            return None

        return face_map, 2, 3

    def _large_prism_faces_to_plot_n_m(self):
        shift = self.n//2

        def face_map(i, j):
            if i == 0:
                return self.faces[(j - shift)%self.n]
            if i == 1 and j == shift:
                return self.faces[self.n]
            return None

        return face_map, 2, self.n

    def faces_to_plot_n_m(self):
        if self.n == 3:
            return self._tetrahedron_faces_to_plot_n_m()
        return self._large_prism_faces_to_plot_n_m()


class Bipyramid(Shape):
    def __init__(self, n, tolerance=.001):
        """
        makes bipyramid with 2n triangles(3<=n<=6)
        faces 0 to (n-1) are on the top (1 to the right of 0)
        faces n to 2n-1 are on the bottom (n+1 to the right of n)

        face i is above face i+n
        """
        super().__init__(tolerance=tolerance)
        if n == 6:
            raise Exception("Use Mirror(6) instead of Bipyramid(6)")
        if n < 2 or n > 6 or type(n) is not int:
            raise Exception("n=" + str(n) + " is not a valid Bipyramid")
        self.n = n
        for i in range(2*n):
            self.add_face()

        for i in range(n):
            curr = self.faces[i]
            neigh = self.faces[(i + 1)%n]
            curr.add_boundary_paired(neigh, rowtation(np.pi/6), 1, -coltation(np.pi/6), rotation_T(-np.pi/3), coltation(np.pi*5/6))

        for i in range(n):
            curr = self.faces[n + i]
            neigh = self.faces[n + (i + 1)%n]
            curr.add_boundary_paired(neigh, rowtation(-np.pi/6), 1, -coltation(-np.pi/6), rotation_T(np.pi/3), coltation(-np.pi*5/6))

        for i in range(n):
            top = self.faces[i]
            bot = self.faces[i + self.n]
            top.add_boundary_paired(bot, np.array([[0, -1]]), 1, np.array([[0], [1]]), np.identity(2), np.array([[0], [1]]))

    def faces_to_plot_n_m(self):
        def face_map(i, j):
            return self.faces[i*self.n + j]

        return face_map, 2, self.n


class Tetrahedron(Pyramid):
    def __init__(self, tolerance=.001):
        """
        special case of pyramid
        makes tetrahedron where face 0 borders face 1 on the right, 2 on the left and 3 on the bottom
        3 is upside down triangle for display purposes
        each face has a circumcenter of radius 1
        """
        super().__init__(3, tolerance=tolerance)


class Cube(Prism):
    def __init__(self, tolerance=.001):
        """
        special case, square prism
        makes cube where faces 0,1,2,3 are the 'sides' in order (1 is to the right of 0)
        4 is the top, aligned with 0 on the bottom
        5 is the bottom, aligned with 0 on the top
        """

        super().__init__(4, tolerance=tolerance)


class Octahedron(Bipyramid):
    def __init__(self, tolerance=.001):
        """
        special case of bipyramid
        makes octahedron where the top half are faces 0,1,2,3 (1 is to the right of 0)
        bottom half is 4,5,6,7 (4 is below 0, 5 below 1)
        each face has a circumcenter of radius 1
        we will make 0,1,2,3 normal triangles and 4,5,6,7 upside down triangles for visualization purposes
        """
        super().__init__(4, tolerance=tolerance)


class Icosahedron(Shape):
    def __init__(self, tolerance=.001):
        """
        makes icosahedron
        top 5 are faces 0-4 (1 to the right of 0)
        bottom 5 are faces 15-19 (16 to the right of 15, and all upside down)
        middle 10 are faces 5-14 (11 to the right of 10, and alternating upside down)
        (0,1,2,3,4) is above (5,7,9,11,13). (also 5,7,9,11,13 are all down)
        (6,8,10,12,14) is above (15,16,17,18,19)
        """
        super().__init__(tolerance=tolerance)

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
    def __init__(self, tolerance=.05):
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
        super().__init__(tolerance=tolerance)
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
    def __init__(self, n, tolerance=.001):
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
        super().__init__(tolerance=tolerance)
        if n == 2:
            print("WARNING: 2 antiprism will display/act weird because of the extra line face")
        if n < 2 or type(n) is not int:
            raise Exception("n=" + str(n) + " is not a valid Antiprism")
        self.n = n
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
        center = self.n//2

        def _face_map(i, j):
            if i == 0 and j == center:
                return self.faces[0]
            if i == 3 and j == center - ((self.n + 1)%2):
                return self.faces[1]
            if i == 1:
                return self.faces[2 + (j + center)%self.n]
            if i == 2:
                return self.faces[self.n + 2 + (j + center + 1)%self.n]
            return None

        def face_map(i, j):
            if self.n > 2:
                return _face_map(i, j)
            else:
                return _face_map(i + 1, j)

        return face_map, 4 if self.n > 2 else 2, self.n


class ElongatedBipyramid(Shape):
    def __init__(self, n, tolerance=.001):
        """
        makes elongated bipyramid with 2n triangles(3<=n<=6) and n squares
        triangle faces 0 to (n-1) are on the top (1 to the right of 0)
        square faces n to 2n-1 are on the bottom (n+1 to the right of n)
        triangle faces 2n to 3n-1 are on the bottom (n+1 to the right of n)

        face i is above face i+n

        triangle faces have an inscribed circle radius 1
        thus, square side lengths are 2*sqrt(3)
        """
        super().__init__(tolerance=tolerance)
        if n == 6:
            raise Exception("Use Prism(6) instead of ElongatedBipyramid(6)")
        if n < 2 or n > 6 or type(n) is not int:
            raise Exception("n=" + str(n) + " is not a valid Elongated Bipyramid")
        self.n = n
        for i in range(3*n):
            self.add_face()
        I = np.identity(2)
        e1 = np.array([[1], [0]])
        e2 = np.array([[0], [1]])

        square_r = np.sqrt(3)

        for i in range(n):
            curr = self.faces[i]
            neigh = self.faces[(i + 1)%n]
            curr.add_boundary_paired(neigh, rowtation(np.pi/6), 1, -coltation(np.pi/6), rotation_T(-np.pi/3), coltation(np.pi*5/6))

        for i in range(n):
            curr = self.faces[n + i]
            neigh = self.faces[n + (i + 1)%n]
            curr.add_boundary_paired(neigh, e1.T, square_r, -square_r*e1, I, -square_r*e1)

        for i in range(n):
            curr = self.faces[2*n + i]
            neigh = self.faces[2*n + (i + 1)%n]
            curr.add_boundary_paired(neigh, rowtation(-np.pi/6), 1, -coltation(-np.pi/6), rotation_T(np.pi/3), coltation(-np.pi*5/6))

        for i in range(n):
            top = self.faces[i]
            mid = self.faces[n + i]
            bot = self.faces[2*n + i]
            top.add_boundary_paired(mid, -e2.T, 1, e2, I, square_r*e2)
            mid.add_boundary_paired(bot, -e2.T, square_r, square_r*e2, I, e2)

    def faces_to_plot_n_m(self):
        def face_map(i, j):
            return self.faces[i*self.n + j]

        return face_map, 3, self.n


class ElongatedPyramid(Shape):
    def __init__(self, n, tolerance=.001):
        """
        makes elongated pyramid with n triangles(3<=n<=6), n squares, and one n-gon
        triangle faces 0 to (n-1) are on the top (1 to the right of 0)
        square faces n to 2n-1 are on the bottom (n+1 to the right of n)
        face i is above face i+n

        n-gon is face 2n on the bottom (face n is above face 2n)

        square faces have side length 2
        thus, triangle faces have an inscribed circle of 1/sqrt(3)
        n-gon has a face with inscribed r = 1/tan(pi/n)
        """
        super().__init__(tolerance=tolerance)
        if n == 6:
            raise Exception("Use Prism(6) instead of ElongatedPyramid(6)")
        if n < 2 or n > 6 or type(n) is not int:
            raise Exception("n=" + str(n) + " is not a valid Elongated Pyramid")
        self.n = n
        for i in range(2*n + 1):
            self.add_face()
        I = np.identity(2)
        e1 = np.array([[1], [0]])
        e2 = np.array([[0], [1]])
        triangle_r = 1/np.sqrt(3)
        n_gon_r = 1/np.tan(np.pi/n)
        for i in range(n):
            curr = self.faces[i]
            neigh = self.faces[(i + 1)%n]
            curr.add_boundary_paired(neigh, rowtation(np.pi/6), triangle_r, -triangle_r*coltation(np.pi/6), rotation_T(-np.pi/3), triangle_r*coltation(np.pi*5/6))

        for i in range(n):
            curr = self.faces[n + i]
            neigh = self.faces[n + (i + 1)%n]
            curr.add_boundary_paired(neigh, e1.T, 1, -e1, I, -e1)

        for i in range(n):
            top = self.faces[i]
            bot = self.faces[n + i]
            top.add_boundary_paired(bot, -e2.T, triangle_r, triangle_r*e2, I, e2)

        bottom = self.faces[2*n]
        for i in range(n):
            curr = self.faces[n + i]
            theta = np.pi/2 - 2*np.pi*i/n
            # "angle" of boundary of bottom face.
            bottom.add_boundary_paired(curr, rowtation(theta), n_gon_r, -n_gon_r*coltation(theta), rotation_T(2*np.pi*i/n), -e2)

    def faces_to_plot_n_m(self):
        shift = self.n//2

        def face_map(i, j):
            if i <= 1:
                return self.faces[i*self.n + (j - shift)%self.n]
            elif j == shift:
                return self.faces[2*self.n]
            return None

        return face_map, 3, self.n


class Mirror(Shape):
    def __init__(self, n, tolerance=.001):
        """
        makes n-gon 'mirror', two n-gons pasted to each other

        :param n: n-gon
        """
        super().__init__(tolerance=tolerance)
        if n <= 2 or type(n) is not int:
            raise Exception("n=" + str(n) + " is not a valid Mirror")
        self.n = n
        for i in range(2):
            self.add_face()
        front = self.faces[0]
        back = self.faces[1]

        for i in range(n):
            theta = -np.pi/2 + 2*np.pi*i/n
            # "angle" of boundary of front face.
            phi = -np.pi - theta
            # "angle" of boundary of back face.
            front.add_boundary_paired(back, rowtation(theta), 1, -coltation(theta), rotation_T(phi - theta + np.pi), coltation(phi))

    def is_polyhedra(self):
        return False


class NTorus(Shape):
    def __init__(self, n, tolerance=.001):
        """
        makes n-torus, single face

        :param n: dimension of torus
        """
        super().__init__(tolerance=tolerance)
        self.add_face()
        face = self.faces[0]
        face: Face
        for i in range(n):
            ei = np.zeros((n, 1))
            ei[i, 0] = 1
            face.add_boundary_paired(face, ei.T, 1, -ei, np.identity(n), -ei)

    def is_polyhedra(self):
        return False


class LargeNTorus(Shape):
    def __init__(self, n, tolerance=.001):
        """
        makes n-torus, 2^n faces, each dimension is 2 faces long
        """
        super().__init__(tolerance=tolerance)

        def name(i):
            """
            binary name of face i
            """
            name = bin(i)[2:]
            while len(name) < n:
                name = '0' + name
            return name[::-1]

        for i in range(int(2**n)):
            self.add_face(Face(name=name(i), tolerance=self.tol))
        for i in range(int(2**n)):
            nm = name(i)
            for k in range(n):
                neigh_nm = nm[:k] + str((int(nm[k]) + 1)%2) + nm[k + 1:]
                ek = np.zeros((n, 1))
                ek[k][0] = 1

                self.faces[nm].add_boundary_paired(self.faces[neigh_nm], ek.T, 1, -ek, np.identity(n), -ek)

    def is_polyhedra(self):
        return False


class Large2Torus(LargeNTorus):
    def __init__(self, tolerance=.001):
        """
        makes 2-torus, 4 faces for plotting reasons

        0 1
        2 3
        """
        super().__init__(2, tolerance=tolerance)

    def faces_to_plot_n_m(self):
        def face_map(i, j):
            return self.faces[str(j) + str(i)]

        return face_map, 2, 2


if __name__ == "__main__":

    cube = Pyramid(5)
    p = np.array([[0.0], [0.0]])
    # cube.add_point_to_face(p, '00', {'color':'black', 's':20})

    cube.interactive_vornoi_plot(event_key='motion_notify_event', legend=lambda i, j: True)
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

    cube.add_point_to_face(p, fn, {'color': 'red', 's': 20})

    radii = [1 + i/11 for i in range(115)]
    for r, color in zip(radii, rb_gradient(len(radii))):
        A = Arc(p, 0, np.pi*2, r)
        cube.add_arc_end_to_face(A, fn, arc_info={"color": color, 'plot': False})
    cube.add_all_cut_locus_points(point_info={'color': 'black', 's': 2})

    cube.plot_faces(show=True, legend=lambda i, j: i == 1 or i == 2)
