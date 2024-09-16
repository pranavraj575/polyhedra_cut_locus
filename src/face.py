import itertools
import numpy as np

from src.bound import Bound


class Face:
    def __init__(self, name, tolerance, bounds_faces=None, basepoint=None):
        """
        Creates a face
        :param name: name of face
        :param tolerance: tolerance for methods like 'within'
        :param bounds_faces: list of (Bound,Face) to initialize boundaries
            Note: usually start with an empty face, then create bounds later
        :param basepoint: point inside face, if None, uses 0
            this is basically only used for ordering vertices
        """
        self.name = name
        self.bounds = []
        self.tol = tolerance
        self.vertices = None
        self.dimension = None
        self.bound_M = None
        self.bound_b = None
        self.double_face_edge = []

        if bounds_faces is not None:
            for (bound, F) in bounds_faces:
                self.add_boundary(bound, F)
        self.basepoint = basepoint  # if None, will update this when the dimension is set by adding a bound

    def __str__(self):
        return '(Face ' + str(self.name) + ': connected to faces ' + str([F.name for (_, F) in self.bounds]) + ')'

    def add_boundary(self, bound: Bound, F, update=True):
        """
        adds boundary to face F into self
        :param bound: Bound
        :param update: whether to update internal bound arrays and vertices
        :param F: Face
        """
        if F in [Fp for (_, Fp) in self.bounds]:
            self.double_face_edge.append(F)
        self.bounds.append((bound, F))
        self.dimension = bound.check_valid(self.dimension)
        if self.basepoint is None:
            self.basepoint = np.zeros((self.dimension, 1))
        if update:
            self._create_bound_arrays()
            self._create_vertices()

    def _order_vertices(self):
        """
        orders the vertices by angle (so that a path through all of them is the bounds of the face)
            Note: only valid if 2 dimensional
        """
        if self.dimension != 2:
            return
        self.vertices.sort(
            key=lambda v: (2*np.pi + np.arctan2((v[0] - self.basepoint)[1][0], (v[0] - self.basepoint)[0][0]))%(
                    2*np.pi))

    def get_path_and_faces(self):
        """
        grabs path of vertices (v0,v1),(v1,v2),...,(vn,v0)
        takes order 'around' the face
        :return: ((v,v'),(Bound, Face)) list with v,v' column vectors
        """
        if not self.dimension == 2:
            raise Exception("this will not work in !=2 dimensions")
        path = self.get_vertices()
        out = []
        for i in range(len(path)):
            v1, rows = path[i]
            v2, rowsp = path[(i + 1)%len(path)]
            row = None
            for r in rows:
                if r in rowsp:
                    row = r
            out.append(((v1, v2), self.bounds[row]))
            # out.append(((v1, v2), self.bounds[row][1]))
        return out

    def within_bounds(self, p):
        """
        returns if point p is inside face
        :param p: column vector (np array of dimension (self.dimension,1))
        :return: boolean
        """
        if self.bound_M is None:
            self._create_bound_arrays()
        return np.all(self.bound_M@p <= self.bound_b + self.tol)

    def bound_of_face(self, F):
        """
        returns the bound corresponding with face F, None if non existant
        :param F: Face
        :return: Bound
        """
        for (bound, Fp) in self.bounds:
            if Fp.name == F.name:
                return bound
        raise Exception("BOUND NOT FOUND WITH THIS FACE")

    def _create_bound_arrays(self):
        """
        creates arrays M and b of all bounds
        for n bounds, M is of dimension (n,self.dimension) and b is of dimension (n,1)
        a point p is within this face if Mp<=b
        """
        self.bound_M = np.array([bound.m[0] for (bound, _) in self.bounds])
        self.bound_b = np.array([[bound.b] for (bound, _) in self.bounds])

    def _create_vertices(self):
        """
        creates all vertices of the face
        vertices are a list of (vertex: column vector, indices of bounds that create it: tuple)
        """
        self.vertices = []
        if self.bound_M is None:
            self._create_bound_arrays()
        n = len(self.bounds)
        if n < self.dimension:
            # print("WARNING: no vertices since not enough boundaries")
            return
        for rows in itertools.combinations(range(n), self.dimension):
            sub_M = self.bound_M[rows, :]
            sub_b = self.bound_b[rows, :]
            if abs(np.linalg.det(sub_M)) > self.tol:
                vertex = np.linalg.inv(sub_M)@sub_b
                if self.within_bounds(vertex):
                    self.vertices.append((vertex, rows))
        self._order_vertices()

    def get_vertices(self):
        """
        grabs all vertices of the face
        :return: list of (vertex: column vector, indices of bounds that create it: tuple)
        """
        return self.vertices

    def get_closest_point(self, p):
        """
        returns the closest point in the face to p
        :param p: column vector (np array of dimension (self.dimension,1))
        :return: column vector (np array of dimension (self.dimension,1))
        """
        if self.within_bounds(p):
            return p
        q = p.copy()
        small_bound = None  # first bound to p from basepoint
        for (bound, F) in self.bounds:
            if not bound.within(q, tol=self.tol):
                # goes from inside bound to outside bound
                q = bound.grab_intersection(self.basepoint, q)
                small_bound = bound
        m = small_bound.m
        mbar = m/np.linalg.norm(m)
        center = m.T*small_bound.b/np.square(np.linalg.norm(m))

        offset = p - mbar.T*np.dot(mbar, p)
        proj = offset + center

        exiting = self.get_exit_point(center, offset)

        if exiting is None:
            return proj
        return exiting

    def get_exit_point(self, p, v):
        """
        returns the first point that a ray starting from p and going to v exits face
            None if does not exist
        :param p: column vector (np array of dimension (self.dimension,1))
        :param v: column vector (np array of dimension (self.dimension,1))
        :return: column vector (np array of dimension (self.dimension,1)) or None
        """
        q = p + v
        # q represents the point on the line pq that is intersecting the closest boundary
        if self.within_bounds(q):
            return None
        for (bound, F) in self.bounds:
            if bound.within(p, tol=self.tol) and not bound.within(q, tol=self.tol):
                # goes from inside bound to outside bound
                q = bound.grab_intersection(p, v)
        if not self.within_bounds(q):
            # in this case, we have a line that ends outside the face, and never enters the face
            return None
        return q

    def line_within_bounds(self, p, q, ignore_points):
        """
        returns if any part of the line p->q is within the face
        :param p: column vector (np array of dimension (self.dimension,1))
        :param q: column vector (np array of dimension (self.dimension,1))
        :return: boolean
        """
        if ignore_points:
            if not self.within_bounds(p):
                p = self.get_exit_point(q, p - q)
            if p is None:
                # if p is not within bounds and the line q->p does not exit face, this line is not in face
                return False
            if not self.within_bounds(q):
                q = self.get_exit_point(p, q - p)
            if q is None:
                return False
            if q is None:
                return False
            # check if the two endpoints are far enough apart
            return np.linalg.norm(p - q) > self.tol
        else:
            # enough to test whether either of the endpoints are within the bounds or if the line exits the face
            return self.within_bounds(p) or self.within_bounds(q) or (self.get_exit_point(p, q - p) is not None)

    def get_ray_within_bounds(self, p, direction):
        """
        returns the part of the ray starting at p and going in driection in face
        :param p: column vector (np array of dimension (self.dimension,1))
        :return: None or (p',q')
        """
        while self.within_bounds(direction + p):
            direction *= 2
        return self.get_segment_within_bounds(p, direction + p)

    def get_segment_within_bounds(self, p, q):
        """
        returns the part of the line p->q within the face
        :param p: column vector (np array of dimension (self.dimension,1))
        :param q: column vector (np array of dimension (self.dimension,1))
        :return: None or (p',q')
        """
        if self.within_bounds(p) and self.within_bounds(q):
            return (p, q)
        if self.within_bounds(p):
            qp = self.get_exit_point(p, q - p)
            return (p, qp)
        if self.within_bounds(q):
            pp = self.get_exit_point(q, p - q)
            return (pp, q)
        qp = self.get_exit_point(p, q - p)
        pp = self.get_exit_point(q, p - q)
        if qp is None and pp is None:
            return None
        # weird error due to tolerance
        if qp is None:
            qp = pp
        if pp is None:
            pp = qp
        return (pp, qp)

    def add_boundary_paired(self, f2, m, b, s, T, si):
        """
        adds boundary to face f2, and corresponding bound to self
        :param f2: other face that this bound connects to
        :param m: argument of Bound
        :param b: argument of Bound
        :param s: argument of Bound
        :param T: argument of Bound
        :param si: argument of Bound
        """
        B1 = Bound(m, b, s, T, si, dimension=self.dimension, identifier=str(self.name) + str(f2.name))
        self.add_boundary(B1, f2)
        f2.add_boundary(B1.get_inverse_bound(), self)

    def face_paths_to(self, fn, visited_names=None, diameter=None):
        """
        returns all paths to specified face using DFS

        :param fn: name of target face
        :param visited_names: set of faces we have already visited
        :param diameter: longest path of faces to consider (None if infinite)
        :return: (Bound,Face) list of 'edges' and 'next Faces'
        """
        if visited_names is None:
            visited_names = {self.name}
        visited_names.add(self.name)
        if self.name == fn:
            yield []
        elif diameter is not None and diameter <= 0:
            return
        else:
            for (bound, f) in self.bounds:
                if not f.name in visited_names:
                    for path in f.face_paths_to(fn, visited_names=visited_names.copy(),
                                                diameter=None if diameter is None else diameter - 1):
                        yield [(bound, f)] + path
