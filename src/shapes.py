import numpy as np, itertools
from matplotlib import pyplot as plt
from scipy.spatial import Voronoi
from src.my_vornoi import voronoi_plot_2d

TOL = .001


# rotation matrices
def rotation_T(theta):
    """
    :param theta: angle
    :return: 2x2 rotation matrix T by angle theta (np array)
    """
    # T where T @ v is v rotatied by theta
    return np.array([[np.cos(theta), -np.sin(theta)],
                     [np.sin(theta), np.cos(theta)]])


def rotation_of_matrix(T):
    """
    :param T: 2x2 rotation matrix by angle theta (np array)
    :return: theta
    """
    # inverse of rotation_T
    cos = T[0][0]
    sin = T[1][0]
    return np.arctan2(sin, cos)


def coltation(theta):
    """
    returns 2d unit column vector towards theta
    :param theta: angle
    :return: <<cos theta>,<sin theta>>
    """
    return rotation_T(theta)[:, [0]]


def rowtation(theta):
    """
    returns 2d unit row vector towards theta
    :param theta: angle
    :return: <<cos theta,sin theta>>
    """
    return coltation(theta).T


def flatten(L):
    """
    flattens a list of lists
    :param L: list of lists
    :return: L flattened
    """
    out = []
    for l in L:
        out += l
    return out


class Arc:
    def __init__(self, p, low, high, r, distance=None, tolerance=TOL):
        """
        Defines an arc in 2d space
        :param p: point of origin, column vector (np array of dimension (self.dimension,1))
        :param low: low angle in radians
        :param high: high angle in radians (must be >= low)
        :param r: radius, scalar
        :param distance: original distance, or used to keep track of arbitrary data (can be different if we 'diffract')
        :param tolerance: tolerance for intersection and containment
        """
        assert p.shape == (2, 1)
        assert low <= high

        self.tol = tolerance
        self.p = p
        self.r = r
        if distance is None:
            distance = self.r
        self.dist = distance
        diff = high - low
        self.low = low%(2*np.pi)
        self.high = self.low + diff

    def avg_v(self):
        """
        returns vector going towards center of arc
        :return: column vector (np array of dimension (self.dimension,1))
        """
        theta = (self.low + self.high)/2
        return coltation(theta)

    def is_empty(self):
        """
        returns if the arc has empty interior
        :return: boolean
        """
        return (self.high - self.low) <= self.tol

    def angle_of(self, q):
        """
        returns angle of point q wrt arc
        :param q: column vector (np array of dimension (self.dimension,1))
        :return: angle in radians
        """
        x, y = tuple((q - self.p).flatten())
        theta = np.arctan2(y, x)%(np.pi*2)
        return theta

    def _strict_contains_angle(self, theta):
        """
        returns whether theta is in the interior of the arc angles
        :param theta: angle in radians
        :return: boolean
        """
        if (self.low < theta + 2*np.pi and theta + 2*np.pi < self.high):
            return True, theta + 2*np.pi
        if (self.low < theta - self.tol and theta < self.high - self.tol):
            return True, theta
        return False, np.nan

    def _weak_contains_angle(self, theta):
        """
        returns whether theta is in the closure of the arc angles
        :param theta: angle in radians
        :return: boolean
        """
        if (self.low <= theta + 2*np.pi + self.tol and theta + 2*np.pi <= self.high + self.tol):
            return True, theta + 2*np.pi
        if (self.low <= theta + self.tol and theta <= self.high + self.tol):
            return True, theta
        return False, np.nan

    def contains_angle(self, theta, interior=False):
        """
        returns whether theta is within the arc angles
        assumes theta is on [0,2pi)
        :param theta: angle in radians
        :param interior: whether we use exclusive or inclusive bounds
        :return: boolean
        """
        #
        return self._strict_contains_angle(theta) if interior else self._weak_contains_angle(theta)

    def within_arc(self, q, interior=False):
        """
        returns whether a point q is within the arc
        :param q: column vector (np array of dimension (self.dimension,1))
        :param interior: whether we check in interior or simply the closure
        :return: boolean
        """
        # returns whether q is in the arc
        # interior is strict containment
        r = np.linalg.norm(self.p - q)

        if not interior and r == 0:
            return True
        if r > self.r + self.tol:
            return False

        theta = self.angle_of(q)
        return self.contains_angle(theta, interior=interior)[0]

    def split_arc_on_point(self, q):
        """
        returns two arcs that are this arc split based on q
        assumes q is within the arc
        :param q: column vector (np array of dimension (self.dimension,1))
        :return: partition (Arc,Arc), so that the angle to q is an endpoint of each
        """
        # returns two arcs that are this arc split in parts based on p
        if not self.within_arc(q, interior=True):
            raise Exception("split arc called on point not in interior of arc")
        theta = self.angle_of(q)
        if self.low > theta:
            theta = theta + 2*np.pi
        return (Arc(self.p, self.low, theta, self.r, distance=self.dist, tolerance=self.tol),
                Arc(self.p, theta, self.high, self.r, distance=self.dist, tolerance=self.tol))

    def _break_arc(self, a, b):
        """
        breaks arc into pieces based on line segment a-b
        idea is if arc interacts with a-b, split it into pieces that fully go through a-b and pieces that do not

        :param a: column vector (np array of dimension (self.dimension,1))
        :param b: column vector (np array of dimension (self.dimension,1))
        :return: list of arcs
        """
        for q in (a, b):
            if self.within_arc(q, interior=True):
                (A, B) = self.split_arc_on_point(q)
                return A.break_arc(a, b) + B.break_arc(a, b)
        whole_interx = self.circle_line_segment_intersection(a, b)
        split_angles = [self.low]
        for q in whole_interx:
            theta = self.angle_of(q)
            b, theta = self.contains_angle(theta)
            if b:
                split_angles.append(theta)
        split_angles.append(self.high)
        split_angles.sort()
        if len(split_angles) > 2:
            out = []
            for i in range(len(split_angles) - 1):
                th, ph = split_angles[i:i + 2]
                out.append(Arc(self.p, th, ph, self.r, distance=self.dist, tolerance=self.tol))
            return out
        else:
            return [self]

    def break_arc(self, a, b):
        """
        breaks arc into pieces based on line segment a-b, ignoring empty arcs
        idea is if arc interacts with a-b, split it into pieces that fully go through a-b and pieces that do not

        :param a: column vector (np array of dimension (self.dimension,1))
        :param b: column vector (np array of dimension (self.dimension,1))
        :return: list of arcs
        """
        return [A for A in self._break_arc(a, b) if not A.is_empty()]

    def breakable(self, a, b):
        """
        checks whether arc is breakable on a-b
        :param a: column vector (np array of dimension (self.dimension,1))
        :param b: column vector (np array of dimension (self.dimension,1))
        :return: boolean
        """
        return len(self.break_arc(a, b)) > 1

    def circle_line_segment_intersection(self, pt1, pt2):
        """Find the points where the circle (completed arc) intersects a line segment pt1-pt2.
        This can happen at 0, 1, or 2 points.

        :param pt1: column vector (np array of dimension (self.dimension,1))
        :param pt2: column vector (np array of dimension (self.dimension,1))
        :return: list of numpy arrays (col vectors)
        """
        full_line = False
        (p1x, p1y), (p2x, p2y), (cx, cy) = tuple(pt1.flatten()), tuple(pt2.flatten()), tuple(self.p.flatten())
        (x1, y1), (x2, y2) = (p1x - cx, p1y - cy), (p2x - cx, p2y - cy)
        dx, dy = (x2 - x1), (y2 - y1)
        dr = (dx**2 + dy**2)**.5
        big_d = x1*y2 - x2*y1
        discriminant = self.r**2*dr**2 - big_d**2

        if discriminant < 0:  # No intersection between circle and line
            return []
        else:  # There may be 0, 1, or 2 intersections with the segment
            intersections = [
                (cx + (big_d*dy + sign*(-1 if dy < 0 else 1)*dx*discriminant**.5)/dr**2,
                 cy + (-big_d*dx + sign*abs(dy)*discriminant**.5)/dr**2)
                for sign in ((1, -1) if dy < 0 else (-1, 1))]  # This makes sure the order along the segment is correct
            if not full_line:  # If only considering the segment, filter out intersections that do not fall within the segment
                fraction_along_segment = [(xi - p1x)/dx if abs(dx) > abs(dy) else (yi - p1y)/dy for xi, yi in intersections]
                intersections = [pt for pt, frac in zip(intersections, fraction_along_segment) if 0 <= frac <= 1]
            if len(intersections) == 2 and abs(discriminant) <= self.tol:  # If line is tangent to circle, return just one point (as both intersections have same location)
                (x, y) = intersections[0]
                return [np.array([[x], [y]])]
            else:
                return [np.array([[x], [y]]) for (x, y) in intersections]

    def intersects_with(self, A):
        """
        returns points where self intersects with arc A

        can be 0, 1, or 2 elements
        :param A: Arc
        :return: list of column vectors
        """
        A: Arc
        points = []
        p = self.p
        r0 = self.r
        x0, y0 = tuple(p.flatten())
        q = A.p
        r1 = A.r
        x1, y1 = tuple(q.flatten())
        d = np.linalg.norm(p - q)
        if d <= self.tol:  # starting from same point
            return []
        if d <= abs(r0 - r1):  # one circle is within the other
            return []
        if d > r0 + r1:  # non intersecting
            return []
        if abs(r0 + r1 - d) < self.tol:  # if we are exactly here, just take the midpoint
            points.append((p + q)/2)
        else:
            a = (r0**2 - r1**2 + d**2)/(2*d)
            h = np.sqrt(r0**2 - a**2)
            x2 = x0 + a*(x1 - x0)/d
            y2 = y0 + a*(y1 - y0)/d
            x3 = x2 + h*(y1 - y0)/d
            y3 = y2 - h*(x1 - x0)/d

            x4 = x2 - h*(y1 - y0)/d
            y4 = y2 + h*(x1 - x0)/d

            points = [np.array([[x3], [y3]]), np.array([[x4], [y4]])]
        return [pt for pt in points if self.within_arc(pt) and A.within_arc(pt)]

    def __str__(self):
        return "ARC:{p:" + str(tuple(self.p.flatten())) + "; r:" + str(self.r) + "; dist:" + str(self.dist) + "; range:" + str((self.low, self.high)) + "}"


class Bound:
    def __init__(self, m, b, s, T, si, dimension=None, name=None, identifier=''):
        """
        linear boundary of the form mx<=b
        s,T,si represent the linear transformation that puts a point in our face to the neighboring face
        represented as a shift, rotation matrix, and another shift
        for some point x on this face, s + T x + si = x' where x' is the neighbor face coordinates

        :param m: row vector of bound (mx<=b)
        :param b: scalar of bound (mx<=b)
        :param s: first translation to shift bound to origin (for a point x, x' on new face is si+T(x+s))
        :param T: rotation to shift bound to origin (for a point x, x' on new face is si+T(x+s))
        :param si: final translation to shift bound to origin (for a point x, x' on new face is si+T(x+s))
        :param dimension: dimension of bound, if none, it is set, if value inserted, it is checked
        :param identifier: identifier of bound, should mention face names
        :param name: identifier of bound, or the bound that it is paired with, specify if spawned by another bound
        """
        self.m = m
        self.b = b
        self.s = s
        self.T = T
        self.si = si
        self.dimension = self.check_valid(dimension)
        if name is None:
            base_id = str(tuple(self.m.flatten())) + str(self.b) + str(tuple(self.s.flatten())) + str(tuple(self.T.flatten())) + str(tuple(self.si.flatten())) + str(self.dimension)
            name = base_id + identifier
        self.name = name

    def check_valid(self, dimension):
        """
        Checks if self is valid, returns the correct dimension
        :param dimension: proposed dimension to check
        :return: dimension, raises exception if invalid
        """
        if (len(np.shape(self.m)) != 2 or
                len(np.shape(self.s)) != 2 or
                len(np.shape(self.si)) != 2 or
                len(np.shape(self.T)) != 2 or
                np.shape(self.m)[0] != 1 or
                np.shape(self.s)[1] != 1 or
                np.shape(self.si)[1] != 1
        ):
            raise Exception("ERROR: tried putting in bounds of incorrect dimensions")
        if (np.shape(self.T)[0] != np.shape(self.T)[1] or
                np.shape(self.m)[1] != np.shape(self.s)[0] or
                np.shape(self.s)[0] != np.shape(self.si)[0] or
                np.shape(self.si)[0] != np.shape(self.m)[1] or
                np.shape(self.m)[1] != np.shape(self.T)[0]
        ):
            raise Exception("ERROR: tried putting in bounds of inconsistent dimensions")
        if dimension is None:
            dimension = np.shape(self.T)[0]
        else:
            if dimension != np.shape(self.T)[0]:
                raise Exception("ERROR: inconsistent dimensions")
        if not self.within(np.zeros(dimension)):
            raise Exception("ERROR: zero point needs to be within the face")
        return dimension

    def get_shift(self):
        """
        :return: just the shift part of the bound (s,T,si)
        """
        return (self.s, self.T, self.si)

    def within(self, p, tol=0.):
        """
        checks if point p is within bound with some tolerance
        mp<=b+tol
        :param p: column vector (np array of dimension (self.dimension,1))
        :param tol: scalar
        :return: boolean
        """
        return np.dot(self.m, p) <= self.b + tol

    def grab_intersection(self, p, v):
        """
        find where p+vt intersects the bound (mx<=b)

        :param p: column vector (np array of dimension (self.dimension,1))
        :param v: column vector (np array of dimension (self.dimension,1))
        :return: point where p+vt intersects bound, (np array of dimension (self.dimension,1))
        """
        # m(p+vt)=b
        # simplify to mvt=b-mp
        # then t=(b-mp)/(mv)
        t = (self.b - np.dot(self.m, p))/np.dot(self.m, v)
        # note that v cannot be parallel to m, which makes sense for this to even work
        return p + v*t

    def shift_point(self, x):
        """
        shifts point x to equivalent x' on face according to bound
        :param x: column vector (np array of dimension (self.dimension,1))
        :return: column vector (np array of dimension (self.dimension,1))
        """
        return self.T@(x + self.s) + self.si

    def shift_vec(self, v):
        """
        shifts vector v to equivalent v'  according to bound
        :param v: column vector (np array of dimension (self.dimension,1))
        :return: column vector (np array of dimension (self.dimension,1))
        """
        return self.T@v

    def shift_angle(self, theta):
        """
        shifts angle theta to equivalent theta' according to bound
        :param theta: angle in radians
        :return: angle in radians
        """
        return theta + (rotation_of_matrix(self.T))%(2*np.pi)

    def shift_arc(self, A: Arc):
        """
        shifts arc A to equivalent A' according to bound
        :param A: Arc object
        :return: Arc object
        """
        #
        diff = A.high - A.low
        angle = self.shift_angle(A.low)
        return Arc(self.shift_point(A.p), angle, angle + diff, A.r, distance=A.dist, tolerance=A.tol)

    def get_inverse_bound(self):
        """
        :return: inverse bound, from neighboring face to this face
        """

        Ti = np.linalg.inv(self.T)
        m = -self.m@Ti
        b = -self.b - np.dot(self.m, self.s) - np.dot(self.m@Ti, self.si)
        b = b.flatten()[0]
        return Bound(m, b, -self.si, Ti, -self.s, self.dimension, name=self.name)

    def concatenate_with(self, T=None, s=None):
        """
        given some point x that is translated previously like T' x + s',
            concatenate the translation to this one

            T((T' x + s')+s)+si=T T' x + T(s'+s)+si

        :param T: (self.dimension x self.dimension) translation matrix
        :param s: shift
        :return: (self.dimension x self.dimension) tranlstion, (self.dimension x 1) shift
        """
        if T is None: T = np.identity(self.dimension)
        if s is None: s = np.zeros((self.dimension, 1))
        return (self.T@T, self.si + self.T@(self.s + s))


class Face:
    def __init__(self, name, bounds_faces=None, tolerance=TOL):
        """
        Creates a face
        Note: the point 0 should be inside the face
        :param name: name of face
        :param bounds_faces: list of (Bound,Face) to initialize boundaries
            Note: usually start with an empty face, then create bounds later
        :param tolerance: tolerance for methods like 'within'
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

    def __str__(self):
        return '(Face ' + str(self.name) + ': connected to faces ' + str([F.name for (_, _, F, _) in self.bounds]) + ')'

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
        self.vertices.sort(key=lambda v: (2*np.pi + np.arctan2(v[0][1][0], v[0][0][0]))%(2*np.pi))

    def get_plot_bounds(self):
        """
        returns the bounds of the graphical plot of the face
        :return: (xlim, ylim)
        """
        if self.dimension != 2:
            raise Exception("ERROR: cannot graph a non-2d face")
        xmin, xmax = np.inf, -np.inf
        ymin, ymax = np.inf, -np.inf
        for (v, _) in self.get_vertices():
            x, y = v[:, 0]
            xmin = min(x, xmin)
            xmax = max(x, xmax)
            ymin = min(y, ymin)
            ymax = max(y, ymax)
        return ((xmin, xmax), (ymin, ymax))

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
        return all(self.bound_M@p <= self.bound_b + self.tol)

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
        small_bound = None
        for (bound, F) in self.bounds:
            if not bound.within(q, tol=self.tol):
                # goes from inside bound to outside bound
                q = bound.grab_intersection(np.zeros(q.shape), q)
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
        return q

    def points_on_path(self, p, v):
        """
        returns list of (face, initial point, end point) of a ray starting from p, taking vector v

        :param p: column vector (np array of dimension (self.dimension,1))
        :param v: column vector (np array of dimension (self.dimension,1))
        :return: list of (face, initial point, end point)
        """
        q = p + v
        # q is the end point on current face
        q_ = p + v
        # q_ represents the point on the line pq that is intersecting the closest boundary

        if self.within_bounds(q):
            return [(self, p, q)]
        closest_bound = None
        for (bound, F) in self.bounds:
            bound: Bound
            if bound.within(p, tol=self.tol) and not bound.within(q_, tol=self.tol):
                # goes from inside bound to outside bound
                q_ = bound.grab_intersection(p, v)
                closest_bound = (bound, F)
        if closest_bound is None:
            raise Exception("ERROR: line passes outside of face")
        # now we move to the new face with starting point q_ and vector (q-q_)
        # however, we need to translate this into new face coords
        (bound, F) = closest_bound
        new_p = bound.shift_point(q_)
        new_v = bound.shift_vec(q - q_)
        rest_of_path = F.points_on_path(new_p, new_v)
        return [(self, p, q_)] + rest_of_path

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

    def neighbors(self):
        """
        returns all face neighbors of self
        :return: Face list
        """
        return [f for (_, f) in self.bounds]

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
                    for path in f.face_paths_to(fn, visited_names=visited_names.copy(), diameter=None if diameter is None else diameter - 1):
                        yield [(bound, f)] + path

    def get_vertices_of_face_bound(self, fn):
        """
        returns the vertices that specified face is touching

        :param fn: name of specified face
        :return: list of column vectors that face fn touches
        """
        row = None
        for i, (_, F) in enumerate(self.bounds):
            if F.name == fn:
                row = i
        if row is None:
            raise Exception("ERROR: face " + str(fn) + " not a neighbor")
        out = []
        for v, rows in self.get_vertices():
            if row in rows:
                out.append(v)
        return out

    def push_arc_to_faces(self, A: Arc):
        """
        pushes A to the faces it belongs in
        :param A: Arc
        :return: list of (Face, Arc), list of resulting arcs and their home faces
        """
        arr = np.array([A])
        for ((a, b), (bound, f)) in self.get_path_and_faces():
            arr = [B.break_arc(a, b) for B in arr]
            arr = flatten(arr)
        faces = (self.get_correct_next_face(B) for B in arr)

        arr = [[(B, self)] if F is None else F.push_arc_to_faces(self.bound_of_face(F).shift_arc(B)) for (B, F) in zip(tuple(arr), faces)]
        out = []
        for t in arr:
            out += t
        return out

    def get_correct_next_face(self, A: Arc):
        """
        gets next face that a line within A touches
            Note: uses center line, should be unique if A is a result of 'push_arc_to_faces'

        :param A: Arc
        :return: Face, None if A is within self
        """
        v = A.avg_v()*A.r
        arr = self.points_on_path(A.p, v)
        if len(arr) == 1:
            return None
        else:
            return arr[1][0]


class Shape:
    def __init__(self, faces=None, tolerance=TOL):
        """
        A set of faces, representing a shape
        :param faces: initial set of faces
            Note: usually start empty and add faces as we go
        :param tolerance: tolerance for within_face and such
            adds this to autogenreated faces
        """
        if faces is None:
            faces = dict()
        self.tol = tolerance
        self.faces = {face.name: face for face in faces}
        self.points = {face.name: [] for face in self.faces}
        self.arcs = {face.name: [] for face in self.faces}
        self.memoized_face_translations = dict()
        self.seen_bounds = []
        self.extra_legend = None

    def _pick_new_face_name(self):
        """
        picks name for a new face
        :return: int, first unique face name
        """
        i = len(self.faces)
        while i in self.faces:
            i += 1
        return i

    def _face_exists_correctly(self, fn):
        """
        asserts if face fn exists
        :param fn: face name
        :return: boolean
        """
        return all(fn in dic for dic in (self.faces, self.points, self.arcs))

    def add_face(self, face=None):
        """
        adds new face to shape
        :param face: Face, if already defined
            if None, initializes new face and adds it
        """
        if face is None:
            face = Face(self._pick_new_face_name(), tolerance=self.tol)
        self.faces[face.name] = face
        self.reset_face(face.name)

    def reset_face(self, fn):
        """
        Resets points and arcs of a face to empty
        :param fn: face name
        """
        self.points[fn] = []
        self.arcs[fn] = []

    def reset_all_faces(self):
        """
        resets all faces
        """
        for fn in self.faces:
            self.reset_face(fn)

    def add_point_to_face(self, point, fn, point_info):
        """
        adds point to specified face
        :param point: column vector (np array of dimension (self.dimension,1))
        :param fn: face name
        :param point_info: dictionary of extra information to add to point
        """
        assert self._face_exists_correctly(fn)
        self.points[fn].append((point, point_info))

    def add_points_to_face(self, points, fn, point_info):
        """
        adds points to specified face
        :param points: column vector list
        :param fn: face name
        :param point_info: dictionary of extra information to add to each point
        """
        for point in points:
            self.add_point_to_face(point, fn, point_info=point_info)

    def add_arc_to_face(self, A: Arc, fn):
        """
        adds arc to specified face
        :param A: Arc
        :param fn: face name
        """
        assert self._face_exists_correctly(fn)
        self.arcs[fn].append(A)

    def add_arc_end_to_face(self, A, fn, arc_info=None):
        """
        takes an arc, finds its end arcs and adds them to the correct faces
        :param A: Arc
        :param fn: face name
        :param arc_info: dictionary of extra info to add to arc
        """
        assert self._face_exists_correctly(fn)
        src: Face = self.faces[fn]
        arcs = src.push_arc_to_faces(A)
        for (B, F) in arcs:
            self.arcs[F.name].append((B, arc_info))

    def add_all_cut_locus_points(self, point_info=None, conditional_point_info=None):
        """
        adds cut locus points to all faces
        works best if there are a ton of arcs of various lengths

        :param point_info: dictionary of info to add to cut locus points
        :param conditional_point_info: radius->dictionary of info to add to cut locus points, conditional on r
        """
        for fn in self.faces:
            for (p, r) in self.get_cut_locus_points(fn):
                new_point_info = {k: point_info[k] for k in point_info}
                if conditional_point_info is not None:
                    new_point_info.update(conditional_point_info(r))
                new_point_info.update({'locus_pt': True})
                self.add_point_to_face(p, fn, new_point_info)

    def get_cut_locus_arcs(self, fn):
        """
        Gets paired arcs on specified face that intersect to form points on cut locus
        :param fn: face name
        :return: list of pairs (Arc,Arc)
        """
        arcs = self.arcs[fn]
        out = []
        for (A, _), (B, _) in itertools.combinations(arcs, 2):
            if A.dist == B.dist:
                # equality is fine since we should never update dist
                ps = A.intersects_with(B)
                if ps:
                    out += [(p, A.dist, A, B) for p in ps]

        return [(A, B) for (p, r, A, B) in out if self.is_best_seen(p, r, fn)]

    def _get_voronoi_translations(self, source_fn, sink_fn, diameter=None):
        """
        full version of get_voronoi_translations
        """
        source: Face = self.faces[source_fn]
        translations = []
        for path in source.face_paths_to(sink_fn, diameter=diameter):
            T, s = np.identity(source.dimension), np.zeros((source.dimension, 1))
            for (bound, _) in path:
                bound: Bound
                T, s = bound.concatenate_with(T, s)
            translations.append((T, s))
        return translations

    def get_voronoi_translations(self, source_fn, sink_fn, diameter=None):
        """
        memoized _get_voronoi_translations
        Gets translations of p on the source
            considers every possible face path from source face to sink face
        :param source_fn: face name of source
        :param sink_fn: face name of sink
        :param diameter: cap on length of face path to consider, None if infinite
        :return: list of (T,s) translation matrix and shift such that each Tp+s translates p to sink face
        """
        if (source_fn, sink_fn, diameter) not in self.memoized_face_translations:
            self.memoized_face_translations[(source_fn, sink_fn, diameter)] = self._get_voronoi_translations(source_fn, sink_fn, diameter=diameter)
        return self.memoized_face_translations[(source_fn, sink_fn, diameter)]

    def get_voronoi_points_from_face_paths(self, p, source_fn, sink_fn, diameter=None):
        """
        Gets voronoi points spawned by p on the source
            considers every possible face path from source face to sink face
        :param p: column vector (np array of dimension (self.dimension,1))
        :param source_fn: face name of source
        :param sink_fn: face name of sink
        :param diameter: cap on length of face path to consider, None if infinite
        :return: list of column vector voronoi points
        """
        points = []
        for (T, s) in self.get_voronoi_translations(source_fn, sink_fn, diameter=diameter):
            points.append(T@p + s)
        return points

    def plot_voronoi(self, p, source_fn, sink_fn, diameter, ax):
        """
        creates a voronoi plot for the sink face from p on a souce face
        :param p: column vector (np array of dimension (self.dimension,1))
        :param source_fn: face name of source
        :param sink_fn: face name of sink
        :param diameter: cap on length of face path to consider, None if infinite
        :param ax: plot to plot on (pyplot, or ax object)
        :return: whether we were successful
        """
        vp = self.get_voronoi_points_from_face_paths(p, source_fn, sink_fn, diameter=diameter)
        if len(vp) >= 2:
            if len(vp) < 4:
                large = 69*sum(np.linalg.norm(p) for p in vp)
                vp.append(np.ones(vp[0].shape)*large)
                vp.append(-np.ones(vp[0].shape)*large)
            points = np.concatenate(vp, axis=1)
            points = points.T
            vor = Voronoi(points)
            fig = voronoi_plot_2d(vor, ax=ax, show_points=False, show_vertices=False, line_colors='black',
                                  line_width=1, line_alpha=1)
            return True
        return False

    def get_cut_locus_points(self, fn):
        """
        Gets all cut locus points on fn by considering all possible intersections of arcs
        :param fn: face name
        :return: (column vector, radius) list
        """
        arcs = self.arcs[fn]
        intersections = []
        for (A, _), (B, _) in itertools.combinations(arcs, 2):
            if A.dist == B.dist:
                # equality is fine since we should never update dist
                ps = A.intersects_with(B)
                if ps:
                    intersections += [(p, A.dist) for p in ps]

        return [(p, r) for (p, r) in intersections if self.is_best_seen(p, r, fn)]

    def is_best_seen(self, p, r, fn):
        """
        returns if the distnace r to p is the best seen out of all arcs assigned to fn
        :param p: column vector (np array of dimension (self.dimension,1))
        :param r: radius we are comparing to
        :param fn: face name
        :return: whether r is the best radius to p
        """
        for (A, _) in self.arcs[fn]:
            if A.within_arc(p) and A.dist < r:
                # if we find an arc containing p that is smaller than r
                return False
        # otherwise, this is the best
        return True

    def add_path_end_to_face(self, p, v, fn, point_info=None):
        """
        takes a vector starting at p on face fn, goes towards v
            finds its endpoint and adds it to the correct face with point_info attached

        :param p: column vector (np array of dimension (self.dimension,1))
        :param v: column vector (np array of dimension (self.dimension,1))
        :param fn: face name
        :param point_info: dictionary of info to attach to path end point
        """
        assert self._face_exists_correctly(fn)
        src: Face = self.faces[fn]
        points = src.points_on_path(p, v)
        (sink, _, p) = points[-1]
        self.add_point_to_face(p, sink.name, point_info)

    def faces_to_plot_n_m(self):
        """
        gives plotting information
            n,m are dimensions for the number of plots we need (will be plotted on an n x m grid)
            face map of (i,j)->face on plot (i,j)
        :return: (face map, n,m)
        """
        grabbable = []
        for fn in self.faces:
            face = self.faces[fn]
            if face.dimension == 2:
                grabbable.append(face)
        m = int(np.ceil(np.sqrt(len(grabbable))))
        n = int(np.ceil(len(grabbable)/m))

        def face_map(i, j):
            k = i*m + j
            if k >= len(grabbable):
                return None
            return grabbable[k]

        return face_map, n, m

    def draw_arc(self, ax, A: Arc, arc_info, n=20):
        """
        draws arc on face
        :param ax: axis to plot on (plt or something)
        :param A: Arc to plot
        :param arc_info: dictionary of info to attach to arc
        :param n: how many points to approximate arc with
        """

        x0, y0 = tuple(A.p.flatten())
        thetas = A.low + (np.arange(n)/(n - 1))*(A.high - A.low)
        X, Y = x0 + A.r*np.cos(thetas), y0 + A.r*np.sin(thetas)

        if arc_info is None:
            arc_info = dict()
        color = None
        plot = True
        if 'color' in arc_info:
            color = arc_info['color']
        if 'plot' in arc_info:
            plot = arc_info['plot']

        if plot:
            ax.plot(X, Y, color=color)

    def plot_face_boundaries(self, axs, legend):
        """
        plots faces and points/arcs on faces

        :param axs: axis to plot on
        :param legend: (i,j)-> whether to put a legend on plot (i,j)
        """
        face_map, n, m = self.faces_to_plot_n_m()

        def ploot(i, j):
            if m > 1 and n > 1:
                return axs[i, j]
            if m == 1 and n == 1:
                return axs
            return axs[m*i + j]

        if self.extra_legend == None:
            self.extra_legend = {fn: False for fn in self.faces}
            for fn in self.faces:
                for F in self.faces[fn].double_face_edge:
                    self.extra_legend[fn] = True
                    self.extra_legend[F.name] = True

        for i in range(n):
            for j in range(m):
                face = face_map(i, j)
                if face is not None:
                    ploot(i, j).set_title("FACE " + str(face.name))
                    path = face.get_path_and_faces()
                    for ((p1, p2), (bound, f)) in path:
                        (x, y) = tuple(p1.flatten())
                        (xp, yp) = tuple(p2.flatten())
                        label = str(f.name)
                        if self.extra_legend[face.name]:
                            if bound.name not in self.seen_bounds:
                                self.seen_bounds.append(bound.name)
                            label += ' (id:' + str(self.seen_bounds.index(bound.name)) + ')'

                        ploot(i, j).plot([x, xp], [y, yp], label=label, alpha=.5)

                    for (p, point_info) in self.points[face.name]:
                        x, y = tuple(np.array(p).flatten())
                        color = None
                        s = None
                        plot = True
                        if point_info is None:
                            point_info = dict()
                        if 'color' in point_info:
                            color = point_info['color']
                        if 's' in point_info:
                            s = point_info['s']
                        if 'plot' in point_info:
                            plot = point_info['plot']
                        if plot:
                            ploot(i, j).scatter(x, y, color=color, s=s)

                    for (A, arc_info) in self.arcs[face.name]:
                        self.draw_arc(ploot(i, j), A, arc_info=arc_info)
                    if legend(i, j):
                        ploot(i, j).legend()

    def interactive_vornoi_plot(self, figsize=None, legend=lambda i, j: False, diameter=None, event_key='button_press_event'):
        """
        :param figsize: initial figure size (inches)
        :param legend: (i,j)-> whether to put a legend on plot (i,j)
        :param diameter: longest path of faces to consider when creating paths for vornoi plot
            (None if infinite)
        :param event_key: when to update special point
            'motion_notify_event' or 'button_press_event' are tested
            https://matplotlib.org/stable/users/explain/figure/event_handling.html
        """
        plt.rcParams["figure.autolayout"] = True
        face_map, n, m = self.faces_to_plot_n_m()
        fig, axs = plt.subplots(n, m, figsize=figsize)

        def ploot(i, j):
            if m > 1 and n > 1:
                return axs[i, j]
            if m == 1 and n == 1:
                return axs
            return axs[m*i + j]

        list_axes = []
        for i in range(n):
            for j in range(m):
                list_axes.append((ploot(i, j), (i, j)))

        def ploot_inv(ax):
            for (x, (i, j)) in list_axes:
                if x == ax:
                    return (i, j)
            return None

        def mouse_event(event):
            ax = event.inaxes
            if ax is None:
                return
            for i in range(n):
                for j in range(m):
                    ploot(i, j).cla()
            p = np.array([[event.xdata], [event.ydata]])
            (i, j) = ploot_inv(ax)
            fc = face_map(i, j)
            fc: Face
            if fc is None:
                return
            p = fc.get_closest_point(p)

            self.plot_face_boundaries(axs, legend=legend)
            for i in range(n):
                for j in range(m):
                    ploot(i, j).set_xticks([])
                    ploot(i, j).set_yticks([])
            ax.scatter(p[0, 0], p[1, 0], color='purple')

            source_fn = fc.name

            for i in range(n):
                for j in range(m):
                    face = face_map(i, j)
                    if face is not None:
                        xlim, ylim = ploot(i, j).get_xlim(), ploot(i, j).get_ylim()
                        self.plot_voronoi(p, source_fn, face.name, diameter=diameter, ax=ploot(i, j))
                        ploot(i, j).set_xlim(xlim)
                        ploot(i, j).set_ylim(ylim)
            plt.show()

        cid = fig.canvas.mpl_connect(event_key, mouse_event)
        self.plot_face_boundaries(axs, legend=legend)

        for i in range(n):
            for j in range(m):
                ploot(i, j).set_xticks([])
                ploot(i, j).set_yticks([])
        plt.show()

    def plot_faces(self, save_image=None, show=False, figsize=None, legend=lambda i, j: True, voronoi=None):
        """
        plots all faces of graph
        :param save_image: whether to save the image
        :param show: whether to show the plot
        :param figsize: size of figure (inches)
        :param legend: (i,j)-> whether to put a legend on plot (i,j)
        :param voronoi: list of (p, source face, diameter) points to use in the vornoi plot

        """
        face_map, n, m = self.faces_to_plot_n_m()

        fig, axs = plt.subplots(n, m, figsize=figsize)

        def ploot(i, j):
            if m > 1 and n > 1:
                return axs[i, j]
            if m == 1 and n == 1:
                return axs
            return axs[m*i + j]

        for i in range(n):
            for j in range(m):
                ploot(i, j).set_xticks([])
                ploot(i, j).set_yticks([])
        self.plot_face_boundaries(axs, legend=legend)
        for i in range(n):
            for j in range(m):
                face = face_map(i, j)
                if face is not None:
                    xlim, ylim = ploot(i, j).get_xlim(), ploot(i, j).get_ylim()
                    if voronoi is not None:
                        # ignore_locus_points = self.plot_voronoi(face.name, ploot(i, j))
                        (p, source_fn, diameter) = voronoi
                        ignore_locus_points = self.plot_voronoi(p, source_fn, face.name, diameter=diameter, ax=ploot(i, j))
                        ploot(i, j).set_xlim(xlim)
                        ploot(i, j).set_ylim(ylim)

        if save_image is not None:
            plt.savefig(save_image)
        if show:
            plt.show()
        plt.close()
