import numpy as np, itertools
from matplotlib import pyplot as plt
from scipy.spatial import Voronoi
from src.my_vornoi import voronoi_plot_2d
import fractions


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


def project_p_onto_line(p, a, v):
    """
    gets the projection point of p onto a line starting at a and going towards vector v
    :param p: column vector (np array of dimension (self.dimension,1))
    :param a: column vector (np array of dimension (self.dimension,1))
    :param v: column vector (np array of dimension (self.dimension,1))
    :return: column vector (np array of dimension (self.dimension,1))
    """
    # v_norm=v/np.linalg.norm(v)
    # return a+v_norm*(np.dot(v_norm.T,p-a))
    return a + v*(np.dot(v.T, p - a))/np.square(np.linalg.norm(v))


class Arc:
    def __init__(self, p, low, high, r, tolerance, distance):
        """
        Defines an arc in 2d space
        :param p: point of origin, column vector (np array of dimension (self.dimension,1))
        :param low: low angle in radians
        :param high: high angle in radians (must be >= low)
        :param r: radius, scalar
        :param tolerance: tolerance for intersection and containment
        :param distance: original distance, or used to keep track of arbitrary data (can be different if we 'diffract')
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
        return (Arc(self.p, self.low, theta, self.r, tolerance=self.tol, distance=self.dist),
                Arc(self.p, theta, self.high, self.r, tolerance=self.tol, distance=self.dist))

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
                fraction_along_segment = [(xi - p1x)/dx if abs(dx) > abs(dy) else (yi - p1y)/dy for xi, yi in
                                          intersections]
                intersections = [pt for pt, frac in zip(intersections, fraction_along_segment) if 0 <= frac <= 1]
            if len(intersections) == 2 and abs(
                    discriminant) <= self.tol:  # If line is tangent to circle, return just one point (as both intersections have same location)
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
        return "ARC:{p:" + str(tuple(self.p.flatten())) + "; r:" + str(self.r) + "; dist:" + str(
            self.dist) + "; range:" + str((self.low, self.high)) + "}"


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
        # rescale so that |m|=1
        # dividing both m and b by |m| yields this with the same bound (mx<=b iff mx/|m|<=b/|m|)
        scale = np.linalg.norm(m)
        self.m = m/scale
        self.b = b/scale
        self.s = s
        self.T = T
        self.si = si
        self.dimension = self.check_valid(dimension)
        if name is None:
            base_id = str(tuple(self.m.flatten())) + str(self.b) + str(tuple(self.s.flatten())) + str(
                tuple(self.T.flatten())) + str(tuple(self.si.flatten())) + str(self.dimension)
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
            raise Exception("tried putting in bounds of incorrect dimensions")
        if (np.shape(self.T)[0] != np.shape(self.T)[1] or
                np.shape(self.m)[1] != np.shape(self.s)[0] or
                np.shape(self.s)[0] != np.shape(self.si)[0] or
                np.shape(self.si)[0] != np.shape(self.m)[1] or
                np.shape(self.m)[1] != np.shape(self.T)[0]
        ):
            raise Exception("tried putting in bounds of inconsistent dimensions")
        if dimension is None:
            dimension = np.shape(self.T)[0]
        else:
            if dimension != np.shape(self.T)[0]:
                raise Exception("inconsistent dimensions")
        if not self.within(np.zeros(dimension)):
            raise Exception("zero point needs to be within the face")
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
        return Arc(self.shift_point(A.p), angle, angle + diff, A.r, tolerance=A.tol, distance=A.dist)

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

    def get_plot_bounds(self):
        """
        returns the bounds of the graphical plot of the face
        :return: (xlim, ylim)
        """
        if self.dimension != 2:
            raise Exception("cannot graph a non-2d face")
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

    def line_within_bounds(self, p, q):
        """
        returns if the line p->q is within the face
        :param p: column vector (np array of dimension (self.dimension,1))
        :param q: column vector (np array of dimension (self.dimension,1))
        :return: boolean
        """
        # enough to test whether either of the endpoints are within the bounds or if the line exits the face
        return self.within_bounds(p) or self.within_bounds(q) or (self.get_exit_point(p, q - p) is not None)

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
            raise Exception("line passes outside of face")
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
                    for path in f.face_paths_to(fn, visited_names=visited_names.copy(),
                                                diameter=None if diameter is None else diameter - 1):
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
            raise Exception("face " + str(fn) + " not a neighbor")
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

        arr = [[(B, self)] if F is None else F.push_arc_to_faces(self.bound_of_face(F).shift_arc(B)) for (B, F) in
               zip(tuple(arr), faces)]
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
    def __init__(self, tolerance, faces=None):
        """
        A set of faces, representing a shape
        :param tolerance: tolerance for within_face and such
            adds this to autogenreated faces
        :param faces: initial set of faces
            Note: usually start empty and add faces as we go
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
        self.extra_data = dict()
        if not self.is_polyhedra():
            print("WARNING: SHAPE IS NOT BOUNDARY OF POLYHEDRON, CUT LOCUS MAY NOT BE MATHEMATICALLY SOUND")

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

    def is_polyhedra(self):
        return True

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
            bound_path = []
            for (bound, F) in path:
                bound: Bound
                bound_path.append((bound, F))
                T, s = bound.concatenate_with(T, s)
            translations.append((T, s, bound_path))
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
            self.memoized_face_translations[(source_fn, sink_fn, diameter)] = self._get_voronoi_translations(source_fn,
                                                                                                             sink_fn,
                                                                                                             diameter=diameter)
        return self.memoized_face_translations[(source_fn, sink_fn, diameter)]

    def get_voronoi_points_from_face_paths(self, p, source_fn, sink_fn, diameter=None):
        """
        Gets voronoi points spawned by p on the source
            considers every possible face path from source face to sink face
        :param p: column vector (np array of dimension (self.dimension,1))
        :param source_fn: face name of source
        :param sink_fn: face name of sink
        :param diameter: cap on length of face path to consider, None if infinite
        :return: list of column vector voronoi points, list of bounds that connect source to sink
        """
        points = []
        bound_paths = []
        for (T, s, bound_path) in self.get_voronoi_translations(source_fn, sink_fn, diameter=diameter):
            points.append(T@p + s)
            bound_paths.append(bound_path)
        return points, bound_paths

    def point_within_cell(self, v, segments, p=None):
        """
        checks whether v is within the cell bounded by segments
            takes an interior point and checking whether v is on the same side of each bound as this point
        :param v: column vector (np array of dimension (self.dimension,1))
        :param segments: list of (a,b) line segments making up cell
        :param p: interior point to check, if None, just takes average of vertices in segments
        :return: whether v is within cell
        """
        if p is None:
            p = np.zeros(2)  # point that is in center of C
            for a, b in segments:
                p += a + b
            p = p.reshape((2, 1))/(2*len(segments))

        for (a, b) in segments:
            a = a.reshape((len(a), 1))
            b = b.reshape((len(b), 1))

            outside_dir = project_p_onto_line(p, a, b - a) - p
            # this is a vector pointing 'outside' segment (a,b)
            # uses the fact that (a,b) bounds a cell containing p, and that p*-p points out of the cell
            #   (where p* is p's projection onto line extending (a,b))
            vertex_dir = project_p_onto_line(v, a, b - a) - v
            # this is a vector pointing from vertex v to its projection on line extending (a,b)
            if np.dot(outside_dir.T, vertex_dir) >= -self.tol*np.linalg.norm(vertex_dir)*np.linalg.norm(outside_dir):
                # equivalent to outside_dir*vertex_dir/(|outside_dir||vertex_dir|)>=-tol
                # we multiply since vertex_dir may be 0 (outside_dir should not be 0)
                # if this is nonnegative then they point the same direction and v is inside this bound
                continue
            else:
                return False
        return True

    def check_if_valid(self, p, source, bound_path, segments):
        """
        pick a face F and p have vornonoi cell C,
            this method returns true if all paths p to (C intersect F) are 'correct'
                correct meaning there is actually a line using faces on bound_path that goes from p to q for all q in F
        we use that being 'correct' in this sense is a convex property wrt the location of q
            shown in the paper, but basically, the average q of any 2 'correct' q1 and q2 must have a path pq that lies between pq1 and pq2
                thus, since we have a path of faces, the path pq must also lie in that path
        because of this convexity, we only need to check a finite set of points whose convex hull contains (C intersect F)
            the specific points we need are all vertices of C in F, all vertices of F in C,
                and all intersections of boundary lines of C and F

        :param p: column vector (np array of dimension (self.dimension,1))
            Note that p's coordinates should be with the sink face as the origin,
        :param source: source face, p's original face
            (p is written from the POV of the sink face so p will need to be transformed before it literally is in source face)
        :param bound_path: list of (bound, F) representing the path of bounds from source face (with p on it) to sink face
            Must have at least one element, this doesnt make sense if the source is the sink face
        :param segments: list of (a,b) line segments making up voronoi cell of p
        :return: list of points that are relevant on intersection of C and F
        """
        checking_pts = []  # points to check
        (_, sink) = bound_path[-1]
        sink: Face
        source: Face

        DEBUG = np.linalg.norm(p.flatten() - np.array([5.84601932, 0.59591081])) < .001

        QBUG = np.array([0.78747522, -1.18753154])

        # all vertices of C that are in F
        # also all boundary intersections
        for (a, b) in segments:
            a = a.reshape((len(a), 1))
            b = b.reshape((len(b), 1))

            # check if a or b is within F, and add them to points to check
            for q in (a, b):
                if sink.within_bounds(q) and not any(np.array_equal(q, q1) for q1 in checking_pts):
                    checking_pts.append(q)
            # check if either the line ab or ba exits F, and add them to points to check
            for q in (sink.get_exit_point(a, b - a), sink.get_exit_point(b, a - b)):
                if q is not None and not any(np.array_equal(q, q1) for q1 in checking_pts):
                    checking_pts.append(q)
        # all vertices of F that are in C
        for (v, _) in sink.get_vertices():
            # if self.point_within_cell(v, segments, p=p):
            if self.point_within_cell(v, segments, p=None):
                checking_pts.append(v)

        # now go through and check each point
        for q in checking_pts:
            q_orig = q.copy()
            p_temp = p.copy()
            # this is a little annoying since bound path goes from p to q,
            #   but it is much easier to check in the opposite direction
            for (inv_bound, face) in bound_path[::-1]:
                face: Face
                # since bound goes from p to q, we need to invert it to go the other way
                inv_bound: Bound
                bound = inv_bound.get_inverse_bound()
                if not face.within_bounds(q):
                    # if the end that we check is outside of the face, we fail
                    print(p.flatten(), 'invalid with point ', q_orig.flatten())
                    return False

                # set new q to the point where qp exits the current face
                q_temp = face.get_exit_point(q, p_temp - q)
                # now update p and q for the next face
                if q_temp is None:
                    # EDGE CASE: p is on the same face as q
                    # this is a literal edge case, as p is on the boundary of the face
                    # then we can simply set p and q to the same value and continue to the next step
                    # the next check will make sure the boundary that p sits on is actually the correct boundary
                    q_temp = p_temp
                q = bound.shift_point(q_temp)
                p_temp = bound.shift_point(p_temp)
            # here, we do one last check to see if our last q is actually in the source face
            if not source.within_bounds(q):
                print(p.flatten(), 'invalid with point ', q.flatten())
                return False
        return True

    def filter_out_points(self, points, bound_paths, source, sink, do_filter=True):
        """
        repeatedly makes voronoi complices, looks at relevant points, then filters out points that do not pass through correct faces
        Note: this augments points by placing 4 very distant points that do not affect the relevant section of the complex
            (artifact of voronoi complex implementation)
        :param points: list of column vectors
        :param bound_paths: list of (list of (bound, F) representing the path of bounds from source face (with p on it) to sink face)
        :param source: source face
        :param sink: sink face
        :param do_filter: whether to filter the points
                    probably only set to false when surface is not polyhedra (i.e. torus or mirror)
        :return: list of points and bound paths that are relevant, augmented by four bounding points that are very far away
        """

        def augment_point_paths(pts, bnd_paths):
            """
            adds 4 large points so that the vornoi diagram is always defined
                puts them in corners of extremely large bounding box
            :param pts: array of column vector points (must be populated)
            :param bnd_paths: array of paths (will add 'None' to this)
            :return (points, bound_paths), both sorted by angle of point for non-augmented points, and augmented at end
            """
            large = 69*(sum(np.linalg.norm(p) for p in pts) + 1)
            shape = pts[0].shape
            vs = [np.ones(shape)]
            vs.append(vs[0].copy())
            vs[1][0, 0] = -vs[1][0, 0]
            vs.append(-vs[0])
            vs.append(-vs[1])

            together = list(zip(pts, bnd_paths))
            together.sort(key=lambda x: np.arctan2(x[0][1, 0], x[0][0, 0])%(2*np.pi))
            for i in range(4):
                large_pt = large*vs[i]
                # pts.append(large_pt)
                # bnd_paths.append(None)
                together.append((large_pt, None))
            # together = list(zip(pts, bnd_paths))
            # together.sort(key=lambda x: np.arctan2(x[0][1, 0], x[0][0, 0])%(2*np.pi))
            return [p for (p, _) in together], [pth for (_, pth) in together]

        bad_point_found = True
        while bad_point_found:
            bad_point_found = False
            # start with set of points, augment them
            if len(points) == 0:
                # here, we return, as this means we have run out of points to check
                return None, None, None
            vp, bound_paths = augment_point_paths(points, bound_paths)
            points = np.concatenate(vp, axis=1)
            points = points.T
            # find the cell complex of them
            vor = Voronoi(points)
            fig, point_to_segments = voronoi_plot_2d(vor, ax=None)

            # gather the relevant points: the ones whose cells intersect the sink face
            # to check this, we can split into two cases
            #  (1): a line from the vornoi cell lies within the face
            #  (2): the face lies completely within the voronoi cell
            # there are no other ways for this intersection to happen
            relevant_points = []
            relevant_bound_paths = []
            relevant_cells = []

            duplicates = []
            for p_idx in point_to_segments:
                point = points[(p_idx,), :]  # row vector of point that created this
                point_included = False
                if len(point_to_segments[p_idx]) == 0:
                    # edge case: point does not have a voronoi cell
                    # this happens with duplicate or very close points
                    # need to be careful here
                    # if the duplicate is not a 'bad point' then we are fine, as any of them would have the same voronoi complex
                    # (i.e. now there are just more possible face paths to get the same cell)
                    # if the duplicate is a 'bad point', we need to remove just that point and start over
                    # its possible that the duplicates have different face paths,
                    #   and in this case we need to just remove the one that we confirm is bad
                    duplicates.append((point.T, bound_paths[p_idx], point_to_segments[p_idx]))
                    # we will save all duplicates, and if we fail anything, we will add back all of them
                    continue
                for a, b in point_to_segments[p_idx]:
                    # if any line of the point's voronoi cell is in face F, we call this point relevant and continue
                    if sink.line_within_bounds(a.reshape((2, 1)), b.reshape((2, 1))):
                        relevant_points.append(point.T)
                        relevant_bound_paths.append(bound_paths[p_idx])
                        relevant_cells.append(point_to_segments[p_idx])
                        point_included = True
                        break
                if not point_included:
                    # check edge case:
                    # the sink face F is completely within cell C
                    # if this is true, we also add the point
                    # to check this, we only need to test if an arbitrary vertex of F is in C
                    (v, _) = sink.get_vertices()[0]

                    # if self.point_within_cell(v,point_to_segments[p_idx],p=point.T):
                    if self.point_within_cell(v, point_to_segments[p_idx], p=None):
                        # if this is actually what is happening, we consider this relevant and continue
                        relevant_points.append(point.reshape((2, 1)))
                        relevant_bound_paths.append(bound_paths[p_idx])
                        relevant_cells.append(point_to_segments[p_idx])
            if not relevant_points:
                print("ERROR NO RELEVANT POINTS")
                # this should not happen as we only fully skip repeats
                return None, None, None
            if not do_filter:
                points, bound_paths = augment_point_paths(relevant_points, relevant_bound_paths)
                return points, bound_paths, relevant_cells

            # now iterate through relevant points and see if they actually have the property we want
            # i.e. points on their voronoi cell intersect the sink face are actually in the path of faces we say they are
            # if any fail, we remove it and restart loop
            # if all succeed, we return the relevant points, bound paths, and segments
            for idx, (pt, bound_path, cell_segment) in enumerate(
                    zip(relevant_points, relevant_bound_paths, relevant_cells)):
                if not self.check_if_valid(pt, source, bound_path, cell_segment):
                    bad_point_found = True

                    # remove this point and continue the loop
                    relevant_points.pop(idx)
                    relevant_bound_paths.pop(idx)

                    points = relevant_points + [dup_pt for (dup_pt, _, _) in duplicates]
                    bound_paths = relevant_bound_paths + [dup_pth for (_, dup_pth, _) in duplicates]
                    break
            if bad_point_found:
                continue
            # otherwise, no bad point is found, we can just augment and return here
            points, bound_paths = augment_point_paths(relevant_points, relevant_bound_paths)
            return points, bound_paths, relevant_cells

    def plot_unwrapping(self, p, source_fn, sink_fn, diameter, ax,
                        i_to_display=None,
                        orient_string='',
                        do_filter=True,
                        label_diagram=False):
        """
        plots an unwrapping of the cut locus on sink face from point p on source face
        :param p: column vector (np array of dimension (self.dimension,1))
        :param source_fn: face name of source
        :param sink_fn: face name of sink
        :param diameter: cap on length of face path to consider, None if infinite
        :param ax: plot to plot on (pyplot, or ax object)
        :param i_to_display: which path to display, if we are only showing one
        :param orient_string: string to add to face annotation to show orientation
        :param do_filter: Whether to filter voronoi cell points based on correctness of paths
                should probably always be true, unless we are not looking at polyhedra
        :param label_diagram: whether to label points and lines
        :return: whether there is any more to show (i.e. i_to_display is None or larger than the number of paths)
            returns none if not enough points
        """
        if ax is None:
            ax = plt.gca()
        vp, bound_paths = self.get_voronoi_points_from_face_paths(p, source_fn, sink_fn, diameter=diameter)
        source: Face = self.faces[source_fn]
        sink: Face = self.faces[sink_fn]

        def plot_label_face(ax, face, name, center, rot_v, color, linewidth):
            """
            plots and labels face on ax
            """
            theta = np.arctan2(rot_v[1, 0], rot_v[0, 0])
            for i in range(len(face)):
                v1 = face[i]
                v2 = face[(i + 1)%len(face)]
                ax.plot([v1[0], v2[0]], [v1[1], v2[1]], color=color, linewidth=linewidth)
            ax.annotate(str(name) + orient_string, (center[0], center[1]), rotation=np.degrees(theta))

        if len(vp) >= 2:  # if there is only one, the cut locus does not exist on this face
            relevant_points, relevant_bound_paths, relevant_cells = self.filter_out_points(vp, bound_paths, source,
                                                                                           sink, do_filter=do_filter)

            if relevant_points is None:
                return None
            all_trans_shown = []

            labeled = False
            disp_i = -1
            n = len([pth for pth in relevant_bound_paths if pth is not None])
            if i_to_display is not None and i_to_display >= n:
                i_to_display = None  # here, just show all
            special_face = None
            for pt_idx, (pt, path) in enumerate(zip(relevant_points, relevant_bound_paths)):
                # we can graph pt here
                tracker_points = [np.zeros((2, 1)), coltation(0), coltation(np.pi/2)]  # 0, x, y
                if path is not None:
                    disp_i += 1
                    if (i_to_display is not None) and (disp_i != i_to_display):
                        # skip this if we are skipping, and the path is not the correct path
                        continue
                    face_tracking = [[v.copy() for (v, _) in
                                      source.get_vertices()]]  # tracking the vertices of each face and their eventual location
                    center_tracking = [np.zeros((2, 1))]  # tracking the center of each face
                    rot_tracking = [np.array([[1], [0]])]  # tracks the 0 angle of each face
                    face_name_tracking = [source_fn]  # tracks face names

                    for (bound, F) in path:
                        bound: Bound
                        face_tracking = [[bound.shift_point(v) for v in vees] for vees in face_tracking] + [
                            [v.copy() for (v, _) in F.get_vertices()]]
                        center_tracking = [bound.shift_point(v) for v in center_tracking] + [np.zeros((2, 1))]
                        rot_tracking = [bound.shift_vec(v) for v in rot_tracking] + [np.array([[1], [0]])]
                        face_name_tracking = face_name_tracking + [F.name]
                        tracker_points = [bound.shift_point(track) for track in tracker_points]
                    iteration = list(zip(face_tracking, face_name_tracking, center_tracking, rot_tracking))
                    special_face = iteration[-1]
                    all_trans_shown.append(tracker_points + [pt])

                    for face, name, center, rot_v in iteration[:-1]:
                        plot_label_face(ax=ax, face=face, name=name, center=center, rot_v=rot_v, color='blue',
                                        linewidth=1)
                    label = None
                    if not labeled:
                        label = '$p$ copies'
                        labeled = True
                    ax.scatter(pt[0], pt[1], color='purple', label=label, alpha=1, s=40)
                    if label_diagram:
                        label_pt = pt + [[.1], [.1]]
                        ax.annotate('$p^{(' + str(pt_idx) + ')}$', (label_pt[0], label_pt[1]), rotation=0)
            (face, name, center, rot_v) = special_face
            plot_label_face(ax=ax, face=face, name=name, center=center, rot_v=rot_v, color='red', linewidth=2)
            # ax.scatter([0], [0], label='center', alpha=.5, s=80)
            relevant_points = np.concatenate(relevant_points, axis=1)

            vor = Voronoi(relevant_points.T)
            xlim, ylim = ax.get_xlim(), ax.get_ylim()
            voronoi_plot_2d(vor, ax=ax, show_points=False, show_vertices=False, line_colors='black',
                            line_width=2, line_alpha=1, label_lines=label_diagram)
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            ax.legend()
            return all_trans_shown
        return None

    def plot_voronoi(self, p, source_fn, sink_fn, diameter, ax, do_filter=True):
        """
        creates a voronoi plot for the sink face from p on a souce face
        :param p: column vector (np array of dimension (self.dimension,1))
        :param source_fn: face name of source
        :param sink_fn: face name of sink
        :param diameter: cap on length of face path to consider, None if infinite
        :param ax: plot to plot on (pyplot, or ax object)
        :param do_filter: Whether to filter voronoi cell points based on correctness of paths
                should probably always be true, unless we are not looking at polyhedra
        :return: whether we were successful
        """
        vp, bound_paths = self.get_voronoi_points_from_face_paths(p, source_fn, sink_fn, diameter=diameter)

        if len(vp) >= 2:  # if there is only one point, the cut locus does not exist on this face
            relevant_points, relevant_bound_paths, relevant_cells = self.filter_out_points(vp, bound_paths,
                                                                                           self.faces[source_fn],
                                                                                           self.faces[sink_fn],
                                                                                           do_filter=do_filter)
            if relevant_points is None:
                return False
            points = np.concatenate(relevant_points, axis=1)

            points = points.T
            vor = Voronoi(points)
            fig, point_to_segments = voronoi_plot_2d(vor, ax=ax, show_points=False, show_vertices=False,
                                                     line_colors='black',
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
                ploot(i, j).set_xticks([])
                ploot(i, j).set_yticks([])
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

    def interactive_vornoi_plot(self,
                                figsize=None,
                                legend=lambda i, j: False,
                                diameter=None,
                                event_key='button_press_event',
                                source_fn_p=None,
                                show=True,
                                save=None,
                                do_filter=True,
                                font_size=None,
                                ):
        """
        :param figsize: initial figure size (inches)
        :param legend: (i,j)-> whether to put a legend on plot (i,j)
        :param diameter: longest path of faces to consider when creating paths for vornoi plot
            (None if infinite)
        :param event_key: when to update special point
            'motion_notify_event' or 'button_press_event' are tested
            https://matplotlib.org/stable/users/explain/figure/event_handling.html
        :param source_fn_p: if specified, use this face and point as the source
            (face name, column vector)
        :param show: whether to display plot
        :param save: file name to save initial image to
            (none if not saved)
        :param do_filter: Whether to filter voronoi cell points based on correctness of paths
                should probably always be true, unless we are not looking at polyhedra
        :param font_size: font size to use for plot (default if None)
        """
        plt.rcParams["figure.autolayout"] = True
        face_map, n, m = self.faces_to_plot_n_m()
        fig, axs = plt.subplots(n, m, figsize=figsize)
        if font_size is not None:
            plt.rcParams.update({'font.size': font_size})

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

        def full_v_plot_from_point_axis(p, ax):
            (i, j) = ploot_inv(ax)
            fc = face_map(i, j)
            fc: Face
            if fc is None:
                return
            self.plot_face_boundaries(axs, legend=legend)
            ax.scatter(p[0, 0], p[1, 0], color='purple')

            source_fn = fc.name

            for i in range(n):
                for j in range(m):
                    face = face_map(i, j)
                    if face is not None:
                        xlim, ylim = ploot(i, j).get_xlim(), ploot(i, j).get_ylim()
                        self.plot_voronoi(p, source_fn, face.name, diameter=diameter, ax=ploot(i, j),
                                          do_filter=do_filter)
                        ploot(i, j).set_xlim(xlim)
                        ploot(i, j).set_ylim(ylim)

        def mouse_event(event):
            ax = event.inaxes
            if ax is None:
                return
            p = np.array([[event.xdata], [event.ydata]])
            (i, j) = ploot_inv(ax)
            fc = face_map(i, j)
            fc: Face
            if fc is None:
                return
            for i in range(n):
                for j in range(m):
                    ploot(i, j).cla()
            p = fc.get_closest_point(p)
            full_v_plot_from_point_axis(p, ax)
            plt.show()

        if event_key is not None:
            cid = fig.canvas.mpl_connect(event_key, mouse_event)
            self.plot_face_boundaries(axs, legend=legend)
        else:
            fn, p = source_fn_p
            if fn not in self.faces:
                raise Exception("invalid file name specified: " + str(fn))
            if not self.faces[fn].within_bounds(p):
                temp = str(tuple(p.flatten()))
                p = self.faces[fn].get_closest_point(p)
                print("WARNING: point " + temp + ' not in face, taking closest point: ' + str(tuple(p.flatten())))
            I, J = None, None
            for i in range(n):
                for j in range(m):
                    face = face_map(i, j)
                    if face is not None:
                        if face.name == fn:
                            I, J = i, j
            full_v_plot_from_point_axis(p, ploot(I, J))
        if save is not None:
            plt.savefig(save)
        if show:
            plt.show()

    def interactive_unwrap(self,
                           figsize=None,
                           legend=lambda i, j: False,
                           diameter=None,
                           track=True,
                           single_display=True,
                           source_fn_p=None,
                           sink_fn=None,
                           show=True,
                           save=None,
                           orient_string='',
                           do_filter=True,
                           font_size=None,
                           label_diagram=False,
                           ):
        """
        :param figsize: initial figure size (inches)
        :param legend: (i,j)-> whether to put a legend on plot (i,j)
        :param diameter: longest path of faces to consider when creating paths for vornoi plot
            (None if infinite)
        :param track: whether to track cut locus of cursor on movement
        :param single_display: whether to only display one path at once
        :param source_fn_p: if specified, use this face and point as the source
            (face name, column vector)
        :param sink_fn: if specified, use this face as the sink
        :param show: whether to display plot
        :param save: file name to save initial image to
            (none if not saved)
        :param orient_string: string to add onto face annotation to show orientation
        :param do_filter: Whether to filter voronoi cell points based on correctness of paths
                should probably always be true, unless we are not looking at polyhedra
        :param font_size: font size to use for plot (default if None)
        :param label_diagram: whether to label points and lines

        """
        plt.rcParams["figure.autolayout"] = True
        if font_size is not None:
            plt.rcParams.update({'font.size': font_size})
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

        self.extra_data['unwrap_source_fn'] = None
        self.extra_data['p'] = None
        self.extra_data['unwrap_sink_fn'] = None
        if source_fn_p is not None:
            self.extra_data['unwrap_source_fn'], self.extra_data['p'] = source_fn_p
            if not self.faces[self.extra_data['unwrap_source_fn']].within_bounds(self.extra_data['p']):
                temp = str(tuple(self.extra_data['p'].flatten()))
                self.extra_data['p'] = self.faces[self.extra_data['unwrap_source_fn']].get_closest_point(
                    self.extra_data['p'])
                print("WARNING: point " + temp + ' not in face, taking closest point: ' + str(
                    tuple(self.extra_data['p'].flatten())))
            print('source face:', self.extra_data['unwrap_source_fn'])
            print('p:', self.extra_data['p'].flatten())
        if sink_fn is not None:
            self.extra_data['unwrap_sink_fn'] = sink_fn
        self.extra_data['unwrap_counter'] = 0
        self.extra_data['unwrap_source_plotted'] = False

        def spin():

            if self.extra_data['unwrap_source_fn'] is not None and self.extra_data['unwrap_sink_fn'] is not None:
                # if we have finished both
                plt.clf()
                i_to_display = None
                if single_display:
                    i_to_display = self.extra_data['unwrap_counter']
                all_trans_shown = self.plot_unwrapping(self.extra_data['p'], self.extra_data['unwrap_source_fn'],
                                                       self.extra_data['unwrap_sink_fn'],
                                                       diameter=diameter, ax=plt.gca(), i_to_display=i_to_display,
                                                       orient_string=orient_string, do_filter=do_filter,
                                                       label_diagram=label_diagram)
                print('point locations:')
                for zero, xvec, yvec, p in all_trans_shown:
                    print('p copy:', p.flatten())
                    print('\tshift:', zero.flatten())
                    rot_frac = fractions.Fraction(
                        np.arctan2((xvec - zero)[1], (xvec - zero)[0])[0]/np.pi).limit_denominator(1000)
                    print("\trotation:", rot_frac, 'pi')

                    print("\tx vec:", (xvec - zero).flatten())
                    print("\ty vec:", (yvec - zero).flatten())

                plt.xticks([])
                plt.yticks([])
                if single_display and (all_trans_shown is not None):
                    plt.title("click to advance")
                self.extra_data['unwrap_counter'] += 1
            else:
                if not self.extra_data['unwrap_source_plotted'] and self.extra_data['unwrap_source_fn'] is not None:
                    # if we picked a point and havent yet created a plot
                    for i in range(n):
                        for j in range(m):
                            ploot(i, j).cla()
                    self.plot_face_boundaries(axs, legend=legend)
                    for i in range(n):
                        for j in range(m):
                            face = face_map(i, j)
                            if face is not None:
                                xlim, ylim = ploot(i, j).get_xlim(), ploot(i, j).get_ylim()
                                self.plot_voronoi(self.extra_data['p'], self.extra_data['unwrap_source_fn'], face.name,
                                                  diameter=diameter, ax=ploot(i, j), do_filter=do_filter)
                                ploot(i, j).set_xlim(xlim)
                                ploot(i, j).set_ylim(ylim)
                                if self.extra_data['unwrap_source_fn'] == face.name:
                                    ploot(i, j).scatter(self.extra_data['p'][0, 0], self.extra_data['p'][1, 0],
                                                        color='purple')
                    self.extra_data['unwrap_source_plotted'] = True
                if not self.extra_data['unwrap_source_plotted']:
                    plt.suptitle("Click $p$")
                elif self.extra_data['unwrap_sink_fn'] is None:
                    plt.suptitle("Now click sink face")

        def clicked_mouse(event):
            ax = event.inaxes
            p = np.array([[event.xdata], [event.ydata]])
            ij = ploot_inv(ax) if ax is not None else None
            fc = face_map(ij[0], ij[1]) if ij is not None else None
            if self.extra_data['unwrap_source_fn'] is None:
                if fc is None: return
                self.extra_data['unwrap_source_fn'] = fc.name
                self.extra_data['p'] = fc.get_closest_point(p)
                print('source face:', self.extra_data['unwrap_source_fn'])
                print('p:', self.extra_data['p'].flatten())
            elif (self.extra_data['unwrap_sink_fn'] is None and
                  fc is not None and
                  fc.name is not self.extra_data['unwrap_source_fn']):
                self.extra_data['unwrap_sink_fn'] = fc.name
            spin()
            plt.show()

        def moved_mouse(event):
            if self.extra_data['unwrap_source_fn'] is not None:
                return
            ax = event.inaxes
            p = np.array([[event.xdata], [event.ydata]])
            ij = ploot_inv(ax) if ax is not None else None
            fc = face_map(ij[0], ij[1]) if ij is not None else None
            if fc is None: return
            p = fc.get_closest_point(p)
            for i in range(n):
                for j in range(m):
                    ploot(i, j).cla()
            self.plot_face_boundaries(axs, legend=legend)
            ax.scatter(p[0, 0], p[1, 0], color='purple', alpha=.5)

            source_fn = fc.name

            for i in range(n):
                for j in range(m):
                    face = face_map(i, j)
                    if face is not None:
                        xlim, ylim = ploot(i, j).get_xlim(), ploot(i, j).get_ylim()
                        self.plot_voronoi(p, source_fn, face.name, diameter=diameter, ax=ploot(i, j),
                                          do_filter=do_filter)
                        ploot(i, j).set_xlim(xlim)
                        ploot(i, j).set_ylim(ylim)
            plt.show()

        cid = fig.canvas.mpl_connect('button_press_event', clicked_mouse)
        if track:
            cid2 = fig.canvas.mpl_connect('motion_notify_event', moved_mouse)
        self.plot_face_boundaries(axs, legend=legend)
        plt.suptitle("Click $p$")
        spin()
        if save is not None:
            plt.savefig(save)
        if show:
            plt.show()

    def plot_faces(self, save_image=None, show=False, figsize=None, legend=lambda i, j: True, voronoi=None,
                   do_filter=True):
        """
        plots all faces of graph
        :param save_image: whether to save the image
        :param show: whether to show the plot
        :param figsize: size of figure (inches)
        :param legend: (i,j)-> whether to put a legend on plot (i,j)
        :param voronoi: list of (p, source face, diameter) points to use in the vornoi plot
        :param do_filter: Whether to filter voronoi cell points based on correctness of paths
                should probably always be true, unless we are not looking at polyhedra
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
                        ignore_locus_points = self.plot_voronoi(p, source_fn, face.name, diameter=diameter,
                                                                ax=ploot(i, j), do_filter=do_filter)
                        ploot(i, j).set_xlim(xlim)
                        ploot(i, j).set_ylim(ylim)

        if save_image is not None:
            plt.savefig(save_image)
        if show:
            plt.show()
        plt.close()
