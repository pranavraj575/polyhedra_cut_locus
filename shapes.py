import numpy as np, itertools
from matplotlib import pyplot as plt
from collections import defaultdict

TOL = .0001


# rotation matrices
def rotation_T(theta):
    # T where T @ v is v rotatied by theta
    return np.array([[np.cos(theta), -np.sin(theta)],
                     [np.sin(theta), np.cos(theta)]])


def rotation_of_matrix(T):
    # inverse of rotation_T
    cos = T[0][0]
    sin = T[1][0]
    return np.arctan2(sin, cos)


def coltation(theta):
    return rotation_T(theta)[:, [0]]


def rowtation(theta):
    return coltation(theta).T


def flatten(L):
    # flattens a list of lists
    out = []
    for l in L:
        out += l
    return out


class Arc:
    # defines an arc in 2d
    def __init__(self, p, low, high, r, tolerance=TOL):
        assert p.shape == (2, 1)
        assert low <= high

        self.tol = tolerance
        self.p = p
        self.r = r
        diff = high - low
        self.low = low%(2*np.pi)
        self.high = self.low + diff

    def avg_v(self):
        theta = (self.low + self.high)/2
        return coltation(theta)

    def is_empty(self):
        return (self.high - self.low) <= self.tol

    def angle_of(self, q):
        x, y = tuple((q - self.p).flatten())
        theta = np.arctan2(y, x)%(np.pi*2)
        return theta

    def _strict_contains_angle(self, theta):
        if (self.low < theta + 2*np.pi and theta + 2*np.pi < self.high):
            return True, theta + 2*np.pi
        if (self.low < theta and theta < self.high):
            return True, theta
        return False, np.nan

    def _weak_contains_angle(self, theta):
        # assumes theta is on [0,2pi)
        if (self.low <= theta + 2*np.pi + self.tol and theta + 2*np.pi <= self.high + self.tol):
            return True, theta + 2*np.pi
        if (self.low <= theta + self.tol and theta <= self.high + self.tol):
            return True, theta
        return False, np.nan

    def contains_angle(self, theta, interior=False):
        # assumes theta is on [0,2pi)
        return self._strict_contains_angle(theta) if interior else self._weak_contains_angle(theta)

    def within_arc(self, q, interior=False):
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
        # returns two arcs that are this arc split in parts based on p
        if not self.within_arc(q, interior=True):
            raise Exception("split arc called on point not in interior of arc")
        theta = self.angle_of(q)
        if self.low > theta:
            theta = theta + 2*np.pi
        return (Arc(self.p, self.low, theta, self.r, self.tol), Arc(self.p, theta, self.high, self.r, self.tol))

    def _break_arc(self, a, b):
        # breaks arc into pieces based on line segment a-b
        # idea is if arc interacts with a-b, split it into pieces that fully go through a-b and pieces that do not
        # returns list of arcs
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
                out.append(Arc(self.p, th, ph, self.r, self.tol))
            return out
        else:
            return [self]

    def break_arc(self, a, b):
        return [A for A in self._break_arc(a, b) if not A.is_empty()]

    def breakable(self, a, b):
        return len(self.break_arc(a, b)) == 1

    def circle_line_segment_intersection(self, pt1, pt2):
        """ Find the points at which a circle intersects a line-segment.
        This can happen at 0, 1, or 2 points.
        returns list of numpy arrays (col vectors)
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
        A: Arc
        points = []
        p = self.p
        r0 = self.r
        x0,y0=tuple(p.flatten())
        q = A.p
        r1 = A.r
        x1,y1=tuple(q.flatten())
        d = np.linalg.norm(p - q)
        if d<=self.tol: # starting from same point
            return []
        if d<=abs(r0-r1): # one circle is within the other
            return []
        if d>r0+r1: # non intersecting
            return []
        if abs(r0+r1-d) < self.tol:  # if we are exactly here, just take the midpoint
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

            points= [np.array([[x3], [y3]]), np.array([[x4],[y4]])]
        return [pt for pt in points if self.within_arc(pt) and A.within_arc(pt)]

    def __str__(self):
        return "ARC:{p:" + str(tuple(self.p.flatten())) + "; r:" + str(self.r) + "; range:" + str((self.low, self.high)) + "}"


class Face:
    # creates a face
    # note: 0 should be inside the face
    def __init__(self, name, bounds=None, tolerance=.001):
        self.name = name
        self.bounds = []
        self.tol = tolerance
        self.vertices = None
        self.dimension = None
        self.bound_M = None
        self.bound_b = None

        if bounds:
            for (m, b, F, (s, T, si)) in bounds:
                self.add_boundary(m, b, F, (s, T, si))

    def __str__(self):
        return '(Face ' + str(self.name) + ': connected to faces ' + str([F.name for (_, _, F, _) in self.bounds]) + ')'

    def add_boundary(self, m, b, F, shift):
        (s, T, si) = shift
        # adds linear boundary of the form mx<=b
        # note that x is a column vector, m is a row vector, b is a real
        # F represents the neighboring face that this bound touches
        # s,T,si represent the linear transformation that puts a point in our face to the neighboring face
        # represented as a shift, rotation matrix, and another shift
        # s and si are column vectors and T is a rotation matrix
        # s + T x + si = x' where x' is the neighbor face coordinates
        # this is easy since the inverse is -si,T^(-1),-s
        check = True
        if check:
            if (len(np.shape(m)) != 2 or
                    len(np.shape(s)) != 2 or
                    len(np.shape(si)) != 2 or
                    len(np.shape(T)) != 2 or
                    np.shape(m)[0] != 1 or
                    np.shape(s)[1] != 1 or
                    np.shape(si)[1] != 1 or
                    np.shape(T)[0] != np.shape(T)[1] or
                    np.shape(m)[1] != np.shape(s)[0] or
                    np.shape(s)[0] != np.shape(si)[0] or
                    np.shape(si)[0] != np.shape(m)[1] or
                    np.shape(m)[1] != np.shape(T)[0]
            ):
                raise Exception("ERROR: tried putting in bounds of incorrect dimensions")
        if self.dimension is None:
            self.dimension = np.shape(T)[0]
        else:
            if self.dimension != np.shape(T)[0]:
                raise Exception("ERROR: inconsistent dimensions")
        bound = (m, b, F, shift)

        if not self.within_bound(np.zeros(self.dimension), bound):
            raise Exception("ERROR: zero point needs to be within the face")

        self.bounds.append(bound)
        self._create_bound_arrays()

    def grab_intersection(self, p, v, bound):
        (m, b, _, _) = bound
        # find where p+vt intersects the bound (mx<=b)
        # m(p+vt)=b
        # simplify to mvt=b-mp
        # then t=(b-mp)/(mv)
        t = (b - np.dot(m, p))/np.dot(m, v)
        # note that v cannot be parallel to m, which makes sense for this to even work
        return p + v*t

    def _order_vertices(self):
        # orders the vertices by angle (so that a path through all of them is the bounds of the face)
        # only valid if 2d
        if self.dimension != 2:
            return
        self.vertices.sort(key=lambda v:(2*np.pi + np.arctan2(v[0][1][0], v[0][0][0]))%(2*np.pi))

    def get_plot_bounds(self, update=True):
        # returns the bounds of the graphical plot of the face
        # only valid if 2d
        if self.dimension != 2:
            raise Exception("ERROR: cannot graph a non-2d face")
        xmin, xmax = np.inf, -np.inf
        ymin, ymax = np.inf, -np.inf
        for (v, _) in self.get_vertices(update=update):
            x, y = v[:, 0]
            xmin = min(x, xmin)
            xmax = max(x, xmax)
            ymin = min(y, ymin)
            ymax = max(y, ymax)
        return ((xmin, xmax), (ymin, ymax))

    def get_path_and_faces(self, update=True):
        # grabs path of vertices (v0,v1),(v1,v2),...,(vn,v0)
        # takes order 'around' the face
        if not self.dimension == 2:
            raise Exception("this will not work in !=2 dimensions")
        # path=[]
        # for (v, rows) in self.get_vertices(update=update):
        #    path.append(((v[0, 0], v[1, 0]), rows))
        path = self.get_vertices(update=update)
        out = []
        for i in range(len(path)):
            v1, rows = path[i]
            v2, rowsp = path[(i + 1)%len(path)]
            row = None
            for r in rows:
                if r in rowsp:
                    row = r
            out.append(((v1, v2), self.bounds[row][2]))
        return out

    def within_bound(self, p, bound):
        (m, b, _, _) = bound
        return np.dot(m, p) <= b + self.tol

    def within_bounds(self, p):
        # returns if p is in the face
        if self.bound_M is None:
            self._create_bound_arrays()
        return all(self.bound_M@p <= self.bound_b + self.tol)

    def shift_point_with_bound(self, x, bound):
        # shifts point x to equivalent x' on face according to bound
        (m, b, F, (s, T, si)) = bound
        return T@(x + s) + si

    def bound_of_face(self, F):
        for bound in self.bounds:
            if bound[2].name == F.name:
                return bound
        raise Exception("BOUND NOT FOUND WITH THIS FACE")

    def shift_vec_with_bound(self, v, bound):
        # shifts vector v to equivalent v'  according to bound
        (m, b, F, (s, T, si)) = bound
        return T@v

    def shift_angle_with_bound(self, theta, bound):
        (_, _, _, (_, T, _)) = bound
        return theta + (rotation_of_matrix(T))%(2*np.pi)

    def shift_arc_with_bound(self, A, bound):
        diff = A.high - A.low
        angle = self.shift_angle_with_bound(A.low, bound)
        return Arc(self.shift_point_with_bound(A.p, bound), angle, angle + diff, A.r, A.tol)

    def _create_bound_arrays(self):
        self.bound_M = np.array([m[0] for (m, _, _, _) in self.bounds])
        self.bound_b = np.array([[b] for (_, b, _, _) in self.bounds])

    def _create_vertices(self):
        # creates all vertices of the face
        self.vertices = []
        if self.bound_M is None:
            self._create_bound_arrays()
        n = len(self.bounds)
        if n < self.dimension:
            print("WARNING: no vertices since not enough boundaries")
            return
        for rows in itertools.combinations(range(n), self.dimension):
            sub_M = self.bound_M[rows, :]
            sub_b = self.bound_b[rows, :]
            if abs(np.linalg.det(sub_M)) > self.tol:
                vertex = np.linalg.inv(sub_M)@sub_b
                if self.within_bounds(vertex):
                    self.vertices.append((vertex, rows))
        self._order_vertices()

    def get_vertices(self, update=True):
        # grabs all vertices of the face
        if update or self.vertices is None:
            self._create_vertices()
        return self.vertices

    def points_on_path(self, p, v):
        # returns list of (face, initial point, end point) of a ray starting from p, taking vector v
        # if not self.within_bounds(p):
        #    raise Exception("Initial point not in face")
        q = p + v
        # q is the end point on current face
        q_ = p + v
        # q_ represents the point on the line pq that is intersecting the closest boundary

        if self.within_bounds(q):
            return [(self, p, q)]
        closest_bound = None
        for bound in self.bounds:
            if self.within_bound(p, bound) and not self.within_bound(q_, bound):
                # goes from inside bound to outside bound
                q_ = self.grab_intersection(p, v, bound)
                closest_bound = bound
        if closest_bound is None:
            raise Exception("ERROR: line passes outside of face")
        # now we move to the new face with starting point q_ and vector (q-q_)
        # however, we need to translate this into new face coords
        (_, _, F, _) = closest_bound
        new_p = self.shift_point_with_bound(q_, closest_bound)
        new_v = self.shift_vec_with_bound(q - q_, closest_bound)
        rest_of_path = F.points_on_path(new_p, new_v)
        return [(self, p, q_)] + rest_of_path

    def add_boundary_paired(self, f2, m1, b1, s, T, si):
        # adds boundary to face f2, and corresponding bound to self
        self.add_boundary(m1, b1, f2, (s, T, si))
        Ti = np.linalg.inv(T)
        m2 = -m1@Ti
        b2 = -b1 - np.dot(m1, s) - np.dot(m1@Ti, si)
        b2 = b2.flatten()[0]
        f2.add_boundary(m2, b2, self, (-si, Ti, -s))

    def neighbors(self):
        return [f for (_, _, f, _) in self.bounds]

    def get_vertices_of_face_bound(self, fn, update=True):
        # returns the vertices that face fn is touching
        row = None
        for i, (_, _, F, _) in enumerate(self.bounds):
            if F.name == fn:
                row = i
        if row is None:
            raise Exception("ERROR: face " + str(fn) + " not a neighbor")
        out = []
        for v, rows in self.get_vertices(update=update):
            if row in rows:
                out.append(v)
        return out

    def push_arc_to_faces(self, A: Arc):
        # pushes A to the faces it belongs in
        # returns list of [Face, Arc]
        arr = np.array([A])
        for ((a, b), f) in self.get_path_and_faces():
            arr = [B.break_arc(a, b) for B in arr]
            arr = flatten(arr)
        faces = (self.get_correct_next_face(B) for B in arr)

        arr = [[(B, self)] if F is None else F.push_arc_to_faces(self.shift_arc_with_bound(B, self.bound_of_face(F))) for (B, F) in zip(tuple(arr), faces)]
        out = []
        for t in arr:
            out += t
        return out

    def get_correct_next_face(self, A: Arc):
        v = A.avg_v()*A.r
        arr = self.points_on_path(A.p, v)
        if len(arr) == 1:
            return None
        else:
            return arr[1][0]


class Shape:
    def __init__(self, faces=None):
        if faces is None:
            faces = dict()
        self.faces = {face.name:face for face in faces}
        self.points = {face.name:[] for face in self.faces}
        self.arcs = {face.name:[] for face in self.faces}

    def pick_name(self):
        i = len(self.faces)
        while i in self.faces:
            i += 1
        return i

    def _assert_existance(self, fn):
        # checks if face fn exists
        if not fn in self.faces:
            raise Exception("face name does not exist: " + str(fn))
        if not fn in self.points:
            raise Exception("face name does not exist in array points: " + str(fn))

    def add_face(self, face=None):
        if face is None:
            face = Face(self.pick_name())
        self.faces[face.name] = face
        self.reset_face(face.name)

    def reset_face(self, fn):
        self.points[fn] = []
        self.arcs[fn] = []

    def reset_all_faces(self):
        for fn in self.faces:
            self.reset_face(fn)

    def add_point_to_face(self, point, fn):
        self._assert_existance(fn)
        self.points[fn].append(point)

    def add_points_to_face(self, points, fn):
        for point in points:
            self.add_point_to_face(point, fn)

    def add_arc_to_face(self, A: Arc, fn):
        self._assert_existance(fn)
        self.arcs[fn].append(A)

    def add_arc_end_to_face(self, A, fn, arc_info=None):
        # takes an arc
        # finds its endpoints and adds them to the correct face
        self._assert_existance(fn)
        src: Face = self.faces[fn]
        arcs = src.push_arc_to_faces(A)
        for (B, F) in arcs:
            self.arcs[F.name].append((B, arc_info))

    def add_all_cut_locus_points(self, point_info=None,conditional_point_info=None):
        # adds cut locus points to all faces
        # works best if there are a ton of arcs of all different lengths
        for fn in self.faces:
            for (p,r) in self.get_cut_locus_points(fn):
                new_point_info={k:point_info[k] for k in point_info}
                if conditional_point_info is not None:
                    new_point_info.update(conditional_point_info(r))
                self.add_point_to_face((p,new_point_info),fn)

    def get_cut_locus_points(self, fn):
        arcs = self.arcs[fn]
        intersections = []
        for (A,_), (B,_) in itertools.combinations(arcs, 2):
            if A.r==B.r:
                # equality is fine since we never update r
                ps = A.intersects_with(B)
                if ps:
                    intersections += [(p,A.r) for p in ps]

        return [(p,r) for (p,r) in intersections if self.is_best_seen(p,r,fn)]
    def is_best_seen(self,p,r,fn):
        # returns if the distnace r to p is the best seen out of arcs assigned to fn
        for (A,_) in self.arcs[fn]:
            if A.within_arc(p) and A.r<r:
                # if we find an arc containing p that is smaller than r
                return False
        # otherwise, this is the best
        return True
    def add_path_end_to_face(self, p, v, fn, point_info=None):
        # takes a vector starting at p on face fn, goes towards v
        # finds its endpoint and adds it to the correct face with pointinfo attached
        self._assert_existance(fn)
        src: Face = self.faces[fn]
        points = src.points_on_path(p, v)
        (sink, _, p) = points[-1]
        self.add_point_to_face((p, point_info), sink.name)

    def faces_to_plot_n_m(self):
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

    def draw_arc(self, to_plot, A: Arc, arc_info, n=20):
        # draws arc on face
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
            to_plot.plot(X, Y, color=color)

    def plot_faces(self,save_image=None,show=False,figsize=None):

        face_map, n, m = self.faces_to_plot_n_m()

        fig, axs = plt.subplots(n, m,figsize=figsize)
        def ploot(i,j):
            if m>1 and n>1:
                return axs[i, j]
            if m==1 and n==1:
                return axs
            return axs[m*i + j]
        for i in range(n):
            for j in range(m):

                ploot(i, j).set_xticks([])
                ploot(i, j).set_yticks([])

        for i in range(n):
            for j in range(m):
                face = face_map(i, j)
                if face is not None:
                    ploot(i, j).set_title("FACE " + str(face.name))

                    path = face.get_path_and_faces()
                    for ((p1, p2), f) in path:
                        (x, y) = tuple(p1.flatten())
                        (xp, yp) = tuple(p2.flatten())
                        ploot(i, j).plot([x, xp], [y, yp], label=f.name,alpha=.5)
                        # ploot(i, j).annotate(str(f.name),((x+ xp)/2,(y+yp)/2))

                    for (p, point_info) in self.points[face.name]:
                        x, y = tuple(np.array(p).flatten())
                        color = None
                        s = None
                        plot = True
                        if point_info is None:
                            point_info=dict()
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
                    ploot(i, j).legend()

        if save_image is not None:
            plt.savefig(save_image)
        if show:
            plt.show()
        plt.close()