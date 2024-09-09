import numpy as np, itertools
from src.my_vornoi import voronoi_plot_2d
from src.bound import Bound


# rotation matrices
def rotation_T(theta):
    """
    :param theta: angle
    :return: 2x2 rotation matrix T by angle theta (np array)
    """
    # T where T @ v is v rotatied by theta
    return np.array([[np.cos(theta), -np.sin(theta)],
                     [np.sin(theta), np.cos(theta)]])


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

    def line_within_bounds(self, p, q):
        """
        returns if any part of the line p->q is within the face
        :param p: column vector (np array of dimension (self.dimension,1))
        :param q: column vector (np array of dimension (self.dimension,1))
        :return: boolean
        """
        # enough to test whether either of the endpoints are within the bounds or if the line exits the face
        return self.within_bounds(p) or self.within_bounds(q) or (self.get_exit_point(p, q - p) is not None)

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
        return all(fn in dic for dic in (self.faces, self.points))

    def is_polyhedra(self):
        return False

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
            fig, point_to_segments = voronoi_plot_2d(points, ax=None)

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
