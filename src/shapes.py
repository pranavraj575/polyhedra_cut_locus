import numpy as np

from src.my_vornoi import voronoi_diagram_calc
from src.bound import Bound
from src.face import Face


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

    def add_face(self, face=None,face_name=None):
        """
        adds new face to shape
        :param face: Face, if already defined
            if None, initializes new face and adds it
        """
        if face is None:
            if face_name is None:
                face_name=self._pick_new_face_name()
            face = Face(face_name, tolerance=self.tol)
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
                    # print(p.flatten(), 'invalid with point ', q_orig.flatten())
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
                # print(p.flatten(), 'invalid with point ', q.flatten())
                return False
        return True

    def filter_out_points(self, points, bound_paths, source, sink, do_filter=True, ignore_points_on_locus=False):
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
        :param ignore_points_on_locus: whether to ignore single points on cut locus
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
            point_to_segments = dict()
            for point_pair, (seg_type, (a, b)) in voronoi_diagram_calc(points=points).items():
                for pt in point_pair:
                    if pt not in point_to_segments:
                        point_to_segments[pt] = list()
                    if seg_type == 'segment':
                        point_to_segments[pt].append((a, b))
                    elif seg_type == 'ray':
                        point_to_segments[pt].append((a, a + b))
                    else:
                        raise NotImplementedError

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
                    if sink.line_within_bounds(a.reshape((2, 1)),
                                               b.reshape((2, 1)),
                                               ignore_points=ignore_points_on_locus,
                                               ):
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
