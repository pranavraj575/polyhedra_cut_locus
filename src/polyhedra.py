from matplotlib import pyplot as plt
import os
import fractions
import sympy as sym
import numpy as np

from src.utils import coltation, get_correct_end_points
from src.shapes import Shape
from src.face import Face
from src.bound import Bound
from src.my_vornoi import voronoi_diagram_calc


class ConvexPolyhderon(Shape):
    """
    Geodesics on the surface of convex polyhedra
    """

    def __init__(self, tolerance):
        super().__init__(tolerance)

    def is_polyhedra(self):
        return True

    def get_voronoi_diagram(self,
                            p,
                            source_fn,
                            sink_fn,
                            diameter,
                            do_filter=True,
                            intersect_with_face=True,
                            ignore_points_on_locus=False,
                            ):
        """
        implementaiton of algorithm 3

        considers a point p and a perticular sink face, finds the cut locus on the sink face
        returns voronoi diagram (set of lines), as well as relevant (points, face bounds, and faces)
        """
        # TODO: use this for everything
        vp, bound_paths = self.get_voronoi_points_from_face_paths(p, source_fn, sink_fn, diameter=diameter)

        if len(vp) >= 2:  # if there is only one point, the cut locus does not exist on this face
            relevant_points, relevant_bound_paths, relevant_cells = self.filter_out_points(
                vp,
                bound_paths,
                self.faces[source_fn],
                self.faces[sink_fn],
                do_filter=do_filter,
                ignore_points_on_locus=ignore_points_on_locus,
            )
            if relevant_points is None:
                return None
            points = np.concatenate(relevant_points, axis=1)

            points = points.T
            if intersect_with_face:
                face = self.faces[sink_fn]
            else:
                face = None
            point_pair_to_segment = voronoi_diagram_calc(points=points, face=face)
            return point_pair_to_segment, (relevant_points, relevant_bound_paths, relevant_cells)
        return None

    def plot_voronoi_star_unfolding(self,
                                    p,
                                    source_fn,
                                    ax=None,
                                    do_filter=True,
                                    diameter=None,
                                    ignore_points_on_locus=False,
                                    ):
        """
        unfold fixing the source face
        must check the voronoi plot on every face to do this

        """
        if ax is None:
            ax = plt.gca()

        all_stuff = []
        for sink_fn in self.faces:
            if sink_fn != source_fn:
                voronoi_diagram = self.get_voronoi_diagram(p=p,
                                                           source_fn=source_fn,
                                                           sink_fn=sink_fn,
                                                           diameter=diameter,
                                                           do_filter=do_filter,
                                                           intersect_with_face=True,
                                                           ignore_points_on_locus=ignore_points_on_locus,
                                                           )
                if voronoi_diagram is not None:
                    point_pair_to_segment, (relevant_points, _, _) = voronoi_diagram
                    point_pair_to_segment = {point_pair: (a.reshape((2, 1)), b.reshape((2, 1)))
                                             for point_pair, (_, (a, b)) in point_pair_to_segment.items()
                                             }
                    relevant_points = np.concatenate(relevant_points, axis=1)
                    all_stuff.append((point_pair_to_segment, relevant_points))
        contributing_points = [p.flatten()]
        for sink_fn in self.faces:
            if sink_fn != source_fn:
                voronoi_diagram = self.get_voronoi_diagram(p=p,
                                                           source_fn=source_fn,
                                                           sink_fn=sink_fn,
                                                           diameter=diameter,
                                                           do_filter=do_filter,
                                                           intersect_with_face=True,
                                                           ignore_points_on_locus=ignore_points_on_locus,
                                                           )
                if voronoi_diagram is not None:
                    (point_pair_to_segment,
                     (relevant_points, relevant_bound_paths, relevant_cells)
                     ) = voronoi_diagram
                    # keep track of these to figure out which one is our original one
                    relevant_points = np.concatenate(relevant_points, axis=1)

                    for path in relevant_bound_paths:
                        temp_points = relevant_points.copy()

                        # TODO: maybe track all possible faces and check if any points line up?
                        #  if any do, add the thingy
                        temp_stuff = [({p_p: (a.copy(), b.copy()) for p_p, (a, b) in ppts.items()},
                                       rp.copy())
                                      for (ppts, rp) in all_stuff]

                        possible_stuff = []
                        if path is not None:
                            segments = dict()
                            for point_pair in point_pair_to_segment:
                                (_, (a, b)) = point_pair_to_segment[point_pair]
                                segments[point_pair] = (a.reshape((2, 1)), b.reshape((2, 1)))
                            T, s = np.identity(2), np.zeros((2, 1))

                            # list of (v,2) arrays where each column is a vertex
                            vertex_cycles = []
                            rotations = []
                            centers = []
                            names = []

                            for bnd, F in path[::-1]:  # start with sink face
                                bnd: Bound
                                bnd = bnd.get_inverse_bound()
                                T, s = bnd.concatenate_with(T, s)
                                num_v = len(F.vertices)
                                vertex_cycles.append(np.concatenate(
                                    [F.vertices[i%num_v][0]
                                     for i in range(1 + num_v)]
                                    , axis=1))
                                rotations.append(np.identity(2))  # x vec, y vec
                                centers.append(F.basepoint.copy())
                                names.append(str(F.name))

                                # shift point works on arrays as well
                                vertex_cycles = [bnd.shift_point(vc) for vc in vertex_cycles]
                                rotations = [bnd.shift_vec(rt) for rt in rotations]
                                centers = [bnd.shift_point(c) for c in centers]
                                temp_points = bnd.shift_point(temp_points)
                                segments = {point_pair: (bnd.shift_point(a), bnd.shift_point(b))
                                            for point_pair, (a, b) in segments.items()}
                            distances = np.linalg.norm(temp_points - p, axis=0).flatten()
                            for point in temp_points.T:

                                if all([np.linalg.norm(cont_pt - point) > self.tol for cont_pt in contributing_points]):
                                    contributing_points.append(point)
                            best_indices = np.where(distances <= self.tol)[0]
                            plot_faces = False
                            for point_pair, (a, b) in segments.items():
                                if any([best_idx in point_pair for best_idx in best_indices]):
                                    plot_faces = True

                                    ax.plot((a[0], b[0]), (a[1], b[1]), color='black', alpha=1, lw=3, zorder=9)
                                    # ax.scatter((a[0], b[0]), (a[1], b[1]), color='black')
                            if plot_faces:
                                for vertices_cycle, rot, c, name in zip(vertex_cycles, rotations, centers, names):
                                    ax.plot(vertices_cycle[0], vertices_cycle[1], color='blue', alpha=1, lw=1)
                                    continue
                                    self._plot_label_face(ax=ax,
                                                          face=None,
                                                          name=name,
                                                          center=c,
                                                          rot_v=rot[:, (0,)],  # x vector after rotation
                                                          color='blue',
                                                          linewidth=1,
                                                          plot_face=False,
                                                          )
        source = self.faces[source_fn]
        num_v = len(source.vertices)
        vertices_cycle = [source.vertices[i%num_v][0]
                          for i in range(1 + num_v)]
        vertices_cycle = np.concatenate(vertices_cycle, axis=1)
        ax.plot(vertices_cycle[0], vertices_cycle[1],
                color='red')
        ax.scatter(p[0], p[1], color='purple', s=40, zorder=10)  # TODO: mess with zorder

        xlim, ylim = ax.get_xlim(), ax.get_ylim()

        contributing_points = np.stack(contributing_points[1:])
        # ax.scatter(contributing_points[:, 0],
        #           contributing_points[:, 1],
        #           s=10,
        #           color='purple',
        #           )

        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

    def _plot_label_face(self,
                         ax,
                         face,
                         name,
                         center,
                         rot_v,
                         color,
                         linewidth,
                         orient_string='',
                         plot_face=True,
                         zorder=None,
                         label_zorder=None,
                         ):
        """
        plots and labels face on ax
        """
        if label_zorder is None:
            label_zorder = zorder
        theta = np.arctan2(rot_v[1, 0], rot_v[0, 0])
        if plot_face:
            for i in range(len(face)):
                v1 = face[i]
                v2 = face[(i + 1)%len(face)]
                ax.plot([v1[0], v2[0]], [v1[1], v2[1]], color=color, linewidth=linewidth, zorder=zorder)
        ax.annotate(str(name) + orient_string, (center[0], center[1]), rotation=np.degrees(theta), color=color,
                    zorder=label_zorder)

    def plot_unfolding_from_sink(self,
                                 p,
                                 source_fn,
                                 sink_fn,
                                 diameter=None,
                                 ax=None,
                                 i_to_display=None,
                                 orient_string='',
                                 do_filter=True,
                                 label_diagram=False,
                                 p_label_shift=(0., 0.),
                                 line_label_dist=.3,
                                 point_names=None,
                                 ignore_points_on_locus=False,
                                 ):
        # TODO: maybe do the same thing as above method, calculate cut locus for all faces, paste them together
        """
        finds all points that contribute to the cut locus an a particualr face, plots them, plots the voronoi diagram
        fixes sink face and plots the other faces around it
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
        :param p_label_shift: how to shift the point labels if they exist
        :param point_names: names of the points, list or None

        :return: (all transitions shown (none if no points),
            whether we are done plotting (i.e. i_to_display is None or larger than the number of paths))
        """

        voronoi_diagram = self.get_voronoi_diagram(p=p,
                                                   source_fn=source_fn,
                                                   sink_fn=sink_fn,
                                                   diameter=diameter,
                                                   do_filter=do_filter,
                                                   intersect_with_face=False,
                                                   ignore_points_on_locus=ignore_points_on_locus,
                                                   )
        if voronoi_diagram is None:
            # cut locus does not exist on this face
            return None, True

        (point_pair_to_segment,
         (relevant_points, relevant_bound_paths, relevant_cells)
         ) = voronoi_diagram
        if point_names is None:
            point_names = []
        while len(point_names) < len(relevant_points):
            point_names.append(str(len(point_names)))
        if ax is None:
            ax = plt.gca()
        source: Face = self.faces[source_fn]

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
                    self._plot_label_face(ax=ax,
                                          face=face,
                                          name=name,
                                          center=center,
                                          rot_v=rot_v,
                                          color='blue',
                                          linewidth=1,
                                          orient_string=orient_string,
                                          zorder=7,
                                          label_zorder=11,
                                          )
                label = None
                if not labeled:
                    label = '$p$ copies'
                    labeled = True
                pt_color = 'purple'
                ax.scatter(pt[0], pt[1],
                           color=pt_color,
                           label=label,
                           alpha=1,
                           s=40,
                           zorder=11,
                           )
                if label_diagram:
                    label_pt = pt + [[.1], [.1]]
                    pname = point_names[pt_idx]
                    ax.annotate('$p^{(' + str(pname) + ')}$',
                                (label_pt[0] + p_label_shift[0], label_pt[1] + p_label_shift[1]),
                                rotation=0,
                                color=pt_color,
                                zorder=11,
                                )
        (face, name, center, rot_v) = special_face
        self._plot_label_face(ax=ax,
                              face=face,
                              name=name,
                              center=center,
                              rot_v=rot_v,
                              color='red',
                              linewidth=2,
                              orient_string=orient_string,
                              zorder=8,
                              label_zorder=11,
                              )
        # ax.scatter([0], [0], label='center', alpha=.5, s=80)
        relevant_points = np.concatenate(relevant_points, axis=1)

        xlim, ylim = ax.get_xlim(), ax.get_ylim()

        line_alpha = (.69 if label_diagram else 1)
        existing_labels = []
        for point_pair in point_pair_to_segment:
            segtype, (a, b) = point_pair_to_segment[point_pair]
            if segtype == 'segment':
                ax.plot((a[0], b[0]), (a[1], b[1]),
                        color='black',
                        lw=2,
                        alpha=line_alpha,
                        zorder=10,
                        )
                names = [point_names[pt_idx] for pt_idx in point_pair]
                names.sort()
                a, b = get_correct_end_points(a, b, xlim, ylim)
                if label_diagram and a is not None:
                    midpoint = (a + b)/2
                    tangent = np.array([(b - a)[1], -(b - a)[0]])
                    tangent = tangent/np.linalg.norm(tangent)*line_label_dist
                    if np.dot(tangent, np.array((3, 1))) < 0:  # favor the up right direction (strong right)
                        tangent = -tangent
                    if existing_labels:
                        temp_tan = tangent.copy()
                        for mult in [1] + sum([[i, -i] for i in range(3, 5)], []):
                            tangent = mult*temp_tan
                            dist = min([np.linalg.norm(lb - (midpoint + tangent)) for lb in existing_labels])
                            if dist >= 2*line_label_dist:
                                break

                    text_pt = midpoint + tangent
                    ax.annotate('$\\mathbf{\\ell}^{' + '\\{' + names[0] + ',' + names[1] + '\\}' + '}$',
                                text_pt, rotation=0, color='black',
                                zorder=11,
                                )
                    ax.arrow(
                        text_pt[0],  # x
                        text_pt[1],  # y
                        -tangent[0],  # dx
                        -tangent[1],  # dy
                        width=.03,
                        color='black',
                        alpha=line_alpha,
                        length_includes_head=True,
                        zorder=11,
                    )
                    existing_labels.append(text_pt)
                    # TODO: label lines here
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.legend()
        return all_trans_shown, i_to_display == None

    def plot_voronoi(self,
                     p,
                     source_fn,
                     sink_fn,
                     diameter,
                     ax,
                     do_filter=True,
                     plot_endpoints=False,
                     zorder=None,
                     ignore_points_on_locus=False,
                     ):
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
        voronoi_diagram = self.get_voronoi_diagram(p=p,
                                                   source_fn=source_fn,
                                                   sink_fn=sink_fn,
                                                   diameter=diameter,
                                                   do_filter=do_filter,
                                                   intersect_with_face=True,
                                                   ignore_points_on_locus=ignore_points_on_locus,
                                                   )
        if voronoi_diagram is not None:
            point_pair_to_seg, _ = voronoi_diagram
            for point_pair in point_pair_to_seg:
                _, (p, q) = point_pair_to_seg[point_pair]
                ax.plot((p[0], q[0]), (p[1], q[1]), color='black', lw=2, alpha=1, zorder=zorder)
                if plot_endpoints:
                    ax.scatter((p[0], q[0]), (p[1], q[1]), color='black', alpha=1, s=6.9, zorder=zorder)
            return True
        return False

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

    def plot_face_boundaries(self, axs, legend, legend_kwargs=None):
        """
        plots faces and points/arcs on faces

        :param axs: axis to plot on
        :param legend: (i,j)-> whether to put a legend on plot (i,j)
        """
        face_map, n, m = self.faces_to_plot_n_m()
        if legend_kwargs is None:
            legend_kwargs = dict()
        if "face_name_to_label" in legend_kwargs:
            face_name_to_label = legend_kwargs.pop("face_name_to_label")
        else:
            face_name_to_label = lambda x: x

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
                        label = face_name_to_label(str(f.name))
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

                    if legend(i, j):
                        ploot(i, j).legend(**legend_kwargs)

    def interactive_vornoi_plot(self,
                                figsize=None,
                                legend=lambda i, j: False,
                                diameter=None,
                                event_key='button_press_event',
                                source_fn_p=None,
                                show=True,
                                save=None,
                                save_kwargs=None,
                                do_filter=True,
                                font_size=None,
                                ignore_points_on_locus=False,
                                mark_points=(),
                                legend_kwargs=None,
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
        :param mark_points: points to always mark, list of (face id, x, y, color)
        """
        if save_kwargs is None:
            save_kwargs = dict()
        mark_dict = {}
        for fid, mp, c in mark_points:
            if fid not in mark_dict:
                mark_dict[fid] = []
            mark_dict[fid].append((mp, c))
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
            self.plot_face_boundaries(axs, legend=legend, legend_kwargs=legend_kwargs)
            ax.scatter(p[0, 0], p[1, 0], color='purple')

            source_fn = fc.name

            for i in range(n):
                for j in range(m):
                    face = face_map(i, j)
                    if face is not None:
                        xlim, ylim = ploot(i, j).get_xlim(), ploot(i, j).get_ylim()
                        self.plot_voronoi(p,
                                          source_fn,
                                          face.name,
                                          diameter=diameter,
                                          ax=ploot(i, j),
                                          do_filter=do_filter,
                                          plot_endpoints=False,
                                          zorder=10,
                                          ignore_points_on_locus=ignore_points_on_locus,
                                          )
                        for (mpx, mpy), c in mark_dict.get(str(face.name), []):
                            if c is not None:
                                c = c.replace('\\', '')
                            ploot(i, j).scatter(mpx, mpy, color=c)
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
            self.plot_face_boundaries(axs, legend=legend, legend_kwargs=legend_kwargs)
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
            plt.savefig(save, **save_kwargs)
            print('saving to', save)
        if show:
            plt.show()

    def interactive_unfold(self,
                           figsize=None,
                           legend=lambda i, j: False,
                           diameter=None,
                           track=True,
                           single_display=True,
                           source_fn_p=None,
                           sink_fn=None,
                           show=True,
                           save=None,
                           save_kwargs=None,
                           orient_string='',
                           do_filter=True,
                           font_size=None,
                           label_diagram=False,
                           p_label_shift=(0., 0.),
                           line_label_dist=.3,
                           point_names=None,
                           voronoi_star=False,
                           ignore_points_on_locus=False,
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
        :param p_label_shift: how to shift the point labels if they exist
        :param point_names: names of the points, list or None
        """
        if save_kwargs is None:
            save_kwargs = dict()
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
            """
            :return: whether or not done plotting
            """
            if ((self.extra_data['unwrap_source_fn'] is not None and self.extra_data['unwrap_sink_fn'] is not None) or
                    (voronoi_star and self.extra_data['unwrap_source_fn'] is not None)):
                # if we have finished both
                plt.clf()
                i_to_display = None
                if single_display:
                    i_to_display = self.extra_data['unwrap_counter']
                if voronoi_star:
                    self.plot_voronoi_star_unfolding(p=self.extra_data['p'],
                                                     source_fn=self.extra_data['unwrap_source_fn'],
                                                     ax=plt.gca(),
                                                     do_filter=do_filter,
                                                     diameter=diameter,
                                                     ignore_points_on_locus=ignore_points_on_locus,
                                                     )

                    plt.xticks([])
                    plt.yticks([])
                    if save is not None:
                        plt.savefig(save, **save_kwargs)
                else:
                    all_trans_shown, done_plotting = self.plot_unfolding_from_sink(
                        p=self.extra_data['p'],
                        source_fn=self.extra_data['unwrap_source_fn'],
                        sink_fn=self.extra_data['unwrap_sink_fn'],
                        diameter=diameter,
                        ax=plt.gca(),
                        i_to_display=i_to_display,
                        orient_string=orient_string,
                        do_filter=do_filter,
                        label_diagram=label_diagram,
                        p_label_shift=p_label_shift,
                        line_label_dist=line_label_dist,
                        point_names=point_names,
                        ignore_points_on_locus=ignore_points_on_locus,
                    )
                    print('point locations:')
                    for i, (zero, xvec, yvec, p) in enumerate(all_trans_shown):
                        if point_names is not None and i < len(point_names):
                            pname = point_names[i]
                        else:
                            pname = i
                        print('p copy ' + str(pname) + ':', p.flatten())
                        # since mostly this is of the form sqrt(int), try making this nice
                        shiftt = np.array([int(np.sign(val))*sym.sqrt(round(val**2))
                                           if abs(round(val**2) - val**2) <= 1e9
                                           else val
                                           for val in zero.flatten()])
                        print('\tshift:', shiftt)
                        rot_frac = fractions.Fraction(
                            np.arctan2((xvec - zero)[1], (xvec - zero)[0])[0]/np.pi).limit_denominator(1000)
                        print("\trotation:", rot_frac, 'pi')

                        print("\tx vec:", (xvec - zero).flatten())
                        print("\ty vec:", (yvec - zero).flatten())

                    plt.xticks([])
                    plt.yticks([])
                    if single_display and (all_trans_shown is not None):
                        if save is not None:
                            if not done_plotting:
                                plt.savefig(os.path.join(os.path.dirname(save),
                                                         'point_' + str(i_to_display) + '_' + os.path.basename(save)),
                                            **save_kwargs)
                            if done_plotting:
                                plt.savefig(save, **save_kwargs)
                        plt.title("click to advance")
                    self.extra_data['unwrap_counter'] += 1
                    return done_plotting
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
                                self.plot_voronoi(self.extra_data['p'],
                                                  self.extra_data['unwrap_source_fn'], face.name,
                                                  diameter=diameter,
                                                  ax=ploot(i, j),
                                                  do_filter=do_filter,
                                                  ignore_points_on_locus=ignore_points_on_locus,
                                                  )
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
                return False

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
                        self.plot_voronoi(p,
                                          source_fn,
                                          face.name,
                                          diameter=diameter,
                                          ax=ploot(i, j),
                                          do_filter=do_filter,
                                          ignore_points_on_locus=ignore_points_on_locus,
                                          )
                        ploot(i, j).set_xlim(xlim)
                        ploot(i, j).set_ylim(ylim)
            plt.show()

        cid = fig.canvas.mpl_connect('button_press_event', clicked_mouse)
        if track:
            cid2 = fig.canvas.mpl_connect('motion_notify_event', moved_mouse)
        self.plot_face_boundaries(axs, legend=legend)
        plt.suptitle("Click $p$")
        done_plottin = spin()
        if save is not None:
            plt.savefig(save, **save_kwargs)
            if single_display:
                while not done_plottin:
                    done_plottin = spin()
        if show:
            plt.show()

    def plot_faces(self,
                   save_image=None,
                   show=False,
                   figsize=None,
                   legend=lambda i, j: True,
                   voronoi=None,
                   do_filter=True,
                   ignore_points_on_locus=False,
                   ):
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
                        ignore_locus_points = self.plot_voronoi(p,
                                                                source_fn,
                                                                face.name,
                                                                diameter=diameter,
                                                                ax=ploot(i, j),
                                                                do_filter=do_filter,
                                                                ignore_points_on_locus=ignore_points_on_locus,
                                                                )
                        ploot(i, j).set_xlim(xlim)
                        ploot(i, j).set_ylim(ylim)

        if save_image is not None:
            plt.savefig(save_image)
        if show:
            plt.show()
        plt.close()
