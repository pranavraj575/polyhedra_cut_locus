from src.shapes import *


class ConvexPolyhderon(Shape):
    """
    Geodesics on the surface of convex polyhedra
    """

    def __init__(self, tolerance):
        super().__init__(tolerance)

    def is_polyhedra(self):
        return True

    def plot_unwrapping(self, p, source_fn, sink_fn, diameter, ax,
                        i_to_display=None,
                        orient_string='',
                        do_filter=True,
                        label_diagram=False,
                        p_label_shift=(0., 0.),
                        line_label_dist=.3,
                        point_names=None,
                        ):
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
        :param p_label_shift: how to shift the point labels if they exist
        :param point_names: names of the points, list or None

        :return: (all transitions shown (none if no points),
            whether we are done plotting (i.e. i_to_display is None or larger than the number of paths))
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
            ax.annotate(str(name) + orient_string, (center[0], center[1]), rotation=np.degrees(theta), color=color)

        if len(vp) >= 2:  # if there is only one, the cut locus does not exist on this face
            relevant_points, relevant_bound_paths, relevant_cells = self.filter_out_points(vp, bound_paths, source,
                                                                                           sink, do_filter=do_filter)

            if relevant_points is None:
                return None, True
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
                        plot_label_face(ax=ax,
                                        face=face,
                                        name=name,
                                        center=center,
                                        rot_v=rot_v,
                                        color='blue',
                                        linewidth=1)
                    label = None
                    if not labeled:
                        label = '$p$ copies'
                        labeled = True
                    pt_color = 'purple'
                    ax.scatter(pt[0], pt[1],
                               color=pt_color,
                               label=label,
                               alpha=1,
                               s=40)
                    if label_diagram:
                        label_pt = pt + [[.1], [.1]]
                        if point_names is not None and pt_idx < len(point_names):
                            pname = point_names[pt_idx]
                        else:
                            pname = pt_idx
                        ax.annotate('$p^{(' + str(pname) + ')}$',
                                    (label_pt[0] + p_label_shift[0], label_pt[1] + p_label_shift[1]),
                                    rotation=0,
                                    color=pt_color)
            (face, name, center, rot_v) = special_face
            plot_label_face(ax=ax,
                            face=face,
                            name=name,
                            center=center,
                            rot_v=rot_v,
                            color='red',
                            linewidth=2)
            # ax.scatter([0], [0], label='center', alpha=.5, s=80)
            relevant_points = np.concatenate(relevant_points, axis=1)

            vor = Voronoi(relevant_points.T)
            xlim, ylim = ax.get_xlim(), ax.get_ylim()
            voronoi_plot_2d(vor,
                            ax=ax,
                            show_points=False,
                            show_vertices=False,
                            line_colors='black',
                            line_width=2,
                            line_alpha=(.69 if label_diagram else 1),
                            label_lines=label_diagram,
                            line_label_dist=line_label_dist,
                            point_names=point_names,
                            )
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            ax.legend()
            return all_trans_shown, i_to_display == None
        return None, True

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
                           p_label_shift=(0., 0.),
                           line_label_dist=.3,
                           point_names=None,
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
            if self.extra_data['unwrap_source_fn'] is not None and self.extra_data['unwrap_sink_fn'] is not None:
                # if we have finished both
                plt.clf()
                i_to_display = None
                if single_display:
                    i_to_display = self.extra_data['unwrap_counter']
                all_trans_shown, done_plotting = self.plot_unwrapping(self.extra_data['p'],
                                                                      self.extra_data['unwrap_source_fn'],
                                                                      self.extra_data['unwrap_sink_fn'],
                                                                      diameter=diameter, ax=plt.gca(),
                                                                      i_to_display=i_to_display,
                                                                      orient_string=orient_string, do_filter=do_filter,
                                                                      label_diagram=label_diagram,
                                                                      p_label_shift=p_label_shift,
                                                                      line_label_dist=line_label_dist,
                                                                      point_names=point_names,
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
                                                     'point_' + str(i_to_display) + '_' + os.path.basename(save)))
                        if done_plotting:
                            plt.savefig(save)
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
        done_plottin = spin()
        if save is not None:
            plt.savefig(save)
            if single_display:
                while not done_plottin:
                    done_plottin = spin()
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
