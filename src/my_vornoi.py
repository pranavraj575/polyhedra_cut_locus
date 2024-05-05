import numpy as np

__all__ = ['voronoi_plot_2d']


def _adjust_bounds(ax, points):
    margin = 0.1*points.ptp(axis=0)
    xy_min = points.min(axis=0) - margin
    xy_max = points.max(axis=0) + margin
    ax.set_xlim(xy_min[0], xy_max[0])
    ax.set_ylim(xy_min[1], xy_max[1])


def voronoi_plot_2d(vor, ax=None, **kw):
    """
    Plot the given Voronoi diagram in 2-D

    Parameters
    ----------
    vor : scipy.spatial.Voronoi instance
        Diagram to plot
    ax : matplotlib.axes.Axes instance, optional
        Axes to plot on
    show_points : bool, optional
        Add the Voronoi points to the plot.
    show_vertices : bool, optional
        Add the Voronoi vertices to the plot.
    label_lines : (point, point)->str, optional
        label lines in plot, function is point pair to the name
    line_colors : string, optional
        Specifies the line color for polygon boundaries
    line_width : float, optional
        Specifies the line width for polygon boundaries
    line_alpha : float, optional
        Specifies the line alpha for polygon boundaries
    point_size : float, optional
        Specifies the size of points

    Returns
    -------
    fig : matplotlib.figure.Figure instance
        Figure for the plot
    points_to_segments: dictionary of point index to line segments the point creates
        segments are represetned as
    """

    def within_lim(value, lim):
        return (value <= lim[1] and
                value >= lim[0])

    def within_bounds(point, xlim=None, ylim=None):
        if xlim is None:
            xlim = ax.get_xlim()
        if ylim is None:
            ylim = ax.get_ylim()
        return (within_lim(point[0], xlim) and
                within_lim(point[1], ylim))

    def get_correct_end_points(start, end, xlim, ylim):
        end = get_correct_exit_point(start, end, xlim, ylim)
        if end is None:
            return None, None
        start = get_correct_exit_point(end, start, xlim, ylim)
        if start is None:
            return None, None
        return start, end

    def get_correct_exit_point(start, end, xlim, ylim):
        vec = end - start
        if vec[0] > 0:
            if (start + vec)[0] > xlim[1]:
                if start[0] > xlim[1]:
                    return None
                vec = vec*((xlim[1] - start[0])/vec[0])
        elif vec[0] < 0:
            if (start + vec)[0] < xlim[0]:
                if start[0] < xlim[0]:
                    return None
                vec = vec*((start[0] - xlim[0])/(-vec[0]))
        elif vec[0] == 0:
            if not within_lim(start[0], xlim):
                return None

        if vec[1] > 0:
            if (start + vec)[1] > ylim[1]:
                if start[1] > ylim[1]:
                    return None
                vec = vec*((ylim[1] - start[1])/vec[1])
        elif vec[1] < 0:
            if (start + vec)[1] < ylim[0]:
                if start[1] < ylim[0]:
                    return None
                vec = vec*((start[1] - ylim[0])/(-vec[1]))
        elif vec[1] == 0:
            if not within_lim(start[1], ylim):
                return None
        return vec + start

    from matplotlib.collections import LineCollection

    if vor.points.shape[1] != 2:
        raise ValueError("Voronoi diagram is not 2-D")

    if ax is not None:
        if kw.get('show_points', True):
            point_size = kw.get('point_size', None)
            ax.plot(vor.points[:, 0], vor.points[:, 1], '.', markersize=point_size)
        if kw.get('show_vertices', True):
            ax.plot(vor.vertices[:, 0], vor.vertices[:, 1], 'o')

        line_colors = kw.get('line_colors', 'k')
        line_width = kw.get('line_width', 1.0)
        line_alpha = kw.get('line_alpha', 1.0)
        line_label_dist = kw.get('line_label_dist', .3)

    center = vor.points.mean(axis=0)
    ptp_bound = vor.points.ptp(axis=0)

    finite_segments = []
    infinite_segments = []
    points_to_segments = {i: list() for i in
                          range(vor.npoints)}  # dictionary of point indices to the segments they create
    text_positions = []  # list of label positions so that the next one is far from the previous
    for pointidx, simplex in zip(vor.ridge_points, vor.ridge_vertices):  # iterates through all lines
        # pointidx: two indices of points that create this line
        # simplex: two indices of vertices of the vornoi diagram that create this line
        #       if there is an infinite vertex, there is a -1 here
        simplex = np.asarray(simplex)
        avg_point = (vor.points[pointidx[1]] + vor.points[pointidx[0]])/2
        # midpoint the two points whose bisection forms the line
        # for labeling purposes
        tangeant = vor.points[pointidx[1]] - vor.points[pointidx[0]]
        tangeant = tangeant/np.linalg.norm(tangeant)
        if np.dot((1, 1), tangeant) < 0:
            tangeant = -tangeant

        if np.all(simplex >= 0):
            finite_segments.append(vor.vertices[simplex])
            for idx in pointidx:
                points_to_segments[idx].append((vor.vertices[simplex[0]], vor.vertices[simplex[1]]))
            if ax is not None:
                a, b = get_correct_end_points(vor.vertices[simplex[0]], vor.vertices[simplex[1]], ax.get_xlim(),
                                              ax.get_ylim())
                if a is None:
                    label_point = None
                else:
                    label_point = (a + b)/2

        else:
            i = simplex[simplex >= 0][0]  # finite end Voronoi vertex

            t = vor.points[pointidx[1]] - vor.points[pointidx[0]]  # tangent
            t /= np.linalg.norm(t)
            n = np.array([-t[1], t[0]])  # normal

            midpoint = vor.points[pointidx].mean(axis=0)
            direction = np.sign(np.dot(midpoint - center, n))*n
            if (vor.furthest_site):
                direction = -direction
            far_point = vor.vertices[i] + direction*ptp_bound.max()

            infinite_segments.append([vor.vertices[i], far_point])
            for idx in pointidx:
                points_to_segments[idx].append((vor.vertices[i], far_point))
            if ax is not None:
                a, b = get_correct_end_points(vor.vertices[i], far_point, ax.get_xlim(),
                                              ax.get_ylim())
                if a is None:
                    label_point = None
                else:
                    label_point = (a + b)/2

        if ax is not None and kw.get('label_lines', False) and label_point is not None:
            # push label point off line
            potential_text_points = [
                                        label_point + tangeant*(1 + dx/10)*line_label_dist for dx in range(10)
                                    ] + [
                                        label_point - tangeant*(3 + dx/10)*line_label_dist for dx in range(10)
                                    ]
            if len(text_positions) > 0:
                satisfiable = False
                for i in range(len(potential_text_points)):
                    text_point = potential_text_points[i]
                    if min(np.linalg.norm(np.array(text_positions) - text_point, axis=1)) > 1.5:
                        satisfiable = True
                        break
                if not satisfiable:
                    text_point = max(potential_text_points,
                                     key=lambda pt: min(np.linalg.norm(np.array(text_positions) - pt, axis=1)))
            else:
                text_point = potential_text_points[0]

            text_positions.append(text_point)
            label_names = sorted([str(pointidx[0]), str(pointidx[1])])
            if within_bounds(text_point):
                ax.annotate('$\\mathbf{\\ell}^{' + '\{' + label_names[0] + ',' + label_names[1] + '\}' + '}$',
                            (text_point[0], text_point[1]), rotation=0, color=line_colors)
                diff = (label_point - text_point)
                ax.arrow(
                    text_point[0],  # x
                    text_point[1],  # y
                    diff[0],  # dx
                    diff[1],  # dy
                    width=.03,
                    color=line_colors,
                    alpha=line_alpha,
                    length_includes_head=True,
                )

    fig = None
    if ax is not None:
        ax.add_collection(LineCollection(finite_segments,
                                         colors=line_colors,
                                         lw=line_width,
                                         alpha=line_alpha,
                                         linestyle='solid'))
        ax.add_collection(LineCollection(infinite_segments,
                                         colors=line_colors,
                                         lw=line_width,
                                         alpha=line_alpha,
                                         linestyle='solid'))

        _adjust_bounds(ax, vor.points)
        fig = ax.figure
    return fig, points_to_segments
