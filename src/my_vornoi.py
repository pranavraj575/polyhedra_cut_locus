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

    center = vor.points.mean(axis=0)
    ptp_bound = vor.points.ptp(axis=0)

    finite_segments = []
    infinite_segments = []
    points_to_segments = {i: list() for i in
                          range(vor.npoints)}  # dictionary of point indices to the segments they create
    text_positions = []  # list of label positions so that the next one is far from the previous
    arrows = []
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
            potential = (vor.vertices[simplex[0]] + vor.vertices[simplex[1]])/2
            if ax is not None and (potential[0] <= ax.get_xlim()[1] and
                                   potential[0] >= ax.get_xlim()[0] and
                                   potential[1] <= ax.get_ylim()[1] and
                                   potential[1] >= ax.get_ylim()[0]):
                label_point = potential
            else:
                label_point = avg_point
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

            potential = (vor.vertices[simplex[0]] + vor.vertices[simplex[1]])/2
            if ax is not None and (potential[0] <= ax.get_xlim()[1] and
                                   potential[0] >= ax.get_xlim()[0] and
                                   potential[1] <= ax.get_ylim()[1] and
                                   potential[1] >= ax.get_ylim()[0]):
                label_point = potential
            else:
                label_point = avg_point

        if ax is not None and kw.get('label_lines', False):
            # push label point off line
            potential_text_points = [label_point + tangeant*.3,
                                     label_point + tangeant*.5,
                                     label_point - tangeant*.5,  # - np.array([.5, 0.])
                                     ]
            if len(text_positions) > 0:
                text_point = max(potential_text_points, key=lambda pt:
                min(np.linalg.norm(np.array(text_positions) - pt, axis=1)))

                dist = min(np.linalg.norm(np.array(text_positions) - text_point, axis=1))
                if dist > 1:
                    text_point = potential_text_points[0]
            else:
                text_point = potential_text_points[0]

            text_positions.append(text_point)
            label_names = sorted([str(pointidx[0]), str(pointidx[1])])
            ax.annotate('$\\mathbf{\\ell}^{' + '\{' + label_names[0] + ',' + label_names[1] + '\}' + '}$',
                        (text_point[0], text_point[1]), rotation=0, color=line_colors)
            if (text_point[0] <= ax.get_xlim()[1] and
                    text_point[0] >= ax.get_xlim()[0] and
                    text_point[1] <= ax.get_ylim()[1] and
                    text_point[1] >= ax.get_ylim()[0]):
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
