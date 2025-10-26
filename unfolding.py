import matplotlib.pyplot as plt
from utils.shape_argparser import *

PARSER.description = 'interactively view the unfoldings that create the cut locus of a point on the surface of a polytope'

PARSER.add_argument("--sink-face-name", action='store', required=False, default=None,
                    help="Specify sink face name if inputting a specific sink face")

PARSER.add_argument("--voronoi-star", action='store_true', required=False,
                    help="unfold from source face to create a voronoi star")
PARSER.add_argument("--ignore-points", action='store_true', required=False,
                    help="ignore single points on faces of cut loci, useful for fixing corner cases or repeat paths")

display_group.add_argument("--single-display", action='store_true', required=False,
                           help="display only one path at a time")

display_group.add_argument("--orient", action='store', required=False, default='',
                           help="Specify a string to append onto face name to help with orientation ('_' usually works well)")

display_group.add_argument("--label-unwrapping", action='store_true', required=False,
                           help="whether to label the points and lines on the unwrapping diagram")

display_group.add_argument("--shift-p-label", type=float, required=False, default=(0., 0.), nargs=2,
                           help="shift each p label if necessary", metavar=('x', 'y'))
display_group.add_argument("--label-dist-line", type=float, required=False, default=.3,
                           help="distance to label each line from")
display_group.add_argument("--point-names", action='store', nargs='*', required=False, default=None,
                           help="names of each point", metavar='p0 p1')
display_group.add_argument('--mark', action='append', nargs='*', required=False,
                           help='mark point on cut locus faces', metavar='X Y <COLOR>, <LABEL>')


def get_marks(args):
    marks = []
    if args.mark is None:
        return marks
    for mp in args.mark:
        if len(mp) < 2:
            raise Exception('--mark requires at least 2 args (--mark X Y), invalid:', mp)
        if len(mp) > 6:
            raise Exception('--mark takes at most 6 args (--mark X Y <COLOR> <LABEL> <SHIFT X> <SHIFT Y>), invalid:', mp)
        try:
            float(mp[0])
            float(mp[1])
        except:
            raise Exception('--mark usage is (--mark X Y <COLOR> <LABEL> <SHIFT X> <SHIFT Y>), invalid:', mp)

        mark = {
            'pt': (float(mp[0]), float(mp[1]))
        }
        if len(mp) >= 3:
            mark['color'] = (None if mp[2] == "None" else mp[2])

        if len(mp) >= 4:
            mark['label'] = (None if mp[3] == "None" else mp[3])
        if len(mp) == 6:
            mark['shift'] = (float(mp[4]), float(mp[5]))
        marks.append(mark)
    return marks


args = parse_args(PARSER)
shape = shape_from_args(args)
point_names = args.point_names
marks = get_marks(args)

source_fn_p = get_source_fn_p_from_args(args, shape)
sink_face_name = args.sink_face_name

if sink_face_name is not None:
    temp = sink_face_name
    sink_face_name = check_face_name(sink_face_name, shape)
    if sink_face_name is None:
        raise Exception("invalid face name specified: " + str(temp))

do_filter = shape.is_polyhedra() and not args.no_filter
# TODO: check

save_kwargs = dict()
if args.dpi is not None:
    save_kwargs['dpi'] = args.dpi
shape.interactive_unfold(track=not args.no_tracking,
                         figsize=figsize_from_args(args),
                         single_display=args.single_display,
                         diameter=args.diameter if args.diameter > 0 else None,
                         legend=lambda i, j: args.legend,
                         source_fn_p=source_fn_p,
                         sink_fn=sink_face_name,
                         show=not args.no_show,
                         save=args.save_file,
                         save_kwargs=save_kwargs,
                         orient_string=args.orient,
                         do_filter=do_filter,
                         font_size=args.font_size,
                         label_diagram=args.label_unwrapping,
                         p_label_shift=(args.shift_p_label[0], args.shift_p_label[1]),
                         line_label_dist=args.label_dist_line,
                         point_names=point_names,
                         voronoi_star=args.voronoi_star,
                         ignore_points_on_locus=args.ignore_points,
                         mark_points=marks,
                         )
