from utils.shape_argparser import *

PARSER.add_argument("--sink-face-name", action='store', required=False, default=None,
                    help="Specify sink face name if inputting a specific sink face")

display_group.add_argument("--single-display", action='store_true', required=False,
                           help="display only one path at a time")

display_group.add_argument("--orient", action='store', required=False, default='',
                           help="Specify a string to append onto face name to help with orientation ('_' usually works well)")

display_group.add_argument("--label-unwrapping", action='store_true', required=False,
                           help="whether to label the points and lines on the unwrapping diagram")

display_group.add_argument("--shift-x-p-label", type=float, required=False, default=0.,
                           help="shift each p label if necessary")
display_group.add_argument("--shift-y-p-label", type=float, required=False, default=0.,
                           help="shift each p label if necessary")
display_group.add_argument("--label-dist-line", type=float, required=False, default=.3,
                           help="distance to label each line from")
display_group.add_argument("--point-names", action='store', required=False, default=None,
                           help="comma separated names of each point")

args = parse_args(PARSER)
shape = shape_from_args(args)
point_names = None
if args.point_names is not None:
    point_names = args.point_names.split(',')

source_fn_p = get_source_fn_p_from_args(args, shape)
sink_face_name = args.sink_face_name

if sink_face_name is not None:
    temp = sink_face_name
    sink_face_name = check_face_name(sink_face_name, shape)
    if sink_face_name is None:
        raise Exception("invalid face name specified: " + str(temp))

do_filter = shape.is_polyhedra() and not args.no_filter

shape.interactive_unwrap(track=not args.no_tracking,
                         figsize=figsize_from_args(args),
                         single_display=args.single_display,
                         diameter=args.diameter if args.diameter > 0 else None,
                         legend=lambda i, j: args.legend,
                         source_fn_p=source_fn_p,
                         sink_fn=sink_face_name,
                         show=not args.no_show,
                         save=args.save_file,
                         orient_string=args.orient,
                         do_filter=do_filter,
                         font_size=args.font_size,
                         label_diagram=args.label_unwrapping,
                         p_label_shift=(args.shift_x_p_label, args.shift_y_p_label),
                         line_label_dist=args.label_dist_line,
                         point_names=point_names
                         )
