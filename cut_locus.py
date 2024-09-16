from utils.shape_argparser import *

display_group.add_argument("--center-pt", action='store_true', required=False,
                           help="add center point to faces")
PARSER.add_argument("--ignore-points", action='store_true', required=False,
                    help="ignore single points on faces of cut loci, useful for fixing corner cases or repeat paths")

args = parse_args(PARSER)
shape = shape_from_args(args)

if args.center_pt:
    for fn in shape.faces:
        shape.add_point_to_face(np.zeros((2, 1)), fn, {'color': 'black', 's': 1})

source_fn_p = get_source_fn_p_from_args(args, shape)
event_key = None

if source_fn_p is None:
    event_key = 'button_press_event' if args.no_tracking else 'motion_notify_event'

do_filter = shape.is_polyhedra() and not args.no_filter

shape.interactive_vornoi_plot(diameter=args.diameter if args.diameter > 0 else None,
                              figsize=figsize_from_args(args),
                              event_key=event_key,
                              legend=lambda i, j: args.legend,
                              source_fn_p=source_fn_p,
                              show=not args.no_show,
                              save=args.save_file,
                              do_filter=do_filter,
                              font_size=args.font_size,
                              ignore_points_on_locus=args.ignore_points,
                              )
