from utils.shape_argparser import *

PARSER.add_argument("--sink-face-name", action='store', required=False, default=None,
                    help="Specify sink face name if inputting a specific sink face")

PARSER.add_argument("--single-display", action='store_true', required=False,
                    help="display only one path at a time")

PARSER.add_argument("--no-tracking", action='store_true', required=False,
                    help="stop tracking the cut locus as mouse moves")

PARSER.add_argument("--orient", action='store', required=False, default='',
                    help="Specify a string to append onto face name to help with orientation ('_' usually works well)")

args = PARSER.parse_args()
shape = shape_from_args(args)

source_fn_p = get_source_fn_p_from_args(args, shape)
sink_face_name = args.sink_face_name

if sink_face_name is not None:
    temp = sink_face_name
    sink_face_name = check_face_name(sink_face_name, shape)
    if sink_face_name is None:
        raise Exception("invalid face name specified: " + str(temp))

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
                         do_filter=shape.is_polyhedra(),
                         )
