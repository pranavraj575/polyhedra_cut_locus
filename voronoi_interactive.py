from utils.shape_argparser import *

PARSER.add_argument("--click", action='store_true', required=False,
                    help="click to update point instead of updating with movement")

PARSER.add_argument("--center-pt", action='store_true', required=False,
                    help="add center point to faces")

args = PARSER.parse_args()
shape = shape_from_args(args)

if args.center_pt:
    for fn in shape.faces:
        shape.add_point_to_face(np.zeros((2, 1)), fn, {'color': 'black', 's': 1})

source_fn_p = get_source_fn_p_from_args(args, shape)
event_key = None

if source_fn_p is None:
    event_key = 'button_press_event' if args.click else 'motion_notify_event'

shape.interactive_vornoi_plot(diameter=args.diameter if args.diameter > 0 else None,
                              figsize=figsize_from_args(args),
                              event_key=event_key,
                              legend=lambda i, j: args.legend,
                              source_fn_p=source_fn_p,
                              show=not args.no_show,
                              save=args.save_file,
                              )
