import numpy as np

from src.shape_argparser import *

PARSER.add_argument("--center_pt", action='store_true', required=False,
                    help="toggle whether to add center point to faces")

args = PARSER.parse_args()
if args.shape not in mapping:
    raise Exception("ERROR: haven't made shape '" + args.shape + "' yet, valid arguments are " + str(tuple(s for s in mapping)))
SHAPE = mapping[args.shape]
if args.shape in arg_n:
    shape = SHAPE(args.n)
else:
    shape = SHAPE()
if args.center_pt:
    for fn in shape.faces:
        shape.add_point_to_face(np.zeros((2, 1)), fn, {'color': 'black', 's': 1})
shape.interactive_vornoi_plot(diameter=args.diameter if args.diameter > 0 else None,
                              event_key='button_press_event' if args.click else 'motion_notify_event',
                              legend=lambda i, j: args.legend)
