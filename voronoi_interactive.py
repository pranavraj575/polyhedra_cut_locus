import numpy as np

from src.shape_argparser import *

PARSER.add_argument("--click", action='store_true', required=False,
                    help="click to update point instead of updating with movement")

PARSER.add_argument("--center-pt", action='store_true', required=False,
                    help="add center point to faces")

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

face_name = args.face_name
p = np.array([[args.point_x], [args.point_y]])
source_fn_p = None
if face_name is not None:
    temp = face_name
    face_name = check_face_name(face_name, shape)
    if face_name is None:
        raise Exception("invalid file name specified: " + str(temp))
    source_fn_p = (face_name, p)

event_key = None
if source_fn_p is None:
    event_key = 'button_press_event' if args.click else 'motion_notify_event'
shape.interactive_vornoi_plot(diameter=args.diameter if args.diameter > 0 else None,
                              event_key=event_key,
                              legend=lambda i, j: args.legend,
                              source_fn_p=source_fn_p,
                              )
