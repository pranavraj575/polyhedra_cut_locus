from src.shape_argparser import *
import numpy as np

PARSER.add_argument("--sink-face-name", action='store', required=False, default=None,
                    help="Specify sink face name if inputting a specific sink face")

PARSER.add_argument("--single-display", action='store_true', required=False,
                    help="display only one path at a time")

PARSER.add_argument("--no-tracking", action='store_true', required=False,
                    help="stop tracking the cut locus as mouse moves")

args = PARSER.parse_args()
if args.shape not in mapping:
    raise Exception("ERROR: haven't made shape '" + args.shape + "' yet, valid arguments are " + str(tuple(s for s in mapping)))
SHAPE = mapping[args.shape]
if args.shape in arg_n:
    shape = SHAPE(args.n)
else:
    shape = SHAPE()

source_face_name = args.face_name
p = np.array([[args.point_x], [args.point_y]])
source_fn_p = None

if source_face_name is not None:
    temp = source_face_name
    source_face_name = check_face_name(source_face_name, shape)
    if source_face_name is None:
        raise Exception("invalid file name specified: " + str(temp))
    source_fn_p = (source_face_name, p)
sink_face_name = args.sink_face_name

if sink_face_name is not None:
    temp = sink_face_name
    sink_face_name = check_face_name(sink_face_name, shape)
    if sink_face_name is None:
        raise Exception("invalid file name specified: " + str(temp))

shape.interactive_unwrap(track=not args.no_tracking,
                         single_display=args.single_display,
                         diameter=args.diameter if args.diameter > 0 else None,
                         legend=lambda i, j: args.legend,
                         source_fn_p=source_fn_p,
                         sink_fn=sink_face_name,
                         )
# ax= plt.gca()
# shape=Octahedron()
# shape.plot_unwrapping(np.zeros((2,1))+.1,1,7,None,ax)
# shape.plot_unwrapping(np.array([[-.4], [-.6]]),3,0,None,ax)
# plt.show()
