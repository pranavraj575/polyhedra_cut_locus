
from src.shape_argparser import *

args = PARSER.parse_args()
if args.shape not in mapping:
    raise Exception("ERROR: haven't made shape '" + args.shape + "' yet, valid arguments are " + str(tuple(s for s in mapping)))
SHAPE = mapping[args.shape]
if args.shape in arg_n:
    shape = SHAPE(args.n)
else:
    shape = SHAPE()

shape.interactive_unwrap()
#ax= plt.gca()
#shape=Octahedron()
#shape.plot_unwrapping(np.zeros((2,1))+.1,1,7,None,ax)
#shape.plot_unwrapping(np.array([[-.4], [-.6]]),3,0,None,ax)
#plt.show()
