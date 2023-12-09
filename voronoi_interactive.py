from src.shapes import Shape
from src.shape_creation import (Tetrahedron,
                                Cube,
                                Octahedron,
                                Dodecahedron,
                                Icosahedron,
                                Prism,
                                Antiprism,
                                Pyramid,
                                ElongatedPyramid,
                                Bipyramid,
                                ElongatedBipyramid,
                                Large2Torus
                                )
import argparse
import numpy as np

mapping = {'tetrahedron': Tetrahedron,
           'cube': Cube,
           'octahedron': Octahedron,
           'dodecahedron': Dodecahedron,
           'icosahedron': Icosahedron,
           'pyramid': Pyramid,
           'bipyramid': Bipyramid,
           'longpyramid': ElongatedPyramid,
           'longbipyramid': ElongatedBipyramid,
           'prism': Prism,
           'antiprism': Antiprism,
           'torus': Large2Torus}

arg_n = ('prism', 'antiprism', 'pyramid', 'longpyramid', 'bipyramid', 'longbipyramid')
# shapes that take n as an arg

PARSER = argparse.ArgumentParser()
PARSER.add_argument("-s", "--shape", action='store', required=True,
                    help="Specify which face to display, options are " + str(tuple(s for s in mapping)))
PARSER.add_argument("--n", type=int, required=False, default=3,
                    help="additional argument to specify n, used for " + str(arg_n))

PARSER.add_argument("--click", action='store_true', required=False,
                    help="toggle whether to click to update point or have point update with mouse")
PARSER.add_argument("--legend", action='store_true', required=False,
                    help="toggle whether to put legend on plot")
PARSER.add_argument("--diameter", type=int, required=False, default=-1,
                    help="Specify diameter of search graph (longest possible sequence of faces on a geodesic)")

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
