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
                                Mirror,
                                Large2Torus
                                )
import sys, argparse
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
           'mirror': Mirror,
           'torus': Large2Torus}

arg_n = ('prism', 'antiprism', 'pyramid', 'longpyramid', 'bipyramid', 'longbipyramid', 'mirror')


# shapes that take n as an arg

def check_face_name(proposed_name, shape):
    """
    checks if proposed_name is actually a face name of shape
    converts everything to string for checking
        returns the face name
    :param proposed_name: (potentialy string) face name
    :param shape: Shape
    :return: correct face name
    """
    output = None
    for f in shape.faces:
        if str(proposed_name) == str(f):
            output = f
    return output


def figsize_from_args(args):
    """
    gets figsize from args
    :param args: args object
    :return: (width,height) or None
    """
    return args.display_dims


def get_source_fn_p_from_args(args, shape):
    """
    gets (source face name, point) from args
    :param args: args object
    :param shape: Shape
    :return: (source face name, column vector) or None
    """

    face_name = args.source_face
    p = np.array([[args.point[0]], [args.point[1]]])
    source_fn_p = None
    if face_name is not None:
        temp = face_name
        face_name = check_face_name(face_name, shape)
        if face_name is None:
            raise Exception("invalid face name specified: " + str(temp))
        source_fn_p = (face_name, p)
    return source_fn_p


def shape_from_args(args):
    """
    gets the shape object from args
    :param args: args object
    :return: Shape
    """
    possible = []
    if args.shape in mapping:
        possible = [args.shape]
    else:
        for name in mapping:
            if name.startswith(args.shape):
                possible.append(name)
    if len(possible) > 1:
        raise Exception("shape '" + args.shape + "' is ambiguous, could be " + str(tuple(possible)))
    if len(possible) == 0:
        raise Exception(
            "shape '" + args.shape + "' does not exist, valid arguments are " + str(tuple(s for s in mapping)))
    name = possible[0]
    SHAPE = mapping[name]
    if name in arg_n:
        n = args.n
        if n is None:
            raise Exception("shape '" + name + "' requires argument: [--n N]")

        if args.tolerance is None:
            return SHAPE(n)
        else:
            return SHAPE(n, tolerance=args.tolerance)
    else:
        if args.tolerance is None:
            return SHAPE()
        else:
            return SHAPE(tolerance=args.tolerance)


PARSER = argparse.ArgumentParser()
display_group = PARSER.add_argument_group('display', 'arguments to tweak display output')

PARSER.add_argument("-hv", action='store_true', required=False,
                    help="show help message with ALL the arguments and exit (display stuff included)")

PARSER.add_argument("-s", "--shape", action='store', required=True,
                    help="Specify which face to display, options are " + str(tuple(s for s in mapping)))
PARSER.add_argument("-n", "--n", type=int, required=False, default=None,
                    help="additional argument to specify n, used for " + str(arg_n))

PARSER.add_argument("--diameter", type=int, required=False, default=-1,
                    help="Specify diameter of search graph (longest possible sequence of faces on a geodesic)")
PARSER.add_argument("--no-filter", action='store_true', required=False,
                    help="Turn off filter on points of voronoi complex. " +
                         "This should fix tolerance errors, but may result in invalid points (might want to check with --single-display)")
PARSER.add_argument("--tolerance", type=float, required=False, default=None,
                    help="tolerance for things like intersection and containment, default differs for each shape")

PARSER.add_argument("--no-show", action='store_true', required=False,
                    help="don't show the plot")

PARSER.add_argument("--save-file", action='store', required=False, default=None,
                    help="save initial image as specified file")

display_group.add_argument("--legend", action='store_true', required=False,
                           help="put legend on plot")

display_group.add_argument("--display-dims", type=float, nargs=2, required=False, default=None,
                           help="dimenisions of display in inches", metavar=('WIDTH', 'HEIGHT'))
display_group.add_argument("--font-size", type=int, required=False, default=None,
                           help="font size for plotting")

PARSER.add_argument("--source-face", action='store', required=False, default=None,
                    help="Specify the face name if inputting a specific point")
PARSER.add_argument("--point", type=float, required=False, default=(0., 0.), nargs=2,
                    help="Specify point if inputting a specific point (defaults to (0,0))", metavar=('x', 'y'))

PARSER.add_argument("--no-tracking", action='store_true', required=False,
                    help="Click to interact instead of tracking the cut locus as mouse moves")


def parse_args(parser):
    if any(help_string in sys.argv for help_string in ['-h', '--help', "-hv"]):
        if not any(verboseness in sys.argv for verboseness in ["-hv", ]):
            # we have to remove the display group manually
            parser._action_groups.remove(display_group)
        # do this manually
        parser.print_help()
        quit()
    return parser.parse_args()
