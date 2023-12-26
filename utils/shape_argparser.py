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
    figsize = None
    if args.height_display is not None:
        if args.width_display is None:
            figsize = (args.height_display, args.height_display)
        else:
            figsize = (args.width_display, args.height_display)
    if args.width_display is not None:
        if args.height_display is None:
            figsize = (args.width_display, args.width_display)
        else:
            figsize = (args.width_display, args.height_display)
    return figsize


def get_source_fn_p_from_args(args, shape):
    """
    gets (source face name, point) from args
    :param args: args object
    :param shape: Shape
    :return: (source face name, column vector) or None
    """

    face_name = args.face_name
    p = np.array([[args.point_x], [args.point_y]])
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
        raise Exception("shape '" + args.shape + "' does not exist, valid arguments are " + str(tuple(s for s in mapping)))
    name = possible[0]
    SHAPE = mapping[name]
    if name in arg_n:
        n = args.n
        if n is None:
            raise Exception("shape '" + name + "' requires argument: [--n N]")
        return SHAPE(n)
    else:
        return SHAPE()


PARSER = argparse.ArgumentParser()
PARSER.add_argument("-s", "--shape", action='store', required=True,
                    help="Specify which face to display, options are " + str(tuple(s for s in mapping)))
PARSER.add_argument("-n", "--n", type=int, required=False, default=None,
                    help="additional argument to specify n, used for " + str(arg_n))

PARSER.add_argument("--legend", action='store_true', required=False,
                    help="put legend on plot")
PARSER.add_argument("--diameter", type=int, required=False, default=-1,
                    help="Specify diameter of search graph (longest possible sequence of faces on a geodesic)")

PARSER.add_argument("--no-show", action='store_true', required=False,
                    help="don't show the plot")

PARSER.add_argument("--save-file", action='store', required=False, default=None,
                    help="save initial image as specified file")
PARSER.add_argument("--width-display", type=float, required=False, default=None,
                    help="width of display in inches")
PARSER.add_argument("--height-display", type=float, required=False, default=None,
                    help="height of display in inches")

PARSER.add_argument("--face-name", action='store', required=False, default=None,
                    help="Specify face name if inputting a specific point")
PARSER.add_argument("--point-x", type=float, required=False, default=0.,
                    help="Specify point x if inputting a specific point (defaults to 0)")
PARSER.add_argument("--point-y", type=float, required=False, default=0.,
                    help="Specify point y if inputting a specific point (defaults to 0)")
