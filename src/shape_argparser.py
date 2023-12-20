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


arg_n = ('prism', 'antiprism', 'pyramid', 'longpyramid', 'bipyramid', 'longbipyramid', 'mirror')
# shapes that take n as an arg

PARSER = argparse.ArgumentParser()
PARSER.add_argument("-s", "--shape", action='store', required=True,
                    help="Specify which face to display, options are " + str(tuple(s for s in mapping)))
PARSER.add_argument("--n", type=int, required=False, default=3,
                    help="additional argument to specify n, used for " + str(arg_n))

PARSER.add_argument("--legend", action='store_true', required=False,
                    help="put legend on plot")
PARSER.add_argument("--diameter", type=int, required=False, default=-1,
                    help="Specify diameter of search graph (longest possible sequence of faces on a geodesic)")

PARSER.add_argument("--face-name", action='store', required=False, default=None,
                    help="Specify face name if inputting a specific point")
PARSER.add_argument("--point-x", type=float, required=False, default=0.,
                    help="Specify point x if inputting a specific point (defaults to 0)")
PARSER.add_argument("--point-y", type=float, required=False, default=0.,
                    help="Specify point y if inputting a specific point (defaults to 0)")
