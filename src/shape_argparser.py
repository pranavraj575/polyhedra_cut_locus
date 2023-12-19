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

arg_n = ('prism', 'antiprism', 'pyramid', 'longpyramid', 'bipyramid', 'longbipyramid', 'mirror')
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
