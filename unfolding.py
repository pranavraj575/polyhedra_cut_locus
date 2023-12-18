from src.shapes import Shape
from matplotlib import pyplot as plt
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
ax= plt.gca()
#shape=Octahedron()
#shape.plot_unwrapping(np.zeros((2,1))+.1,1,7,None,ax)
shape=Tetrahedron()
shape.plot_unwrapping(np.array([[-.4], [-.6]]),3,0,None,ax)
plt.show()
