# To implement new shapes

We construct a polyhedron as a union of 2-d faces (the Face class is implemented in `src/face.py`, and is defined by its linear boundaries).

To construct a polyhedron, first initialize all faces (i.e. for a cube, start with 6 empty faces).
Then, consider every pair of faces that share an edge, and define this boundary using the Face.add_boundary_paired method.

A boundary is defined in `src/bound.py` and is the transformation that takes a point on an edge of one face and returns the coordinates of the same point on a neighboring face.
This consists of a translation (from the edge towards the origin), a rotation, then a translation to the edge of the neighboring face.

The file `src/shape_creation.py` has all currently implemented shapes.
Any new implementations should probably inherit the ConvexPolyhderon class (or at least the Shape class).
