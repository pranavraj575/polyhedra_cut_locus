import sympy as sym
import numpy as np
from utils.symbolic_utils import *

x = sym.Symbol('x')
y = sym.Symbol('y')

p = np.array([x, y])

# representing p on face 3 and inspecting face 0

p0 = (2*sym.sqrt(3), 0) + sym_rotation_T(sym.pi)@p
p1 = (2*sym.sqrt(3), 4) + p
p2 = (-2*sym.sqrt(3), 0) + sym_rotation_T(sym.pi)@p
p3 = (0, -2) + p

l01 = line_that_bisexts(p0, p1)
l12 = line_that_bisexts(p1, p2)
x012 = line_intersection(l01, l12)
print(sym.simplify(x012))
