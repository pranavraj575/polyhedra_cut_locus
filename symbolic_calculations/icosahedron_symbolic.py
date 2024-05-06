import sympy as sym
import numpy as np
from utils.symbolic_utils import *

x = sym.Symbol('p_1')
y = sym.Symbol('p_2')

p = np.array([x, y])

# representing p on face 17 and inspecting face 0

p0 = (5*sym.sqrt(3), 1) + p
p1 = (3*sym.sqrt(3), 7) + sym_rotation_T(2*sym.pi/3)@p
p2 = (1*sym.sqrt(3), 9) + sym_rotation_T(sym.pi)@p
p3 = (-1*sym.sqrt(3), 9) + sym_rotation_T(sym.pi)@p
p4 = (-3*sym.sqrt(3), 7) + sym_rotation_T(-2*sym.pi/3)@p
p5 = (-5*sym.sqrt(3), 1) + p
p6 = (-5*sym.sqrt(3), -3) + sym_rotation_T(sym.pi/3)@p
p7 = (-4*sym.sqrt(3), -6) + sym_rotation_T(sym.pi/3)@p
p8 = (-2*sym.sqrt(3), -8) + sym_rotation_T(2*sym.pi/3)@p
p9 = (2*sym.sqrt(3), -8) + sym_rotation_T(-2*sym.pi/3)@p
p10 = (4*sym.sqrt(3), -6) + sym_rotation_T(-sym.pi/3)@p
p11 = (5*sym.sqrt(3), -3) + sym_rotation_T(-sym.pi/3)@p

points = [p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11]
N = len(points)

L = get_all_bisecting_lines(points)
X = get_all_triple_points(points)

print(exp_simplify(points))
for p in exp_simplify(points):
    print(latexify(p))
print()


STARTS = [[['-' if j >= i else exp_simplify(item[0][dim]) for (j, item) in enumerate(row)]
           for (i, row) in enumerate(L)] for dim in range(2)]
VECTORS = [[['-' if j >= i else exp_simplify(item[1][dim]) for (j, item) in enumerate(row)]
            for (i, row) in enumerate(L)] for dim in range(2)]
for mats in STARTS, VECTORS:
    for dim in range(2):
        for i in range(N):
            mats[dim][i][i] = 0


print('x dimension intial point')
print(matrixify(STARTS[0]))

print('y dimension intial point')
print(matrixify(STARTS[1]))

print('x dimension direction')
print(matrixify(VECTORS[0]))
print('y dimension direction')
print(matrixify(VECTORS[1]))

print()
