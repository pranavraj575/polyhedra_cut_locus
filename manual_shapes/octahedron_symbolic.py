import sympy as sym
import numpy as np
from utils.symbolic_utils import *

x = sym.Symbol('x')
y = sym.Symbol('y')

p = np.array([x, y])

# representing p on face 6 and inspecting face 0

p0 = (3*sym.sqrt(3), 1) + p
p1 = (2*sym.sqrt(3), 4) + sym_rotation_T(2*sym.pi/3)@p
p2 = (-2*sym.sqrt(3), 4) + sym_rotation_T(-2*sym.pi/3)@p
p3 = (-3*sym.sqrt(3), 1) + p
p4 = (-sym.sqrt(3), -5) + sym_rotation_T(2*sym.pi/3)@p
p5 = (sym.sqrt(3), -5) + sym_rotation_T(-2*sym.pi/3)@p

points=[p0,p1,p2,p3,p4,p5]

k=len(points)

L=[[None for _ in range(k)] for _ in range(k)]
for i in range(k):
    for j in range(k):
        if i!=j:
            L[i][j]=line_that_bisexts(points[i],points[j])
print(L[0][1])