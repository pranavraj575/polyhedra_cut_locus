from utils.symbolic_utils import *

x = sym.Symbol('p_1')
y = sym.Symbol('p_2')

p = np.array([x, y])

# representing p on face 6 and inspecting face 0

p0 = (3*sym.sqrt(3), 1) + p
p1 = (2*sym.sqrt(3), 4) + sym_rotation_T(2*sym.pi/3)@p
p2 = (-2*sym.sqrt(3), 4) + sym_rotation_T(-2*sym.pi/3)@p
p3 = (-3*sym.sqrt(3), 1) + p
p4 = (-sym.sqrt(3), -5) + sym_rotation_T(2*sym.pi/3)@p
p5 = (sym.sqrt(3), -5) + sym_rotation_T(-2*sym.pi/3)@p

right_face = (0, 2), sym_rotation_T(-sym.pi/3)@(1, 0)
left_face = (0, 2), sym_rotation_T(-2*sym.pi/3)@(1, 0)
bottom_face = (0, -1), sym_rotation_T(0)@(1, 0)

points = [p0, p1, p2, p3, p4, p5]
N = len(points)

L = get_all_bisecting_lines(points)
X = get_all_triple_points(points)
print(exp_simplify(points))
for p in exp_simplify(points):
    print(latexify(p))
print()

for i in range(N):
    j = (i + 1)%N
    print('line L' + str(i) + str(j))
    start, vector = exp_simplify(L[i][j])
    print(latexify(start), '+', latexify(vector) + 't')
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
print()

print('x dimension direction')
print(matrixify(VECTORS[0]))
print('y dimension direction')
print(matrixify(VECTORS[1]))
print()

print('calculation of x^{01}')
x01 = line_intersection(right_face, L[0][1])
print(exp_simplify(x01))
print(latexify(exp_simplify(x01)))
print()

print('calculation of x^{23}')
x23 = line_intersection(left_face, L[2][3])
print(exp_simplify(x23))
print(latexify(exp_simplify(x23)))
print()

print('calculation of x^{45}')
x45 = line_intersection(bottom_face, L[4][5])
print(exp_simplify(x45))
print(latexify(exp_simplify(x45)))
print()

for a, b, c in (list(combinations([0, 1, 2, 3], 3)) +
                list(combinations([0, 1, 4, 5], 3)) +
                list(combinations([2, 3, 4, 5], 3))):
    abc = sorted([a, b, c])
    tup = str((a, b, c))[1:-1]
    print()
    print('equation of x^{' + tup + '}')
    print(exp_simplify(X[a][b][c]))
    print('\\item $x^{\{' + tup + '\}}(p)=\n\t' + latexify(exp_simplify(X[a][b][c])) + '$')
