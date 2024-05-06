from utils.symbolic_utils import *

x = sym.Symbol('p_1')
y = sym.Symbol('p_2')

p = np.array([x, y])

# representing p on face 3 and inspecting face 0

p0 = (2*sym.sqrt(3), 0) + sym_rotation_T(sym.pi)@p
p1 = (2*sym.sqrt(3), 4) + p
p2 = (-2*sym.sqrt(3), 0) + sym_rotation_T(sym.pi)@p
p3 = (0, -2) + p

points = [p0, p1, p2, p3]
N = len(points)

L = get_all_bisecting_lines(points)

X = get_all_triple_points(points)
print(exp_simplify(points))

for i in range(N):
    j = (i + 1)%N
    print('line L' + str(i) + str(j))
    start, vector = exp_simplify(L[i][j])
    print(start, '+', str(vector) + 't')
    print()

print()
right_face = (0, 2), sym_rotation_T(-sym.pi/3)@(1, 0)
print('calculation of x^{01}')
x01 = line_intersection(right_face, L[0][1])
print(exp_simplify(x01))
print(latexify(exp_simplify(x01)))
print()

for a, b, c in combinations([0, 1, 2, 3], 3):
    abc = sorted([a, b, c])
    tup = str((a, b, c))[1:-1]
    print()
    print('equation of x^{' + tup + '}')
    print(exp_simplify(X[a][b][c]))
    print('\\item $x^{\{' + tup + '\}}(p)=\n\t' + latexify(exp_simplify(X[a][b][c])) + '$')
