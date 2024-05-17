import sympy as sym
from itertools import chain

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
# transpose VECTORS because we want to do L[i][j] for i<j and also have a lower triangular matrix
VECTORS = [[['-' if j <= i else exp_simplify(item[1][dim]) for (j, item) in enumerate(row)]
            for (i, row) in enumerate(L)] for dim in range(2)]
VECTORS = [[[VECTORS[dim][j][i] for j in range(N)] for i in range(N)] for dim in range(2)]

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

print('calculation of x^{12}')
x12 = line_intersection(right_face, L[1][2])
print(exp_simplify(x12))
print(latexify(exp_simplify(x12)))
print()

print('calculation of x^{23}')
x23 = line_intersection(left_face, L[2][3])
print(exp_simplify(x23))
print(latexify(exp_simplify(x23)))
print()

print('calculation of x^{34}')
x34 = line_intersection(left_face, L[3][4])
print(exp_simplify(x34))
print(latexify(exp_simplify(x34)))
print()

print('calculation of x^{45}')
x45 = line_intersection(bottom_face, L[4][5])
print(exp_simplify(x45))
print(latexify(exp_simplify(x45)))
print()

print('calculation of x^{05}')
x05 = line_intersection(right_face, L[0][5])
print(exp_simplify(x05))
print(latexify(exp_simplify(x05)))
print()

for a, b, c in combinations(range(6), 3):
    abc = sorted([a, b, c])
    tup = str((a, b, c))[1:-1]
    print()
    print('equation of x^{' + tup + '}')
    print(exp_simplify(X[a][b][c]))
    print(latexify(exp_simplify(X[a][b][c])))
    # print('\\item $x^{\{' + tup + '\}}(p)=\n\t' + latexify(exp_simplify(X[a][b][c])) + '$')

print('sympy equality solver does not function well, use mathematica or wolfram alpha beyond this point')

print('inspecting line l12')
line = (1, 2)
# pt=x12

print('boundaries:')

# equality calculations from wolfram alpha
# equality_calculations[(i,j)], correspoinding with line Lij (i<j) is a dictionary D such that
# D[(a,b)] (a<b) is the set where equality of x^{a,i,j} = x^{b,i,j} holds
equality_calculations = {
    (1, 2):
        {
            (0, 3): 'x=0; (y-4)^2+x^2=12',
        }
}
for a, b in combinations(range(6), 2):
    if a not in line and b not in line:
        print('find where ' +
              'x^{' + str(a) + str(line[0]) + str(line[1]) + '}' +
              ' and ' + 'x^{' + str(b) + str(line[0]) + str(line[1]) + '}' +
              ' are equal using wolfram alpha or other software')
        for dim in range(2):
            eq = equality(X[a][line[0]][line[1]][dim], X[b][line[0]][line[1]][dim])
            print(eq.lhs, '=', eq.rhs)
        print()
print()

quit()
print()
print('equations where all x^{a,b,c} are equal:')
print(sym.solve(eq_all_points_equal([X[a][b][c] for a, b, c in combinations(range(6), 3)])))
print()

print()
print('equations where x^{0,1,2}=x^{0,1,3}=x^{0,2,3}=x^{1,2,3}')
print(sym.solve(eq_all_points_equal([X[a][b][c] for a, b, c in combinations(range(4), 3)])))
print()
print()

print()
print('equations where x^{0,1,2}=x^{0,1,3}=x^{0,2,3}=x^{1,2,3} and x^{0,3,4}=x^{0,3,5}=x^{0,4,5}=x^{3,4,5}')
eq0123 = eq_all_points_equal([X[a][b][c] for a, b, c in combinations(range(4), 3)])
eq0345 = eq_all_points_equal([X[a][b][c] for a, b, c in combinations([0, 3, 4, 5], 3)])
print(sym.solve(chain(eq0123, eq0345)))
print('equations those equalties hold and x^{0,1,2} is above x^{0,3,4}')
print(sym.solve(
    chain((eq.subs({x: 0}) for eq in eq0123),
          (eq.subs({x: 0}) for eq in eq0345),
          (sym.GreaterThan(X[0][1][2][1], X[3][4][5][1]).subs({x: 0}),)
          )
)
)
print()

line = [1, 2]
pt = x12
subject = 0
print('top point')
print(latexify(exp_simplify(np.dot(X[subject][line[0]][line[1]] - pt, L[line[0]][line[1]][1]))))
print('other points')
for other in range(6):
    if other not in [subject] + line:
        print(latexify(exp_simplify(np.dot(X[other][line[0]][line[1]] - pt, L[line[0]][line[1]][1]))))
