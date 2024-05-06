import sympy as sym
import numpy as np
from itertools import combinations
from collections.abc import Iterable


def sym_rotation_T(theta):
    """
    :param theta: angle
    :return: 2x2 rotation matrix T by angle theta (np array)
    """
    # T where T @ v is v rotatied by theta
    return np.array([[sym.cos(theta), -sym.sin(theta)],
                     [sym.sin(theta), sym.cos(theta)]])


def sym_2x2_inv(arr):
    """
    inverts a symbolic 2x2 array
    """
    return np.array([[arr[1, 1], -arr[0, 1]], [-arr[1, 0], arr[0, 0]]])/(arr[0, 0]*arr[1, 1] - arr[0, 1]*arr[1, 0])


def line_that_bisexts(p, q):
    """
    returns line between p and q
    returns the starting point and the direction vector
    """
    avg = (p + q)/2
    normal = q - avg
    tangerine = sym_rotation_T(sym.pi/2)@normal

    return avg, tangerine  # /sym.sqrt(tangerine[0]**2 + tangerine[1]**2)


def line_intersection(l1, l2):
    """
    returns interseciton point of l1 and l2
        lines are (point, direction vector)
    l1: p1 + v1*t1
    l2: p2 + v2*t2

    p1 + v1*t1 = p2 + v2*t2
    v1*t1 - v2*t2 = p2 - p1

    [v1 -v2][t1, t2]^T = p2 - p1

    [t1, t2]^T = [v1, -v2]^{-1} (p2 - p1)
    """
    p1, v1 = l1
    p2, v2 = l2

    V = np.array([[v1[0], -v2[0]], [v1[1], -v2[1]]])
    T = sym_2x2_inv(V)@(p2 - p1)

    t1 = T[0]
    t2 = T[1]

    return p1 + v1*t1


def get_all_bisecting_lines(points):
    """
    returns an NxN array of lines that bisect the points
    arr[i][j] is the line bisecting points[i] and points[j]
    arr[i][i] is None
    """
    N = len(points)
    L = [[None for _ in range(N)] for _ in range(N)]

    for i in range(N):
        for j in range(N):
            if i != j:
                L[i][j] = line_that_bisexts(points[i], points[j])
    return L


def get_all_triple_points(points):
    """
    gets all x^{i,j,k} where x^{i,j,k} is equidistant from points i, j, and k
    returns an NxNxN array where arr[i][j][k] is the triple point of points i, j, and k
    arr[i,j,k] is None if any of the indices are equal
    """
    N = len(points)
    L = get_all_bisecting_lines(points)
    X = [[[None for _ in range(N)] for _ in range(N)] for _ in range(N)]
    for i in range(N):
        for j in range(N):
            for k in range(N):
                if not any(a == b for (a, b) in combinations([i, j, k], 2)):
                    X[i][j][k] = line_intersection(L[i][j], L[j][k])
    return X


def latexify(equation):
    """
    returns latex friendly equation
    """
    # if this is a tuple, latexify each element and tuple them together
    if (not type(equation) == str) and isinstance(equation, Iterable):
        return '\\left(' + (', '.join([latexify(thing) for thing in equation])) + '\\right)'

    if not type(equation) == str:
        equation = str(equation)

    def first_term(equation: str):
        paren_count = 0
        i = 1
        next_space = False
        if equation[0] == '(':
            paren_count += 1
        else:
            next_space = True
        while i < len(equation):
            char = equation[i]
            if char == '(':
                paren_count += 1
            elif char == ')':
                paren_count -= 1
            if paren_count == 0:
                if next_space:
                    if char == ' ':
                        break
                else:
                    i += 1
                    break
            i += 1
        return equation[:i]

    def last_term(equation: str):
        equationp = equation[::-1]
        replacement_temp = ''
        while replacement_temp in equationp:
            replacement_temp += str(np.random.random())
        equationp = equationp.replace('(', replacement_temp)
        equationp = equationp.replace(')', '(')
        equationp = equationp.replace(replacement_temp, ')')
        last = first_term(equationp)
        last = last.replace('(', replacement_temp)
        last = last.replace(')', '(')
        last = last.replace(replacement_temp, ')')
        return last[::-1]

    def sqrtify(equation: str):
        while equation.count('sqrt') > equation.count('\sqrt'):
            i = 0
            # we want the first occurance of sqrt in equation[i:] to be unlatexified
            while equation.find('sqrt', i) == equation.find('\\sqrt', i) + 1:
                if equation.find('\\sqrt', i) < 0:
                    break
                i = equation.find('sqrt', i) + 4
            i = equation.find('sqrt', i) + 4
            # now equation[i:] is (...
            # the first non replaced sqrt(...
            first = first_term(equation[i:])
            end = equation[i + len(first):]
            # first should start with 'sqrt(' and end with ')'
            equation = equation[:i - 4] + '\\sqrt{' + strip(first) + '}' + end
        return equation

    def fracify(equation):
        """
        make all fractions a/b into \frac{a}{b}
        """
        while '/' in equation:
            i = equation.find('/')
            nom = last_term(equation[:i])
            start = equation[:i - len(nom)]
            denom = first_term(equation[i + 1:])
            end = equation[i + 1 + len(denom):]

            equation = start + '\\frac{' + strip(nom) + '}{' + strip(denom) + '}' + end
        return equation

    def strip(equation: str):
        """
        strips outer parentheses
        """
        while equation.startswith('(') and equation.endswith(')'):
            equation = equation[1:-1]
        return equation

    def depower(equation: str):
        """
        deals with x**y
        """
        while '**' in equation:
            i = equation.find('**')
            pwr = first_term(equation[i + 2:])
            last = equation[i + 2 + len(pwr):]
            equation = equation[:i] + '^{' + strip(pwr) + '}' + last
        return equation

    equation = sqrtify(equation)
    equation = fracify(equation)
    equation = depower(equation)
    equation = strip(equation)
    return equation.replace('*', '')


def exp_simplify(expression):
    if isinstance(expression, Iterable):
        return [exp_simplify(thing) for thing in expression]
    if expression is None:
        return None
    return sym.simplify(sym.expand(sym.simplify(expression)))


def matrixify(arr):
    """
    returns latex version of an array
    arr must be 2d
    """
    rope = '\\begin{bmatrix}'
    for row in arr:
        rope += '\n\t'
        line = ' & '.join(['' if item is None else latexify(item) for item in row])
        rope += line
        rope += '\\\\'
    rope += '\n\\end{bmatrix}'
    return rope
