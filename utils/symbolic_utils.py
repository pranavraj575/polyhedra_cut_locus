import sympy as sym
import numpy as np
from itertools import combinations
from collections.abc import Iterable
from itertools import chain


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

    lparens = '({'
    rparens = ')}'

    if not type(equation) == str:
        equation = str(equation)

    def first_term(equation: str):
        divisions = ' *+/'
        paren_count = 0
        i = 0
        while i < len(equation):
            char = equation[i]
            if char in lparens:
                paren_count += 1
            elif char in rparens:
                paren_count -= 1
            if paren_count < 0:
                break
            elif paren_count == 0 and char in divisions:
                break

            i += 1
        return equation[:i]

    def swap(equation: str, a, b):
        equationp = ''
        for char in equation:
            if char == a:
                equationp += b
            elif char == b:
                equationp += a
            else:
                equationp += char
        return equationp

    def last_term(equation: str):
        equationp = equation[::-1]
        for left, right in zip(lparens, rparens):
            equationp = swap(equationp, left, right)
        last = first_term(equationp)
        for left, right in zip(lparens, rparens):
            last = swap(last, left, right)
        return last[::-1]

    def sqrtify(equation: str):
        while equation.count('sqrt') > equation.count('\\sqrt'):
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
            if first_term(equation) == equation:
                equation = equation[1:-1]
            else:
                break
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
    equation = depower(equation)
    equation = fracify(equation)
    equation = strip(equation)
    # no need for cdot
    return equation.replace('*', ' ')
    #return equation.replace('*', '\cdot ')


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


def equality(a, b):
    an, ad = sym.fraction(exp_simplify(a))
    bn, bd = sym.fraction(exp_simplify(b))
    return sym.Eq(an*bd, bn*ad)


def eq_points_equal(p1, p2):
    return (equality(p1[dim], p2[dim]) for dim in range(len(p1)))


def eq_all_points_equal(points):
    """
    returns iterable of equations that is where all the points have equality
    p_i=p_{i+1} is sufficient
    """
    return chain(*(eq_points_equal(points[i], points[(i + 1)%len(points)]) for i in range(len(points))))


def roots(equation: sym.core.Expr, var_to_use=None):
    def coefs(eq: sym.core.Expr, var_idx):
        terms, vars = eq.as_terms()
        var = vars[var_idx]

        coeff_list = [0 for _ in range(1 + sym.polys.polytools.degree(eq, var))]

        for term, (_, degrees, _) in terms:
            degree = degrees[var_idx]
            coeff_list[degree] += term/(var**degree)

        return coeff_list

    def linear_roots(eq: sym.core.Expr, var_idx):
        c, b = coefs(eq, var_idx)  # should only be called on linear eqs
        return (-c/b,)

    def quad_roots(eq: sym.core.Expr, var_idx):
        c, b, a = coefs(eq, var_idx)
        thing = b**2 - 4*a*c
        if thing == 0:  # one soln
            return (-b/(2*a),)
        return ((-b + sym.sqrt(thing))/(2*a), (-b - sym.sqrt(thing))/(2*a))

    def cubic_roots(eq: sym.core.Expr, var_idx):
        d, c, b, a = coefs(eq, var_idx)

        p = (3*a*c - b**2)/(3*a**2)
        q = (2*b**3 - 9*a*b*c + 27*a**2*d)/(27*a**3)
        delta = (q**2/4 + p**3/27)

        u = sym.root((-q/2 + sym.sqrt(delta)), 3)
        v = sym.root((-q/2 - sym.sqrt(delta)), 3)
        x1 = u + v - b/(3*a)

        # these could be complex
        # if delta = 0, these are the same
        x2 = -(u + v)/2 - b/(3*a) + (u - v)*sym.sqrt(3)/(2*sym.I)
        x3 = -(u + v)/2 - b/(3*a) - (u - v)*sym.sqrt(3)/(2*sym.I)
        return (x1, x2, x3)

    def poly_roots(eq: sym.core.Expr, var_idx):
        terms, vars = eq.as_terms()
        var = vars[var_idx]
        degree = sym.polys.polytools.degree(eq, var)
        if degree == 0:
            return None  # this variable is not in equation
        elif degree == 1:
            return linear_roots(eq, var_idx)
        elif degree == 2:
            return quad_roots(eq, var_idx)
        elif degree == 3:
            return cubic_roots(eq, var_idx)
        else:
            raise Exception('degree ' + str(degree) + ' equation appeared')

    num, denom = sym.fraction(exp_simplify(equation))

    constant, factor_list = sym.factor_list(exp_simplify(num))
    if len(factor_list) > 1:
        for fac, multiple in factor_list:
            for soln in roots(fac, var_to_use=var_to_use):
                yield soln
        return

    eq = sym.expand(num)
    terms, vars = eq.as_terms()
    if var_to_use is None or (var_to_use not in vars):
        var_to_use = vars[0]
    var_idx = vars.index(var_to_use)
    results = poly_roots(eq, var_idx=var_idx)
    if results is None:
        return [{}]  # any values work
    # degree = sym.polys.polytools.degree(equation, var)
    for result in results:
        yield {var_to_use: result}


def simul_roots_solver(L):
    tree = []
    for equation in L:
        tree.append(list(roots(equation)))

    def get_vars(thing):
        try:
            return thing.as_terms()[1]
        except:
            return []

    # now we must see where these intersect
    def intersect_solver(sol1, sol2):
        for dic1 in sol1:
            for dic2 in sol2:
                variable1, value1 = list(dic1.items())[0]
                variable2, value2 = list(dic2.items())[0]

                if variable1 == variable2:
                    inter = value1 - value2
                    if sym.simplify(inter) == 0:
                        yield {variable1: value1}
                    else:
                        for intersect in sym.solve(inter):
                            yield {variable1: intersect}
                else:
                    vars1, vars2 = get_vars(value1), get_vars(value2)
                    if vars1 is []:
                        yield {variable1: value1, variable2: value2.subs({variable1: value1})}
                    if vars2 is []:
                        yield {variable2: value2, variable1: value1.subs({variable2: value2})}
    while len(tree) > 1:
        sol1, sol2 = tree.pop(), tree.pop()
        tree.append(list(intersect_solver(sol1, sol2)))
    return tree[0]


