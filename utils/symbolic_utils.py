import sympy as sym
import numpy as np


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
    return avg, tangerine


def line_intersection(l1, l2):
    """
    returns interseciton point of l1 and l2
        lines are (point, direction vector)
    l1: p1 + v1*t1
    l2: p2 + v2*t2

    p1 + v1*t1 = p2 + v2*t2
    v1*t1 - v2*t2 = p2 - p1

    [v1 v2][t1, t2]^T = p2 - p1

    [t1, t2]^T = [v1 v2]^{-1} (p2 - p1)
    """
    p1, v1 = l1
    p2, v2 = l2

    V = np.array([[v1[0], v2[0]], [v1[1], v2[1]]])
    T = sym_2x2_inv(V)@(p2 - p1)

    t1 = T[0]

    return p1 + v1*t1
