import numpy as np

# rotation matrices
def rotation_T(theta):
    """
    :param theta: angle
    :return: 2x2 rotation matrix T by angle theta (np array)
    """
    # T where T @ v is v rotatied by theta
    return np.array([[np.cos(theta), -np.sin(theta)],
                     [np.sin(theta), np.cos(theta)]])


def coltation(theta):
    """
    returns 2d unit column vector towards theta
    :param theta: angle
    :return: <<cos theta>,<sin theta>>
    """
    return rotation_T(theta)[:, [0]]


def rowtation(theta):
    """
    returns 2d unit row vector towards theta
    :param theta: angle
    :return: <<cos theta,sin theta>>
    """
    return coltation(theta).T


def flatten(L):
    """
    flattens a list of lists
    :param L: list of lists
    :return: L flattened
    """
    out = []
    for l in L:
        out += l
    return out
