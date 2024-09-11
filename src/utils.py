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





def within_lim(value, lim):
    return (value <= lim[1] and
            value >= lim[0])

def within_bounds(point, xlim, ylim):
    return (within_lim(point[0], xlim) and
            within_lim(point[1], ylim))

def get_correct_end_points(start, end, xlim, ylim):
    end = get_correct_exit_point(start, end, xlim, ylim)
    if end is None:
        return None, None
    start = get_correct_exit_point(end, start, xlim, ylim)
    if start is None:
        return None, None
    return start, end

def get_correct_exit_point(start, end, xlim, ylim):
    vec = end - start
    if vec[0] > 0:
        if (start + vec)[0] > xlim[1]:
            if start[0] > xlim[1]:
                return None
            vec = vec*((xlim[1] - start[0])/vec[0])
    elif vec[0] < 0:
        if (start + vec)[0] < xlim[0]:
            if start[0] < xlim[0]:
                return None
            vec = vec*((start[0] - xlim[0])/(-vec[0]))
    elif vec[0] == 0:
        if not within_lim(start[0], xlim):
            return None

    if vec[1] > 0:
        if (start + vec)[1] > ylim[1]:
            if start[1] > ylim[1]:
                return None
            vec = vec*((ylim[1] - start[1])/vec[1])
    elif vec[1] < 0:
        if (start + vec)[1] < ylim[0]:
            if start[1] < ylim[0]:
                return None
            vec = vec*((start[1] - ylim[0])/(-vec[1]))
    elif vec[1] == 0:
        if not within_lim(start[1], ylim):
            return None
    return vec + start