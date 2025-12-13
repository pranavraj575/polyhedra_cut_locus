import numpy as np


class Bound:
    def __init__(self, m, b, s, T, si, dimension=None, name=None, identifier=''):
        """
        linear boundary of the form mx<=b
        s,T,si represent the linear transformation that puts a point in our face to the neighboring face
        represented as a shift, rotation matrix, and another shift
        for some point x on this face, T (x + s) + si = x' where x' is the neighbor face coordinates

        :param m: row vector of bound (mx<=b)
        :param b: scalar of bound (mx<=b)
        :param s: first translation to shift bound to origin (for a point x, x' on new face is si+T(x+s))
        :param T: rotation to shift bound to origin (for a point x, x' on new face is si+T(x+s))
        :param si: final translation to shift bound to origin (for a point x, x' on new face is si+T(x+s))
        :param dimension: dimension of bound, if none, it is set, if value inserted, it is checked
        :param identifier: identifier of bound, should mention face names
        :param name: identifier of bound, or the bound that it is paired with, specify if spawned by another bound
        """
        # rescale so that |m|=1
        # dividing both m and b by |m| yields this with the same bound (mx<=b iff mx/|m|<=b/|m|)
        scale = np.linalg.norm(m)
        self.m = m/scale
        self.b = b/scale
        self.s = s
        self.T = T
        self.si = si
        self.dimension = self.check_valid(dimension)
        if name is None:
            base_id = str(tuple(self.m.flatten())) + str(self.b) + str(tuple(self.s.flatten())) + str(
                tuple(self.T.flatten())) + str(tuple(self.si.flatten())) + str(self.dimension)
            name = base_id + identifier
        self.name = name

    def __str__(self):
        return f'Bound(m={self.m.flatten()}, b={self.b}, s={self.s.flatten()}, si={self.si.flatten()}, T={self.T})'

    def check_valid(self, dimension):
        """
        Checks if self is valid, returns the correct dimension
        :param dimension: proposed dimension to check
        :return: dimension, raises exception if invalid
        """
        for twodcheck, nm in ((self.m, 'm'),
                              (self.s, 's'),
                              (self.si, 'si'),
                              (self.T, 'T'),
                              ):
            if len(np.shape(twodcheck)) != 2:
                raise Exception("bound should be a 2d array", nm, ':', twodcheck)

        if (np.shape(self.m)[0] != 1 or
                np.shape(self.s)[1] != 1 or
                np.shape(self.si)[1] != 1
        ):
            raise Exception("tried putting in bounds of incorrect dimensions")
        if (np.shape(self.T)[0] != np.shape(self.T)[1] or
                np.shape(self.m)[1] != np.shape(self.s)[0] or
                np.shape(self.s)[0] != np.shape(self.si)[0] or
                np.shape(self.si)[0] != np.shape(self.m)[1] or
                np.shape(self.m)[1] != np.shape(self.T)[0]
        ):
            raise Exception("tried putting in bounds of inconsistent dimensions")
        if dimension is None:
            dimension = np.shape(self.T)[0]
        else:
            if dimension != np.shape(self.T)[0]:
                raise Exception("inconsistent dimensions")
        if not self.within(np.zeros(dimension)):
            raise Exception("zero point needs to be within the face")
        return dimension

    def within(self, p, tol=0.):
        """
        checks if point p is within bound with some tolerance
        mp<=b+tol
        :param p: column vector (np array of dimension (self.dimension,1))
        :param tol: scalar
        :return: boolean
        """
        return np.dot(self.m, p) <= self.b + tol

    def grab_intersection(self, p, v):
        """
        find where p+vt intersects the bound (mx<=b)

        :param p: column vector (np array of dimension (self.dimension,1))
        :param v: column vector (np array of dimension (self.dimension,1))
        :return: point where p+vt intersects bound, (np array of dimension (self.dimension,1))
        """
        # m(p+vt)=b
        # simplify to mvt=b-mp
        # then t=(b-mp)/(mv)
        t = (self.b - np.dot(self.m, p))/np.dot(self.m, v)
        # note that v cannot be parallel to m, which makes sense for this to even work
        return p + v*t

    def shift_point(self, x):
        """
        shifts point x to equivalent x' on face according to bound
        :param x: column vector (np array of dimension (self.dimension,1))
        :return: column vector (np array of dimension (self.dimension,1))
        """
        return self.T@(x + self.s) + self.si

    def shift_vec(self, v):
        """
        shifts vector v to equivalent v'  according to bound
        :param v: column vector (np array of dimension (self.dimension,1))
        :return: column vector (np array of dimension (self.dimension,1))
        """
        return self.T@v

    def get_inverse_bound(self):
        """
        :return: inverse bound, from neighboring face to this face
        """

        Ti = np.linalg.inv(self.T)
        m = -self.m@Ti
        b = -self.b - np.dot(self.m, self.s) - np.dot(self.m@Ti, self.si)
        b = b.flatten()[0]
        return Bound(m, b, -self.si, Ti, -self.s, self.dimension, name=self.name)

    def concatenate_with(self, T=None, s=None):
        """
        given some point x that is translated previously like T' x + s',
            concatenate the translation to this one

            T((T' x + s')+s)+si=T T' x + T(s'+s)+si

        :param T: (self.dimension x self.dimension) translation matrix
        :param s: shift
        :return: (self.dimension x self.dimension) tranlstion, (self.dimension x 1) shift
        """
        if T is None: T = np.identity(self.dimension)
        if s is None: s = np.zeros((self.dimension, 1))
        return (self.T@T, self.si + self.T@(self.s + s))
