'''
Created on Mar 2, 2016

@author: Francisco Dominguez
'''
from __future__ import division
import numpy as np


class Algebra:
    def matrix(self):
        pass

    def vector(self):
        pass

    def exp(self):
        pass

    def __add__(self, other):
        pass


class Group:
    def matrix(self):
        pass

    def log(self):
        pass

    def __mul__(self, other):
        pass


class SO2(Group):
    """2D rotation Lie group"""

    def __init__(self, matrix):
        """Group element constructor
        :param matrix: 2D rotation matrix (2x2)
        :type matrix: numpy.ndarray
        """
        self.M = matrix

    def __mul__(self, other):
        """Returns the group element for the equivalent rotation as applying the given ones one after the other
        :param other: element to multiply this object by
        :type other: SO2
        :return: equivalent group element
        :rtype: SO2
        """
        return SO2(self.M.dot(other.matrix()))

    def matrix(self):
        """
        :return: this group element represented as matrix
        :rtype: numpy.ndarray (2x2)
        """
        return self.M

    def log(self):
        """
        :return: algebra element associated with this group element
        :rtype: so2
        """
        return so2(theta=np.arctan2(self.M[1, 0], self.M[0, 0]))

    def __J(self):
        return so2.G * self.M


class so2(Algebra):
    """2D rotation Lie algebra"""
    G = np.matrix([[0, -1],
                   [1, 0]])

    def __init__(self, **kwargs):
        """Algebra element constructor
        :param kwargs: must be `theta = float`, `matrix = ndarray` with dimensions 2x2 or `vector = ndarray` with dimensions 1x1
        :type kwargs
        :raises TypeError, if the needed argument (theta, matrix or vector) is not given
        """
        if "theta" in kwargs:
            self.angle = kwargs["theta"]
        elif "matrix" in kwargs:
            self.angle = kwargs["matrix"][1, 0]
        elif "vector" in kwargs:
            self.angle = kwargs["vector"][0]
        else:
            raise TypeError("Argument must be theta, matrix or vector")

        if self.angle >= np.pi or self.angle <= -np.pi:
            self.angle = np.arctan2(np.sin(self.angle), np.cos(self.angle))

    def __add__(self, other):
        """Returns the algebra element for the equivalent rotation as applying the given ones one after the other
        :param other: element to add this element to
        :type other: so2
        :return: equivalent algebra element
        :rtype: so2
        """
        R1 = self.exp()
        R2 = other.exp()
        R = R1 * R2
        return R.log()

    def exp(self):
        """
        :return: group element associated with this algebra element
        :rtype: SO2
        """
        theta = self.angle
        cs = np.cos(theta)
        sn = np.sin(theta)
        R = np.matrix([[cs, -sn],
                       [sn, cs]])
        return SO2(matrix=R)

    def matrix(self):
        """
        :return: this algebra element represented as skew matrix (2x2)
        :rtype: ndarray
        """
        return self.angle * so2.G

    def vector(self):
        """
        :return: this algebra element represented as vector (1x1)
        :rtype: ndarray
        """
        return np.array([self.angle])

    def theta(self):
        """
        :return: the rotation in radians counterclock-wise
        :rtype: float
        """
        return self.angle


# TODO SE2/se2: define if A*B/a+b should be A(B(point)) or B(A(point))
class SE2(Group):
    """2D transformation Lie group"""

    def __init__(self, matrix):
        """Group element constructor
        :param matrix: 2D transformation matrix (3x3)
        :type matrix: numpy.ndarray
        """
        self.M = matrix

    def matrix(self):
        """
        :return: this group element represented as matrix
        :rtype: numpy.ndarray (3x3)
        """
        return self.M

    def log(self):
        """
        :return: algebra element associated with this group element
        :rtype: se2
        """
        w = SO2(self.M[0:2, 0:2])
        t = self.M[0:2, 2]
        theta = w.log().theta()
        cs = np.cos(theta)
        sn = np.sin(theta)
        A = sn / theta if theta != 0 else 1
        B = (1 - cs) / theta if theta != 0 else 0
        V1 = 1 / (A ** 2 + B ** 2) * np.matrix([[A, B],
                                                [-B, A]])
        Vt = V1.dot(t)
        return se2(vector=np.array([theta, Vt[0, 0], Vt[0, 1]]))

    def __mul__(self, other):
        """Returns the group element for the equivalent transformation as applying the given ones one after the other
        :param other: element to multiply this object by
        :type other: SE2
        :return: equivalent group element
        :rtype: SE2
        """
        return SE2(self.M.dot(other.matrix()))


class se2(Algebra):
    """2D transformation Lie algebra"""
    G1 = np.matrix([[0, -1, 0],
                    [1, 0, 0],
                    [0, 0, 0]])
    G2 = np.matrix([[0, 0, 1],
                    [0, 0, 0],
                    [0, 0, 0]])
    G3 = np.matrix([[0, 0, 0],
                    [0, 0, 1],
                    [0, 0, 0]])

    def __init__(self, **kwargs):
        """Algebra element constructor
        :param kwargs: must be `matrix = ndarray` with dimensions 3x3 or `vector = ndarray` with dimensions 1x3
        :type kwargs
        :raises TypeError, if the needed argument (matrix or vector) is not given
        """
        if "matrix" in kwargs:
            self.w = np.array([0, 0, 0])
            self.w[0] = kwargs["matrix"][1, 0]
            self.w[1:3] = kwargs["matrix"][0:2, 2]
        elif "vector" in kwargs:
            self.w = kwargs["vector"]
        else:
            raise TypeError("Argument must be matrix or vector")

        if self.w[0] >= np.pi or self.w[0] <= -np.pi:
            self.w[0] = np.arctan2(np.sin(self.w[0]), np.cos(self.w[0]))

    def __add__(self, other):
        """Returns the algebra element for the equivalent transformation as applying the given ones one after the other
        :param other: element to add this element to
        :type other: se2
        :return: equivalent algebra element
        :rtype: se2
        """
        R1 = self.exp()
        R2 = other.exp()
        R = R1 * R2
        return R.log()

    def exp(self):
        """
        :return: group element associated with this algebra element
        :rtype: SE2
        """
        theta = self.w[0]
        cs = np.cos(theta)
        sn = np.sin(theta)
        w = so2(theta=theta)

        t = self.w[1:3]
        R = w.exp().matrix()
        T = np.eye(3, 3)
        T[0:2, 0:2] = R

        # TODO Donde esta definido?
        V = 1 / theta * np.matrix([[sn, -(1 - cs)],
                                   [1 - cs, sn]]) if theta != 0 else np.eye(2)
        tmp = V.dot(t)
        T[0:2, 2] = tmp
        return SE2(T)

    def matrix(self):
        """
        :return: this algebra element represented as matrix (3x3)
        :rtype: ndarray
        """
        return self.w[0] * se2.G1 + self.w[1] * se2.G2 + self.w[2] * se2.G3

    def vector(self):
        """
        :return: this algebra element represented as vector (1x3)
        :rtype: ndarray
        """
        return self.w


class SO3(Group):
    """3D rotation Lie group"""

    def __init__(self, matrix):
        """Group element constructor
        :param matrix: 3D rotation matrix (3x3)
        :type matrix: numpy.ndarray
        """
        self.M = matrix

    def log(self):
        """
        :return: algebra element associated with this group element
        :rtype: so3
        """
        error = 0.0000001
        if np.linalg.norm(np.eye(3) - self.M) < error:
            # case theta == 0
            return so3(vector=np.array([0, 0, 0]))
        elif np.abs(np.trace(self.M) + 1) < error:
            # case theta == pi
            if self.M[0, 0] > 0:
                w = 1 / np.sqrt(2 * (1 + self.M[0, 0])) * np.array([1 + self.M[0, 0], self.M[1, 0], self.M[2, 0]])  # wx
            elif self.M[1, 1] > error:
                w = 1 / np.sqrt(2 * (1 + self.M[1, 1])) * np.array([self.M[0, 1], 1 + self.M[1, 1], self.M[2, 1]])  # wy
            else:
                w = 1 / np.sqrt(2 * (1 + self.M[2, 2])) * np.array([self.M[0, 2], self.M[1, 2], 1 + self.M[2, 2]])  # wz
            w = np.pi * w
            return so3(vector=w)
        else:
            cs = (np.trace(self.M) - 1) / 2
            theta = np.arccos(cs)
            sn = np.sin(theta)
            logR = theta / (2 * sn) * (self.M - self.M.T)
            return so3(vector=np.array([logR[2, 1], logR[0, 2], logR[1, 0]]))

    def matrix(self):
        """
        :return: this group element represented as matrix
        :rtype: numpy.ndarray (3x3)
        """
        return self.M

    def __mul__(self, other):
        """Returns the group element for the equivalent rotation as applying the given ones one after the other
        :param other: element to multiply this object by
        :type other: SO3
        :return: equivalent group element
        :rtype: SO3
        """
        return SO3(self.M.dot(other.matrix()))


class so3(Algebra):
    """3D rotation Lie algebra"""
    G1 = np.matrix([[0, 0, 0],
                    [0, 0, -1],
                    [0, 1, 0]])
    G2 = np.matrix([[0, 0, 1],
                    [0, 0, 0],
                    [-1, 0, 0]])
    G3 = np.matrix([[0, -1, 0],
                    [1, 0, 0],
                    [0, 0, 0]])

    def __init__(self, **kwargs):
        """Algebra element constructor
        :param kwargs: must be `matrix = ndarray` with dimensions 3x3 or `vector = ndarray` with dimensions 1x3
        :type kwargs
        :raises TypeError, if the needed argument (matrix or vector) is not given
        """
        if "vector" in kwargs:
            self.w = kwargs["vector"]
        elif "matrix" in kwargs:
            m = kwargs["matrix"]
            self.w = np.array([0, 0, 0])
            self.w[0] = m[2, 1]
            self.w[1] = m[0, 2]
            self.w[2] = m[1, 0]
        else:
            raise TypeError("Argument must be matrix or vector")

    def __add__(self, other):
        """Returns the algebra element for the equivalent rotation as applying the given ones one after the other
        :param other: element to add this element to
        :type other: so3
        :return: equivalent algebra element
        :rtype: so3
        """
        R1 = self.exp()
        R2 = other.exp()
        R = R1 * R2
        return R.log()

    def magnitude(self):
        """Given the vector w=(wx,wy,wz) returns its norm, which is how many radians this rotation makes around the normalized vector
        :return: radians this rotation makes
        :rtype: float
        """
        return np.linalg.norm(self.w)

    def exp(self):
        """
        :return: group element associated with this algebra element
        :rtype: SO3
        """
        wx = self.matrix()
        theta = np.linalg.norm(self.w)
        cs = np.cos(theta)
        sn = np.sin(theta)
        I = np.eye(3)
        a = (sn / theta) if theta != 0 else 1
        b = ((1 - cs) / theta ** 2) if theta != 0 else 1 / 2.0
        R = I + a * wx + b * wx.dot(wx)
        return SO3(R)

    def vector(self):
        """
        :return: this algebra element represented as vector (1x3)
        :rtype: ndarray
        """
        return self.w

    def matrix(self):
        """
        :return: this algebra element represented as skew matrix (3x3)
        :rtype: ndarray
        """
        return so3.G1 * self.w[0] + so3.G2 * self.w[1] + so3.G3 * self.w[2]


class SE3(Group):
    """3D transformation Lie group"""

    def __init__(self, matrix):
        """Group element constructor
        :param matrix: 3D transformation matrix (4x4)
        :type matrix: numpy.ndarray
        """
        self.M = matrix

    def __mul__(self, other):
        """Returns the group element for the equivalent transformation as applying the given ones one after the other
        :param other: element to multiply this object by
        :type other: SE3
        :return: equivalent group element
        :rtype: SE3
        """
        return SE3(self.M.dot(other.matrix()))

    def matrix(self):
        """
        :return: this group element represented as matrix
        :rtype: numpy.ndarray (4x4)
        """
        return self.M

    def log(self):
        """
        :return: algebra element associated with this group element
        :rtype: se3
        """
        R = self.M[0:3, 0:3]
        w = SO3(R).log()
        t = self.M[0:3, 3]

        I = np.eye(3)
        theta = w.magnitude()
        wx = w.matrix()
        wx2 = wx.dot(wx)
        a = (1 / theta**2) * (1 - (theta * np.sin(theta)) / (2 * (1 - np.cos(theta)))) if theta != 0 else 1 / 12
        V = I - 1 / 2.0 * wx + a * wx2

        v = V.dot(t)

        return se3(vector=np.append(w.vector(), v))


class se3(Algebra):
    """3D transformation Lie algebra"""
    G1 = np.matrix([[0, 0, 0, 0],
                    [0, 0, -1, 0],
                    [0, 1, 0, 0],
                    [0, 0, 0, 0]])

    G2 = np.matrix([[0, 0, 1, 0],
                    [0, 0, 0, 0],
                    [-1, 0, 0, 0],
                    [0, 0, 0, 0]])

    G3 = np.matrix([[0, -1, 0, 0],
                    [1, 0, 0, 0],
                    [0, 0, 0, 0],
                    [0, 0, 0, 0]])

    G4 = np.matrix([[0, 0, 0, 1],
                    [0, 0, 0, 0],
                    [0, 0, 0, 0],
                    [0, 0, 0, 0]])

    G5 = np.matrix([[0, 0, 0, 0],
                    [0, 0, 0, 1],
                    [0, 0, 0, 0],
                    [0, 0, 0, 0]])

    G6 = np.matrix([[0, 0, 0, 0],
                    [0, 0, 0, 0],
                    [0, 0, 0, 1],
                    [0, 0, 0, 0]])

    def __init__(self, **kwargs):
        """Algebra element constructor
        :param kwargs: must be `matrix = ndarray` with dimensions 4x4 or `vector = ndarray` with dimensions 1x6
        :type kwargs
        :raises TypeError, if the needed argument (matrix or vector) is not given
        """
        if "vector" in kwargs:
            self.w = kwargs["vector"]
        elif "matrix" in kwargs:
            m = kwargs["matrix"]
            self.w = np.array([0, 0, 0, 0, 0, 0])
            self.w[0] = m[2, 1]
            self.w[1] = m[0, 2]
            self.w[2] = m[1, 0]
            self.w[3:6] = m[0:3, 3]
        else:
            raise TypeError("Argument must be matrix or vector")

    def __add__(self, other):
        """Returns the algebra element for the equivalent transformation as applying the given ones one after the other
        :param other: element to add this element to
        :type other: se3
        :return: equivalent algebra element
        :rtype: se3
        """
        R1 = self.exp()
        R2 = other.exp()
        R = R1 * R2
        return R.log()

    def vector(self):
        """
        :return: this algebra element represented as vector (1x6)
        :rtype: ndarray
        """
        return self.w

    def matrix(self):
        """
        :return: this algebra element represented as matrix (4x4)
        :rtype: ndarray
        """
        x = np.zeros((4, 4))
        for i in range(6):
            x += self.w[i] * eval("G" + str(i + 1))
        return x

    def exp(self):
        """
        :return: group element associated with this algebra element
        :rtype: SE3
        """
        w = so3(vector=self.w[0:3])
        R = w.exp().matrix()
        t = self.w[3:6]

        theta = w.magnitude()
        I = np.eye(3)
        wx = w.matrix()
        wx2 = wx.dot(wx)
        A = (1 - np.cos(theta)) / (theta ** 2) if theta != 0 else 1 / 2.0
        B = (theta - np.sin(theta)) / (theta ** 3) if theta != 0 else 1 / 6.0
        V = I + A * wx + B * wx2
        t = V.dot(t)

        T = np.eye(4)
        T[0:3, 0:3] = R
        T[0:3, 3] = t
        return SE3(T)

    def oldV(self):
        # TODO outdated?
        w = so3(self.w[3:5, 0])
        wx = w.matrix()
        theta = w.magnitude()
        cs = np.cos(theta)
        sn = np.sin(theta)
        I3 = np.eye(3, 3)
        V = I3 + (1 - cs) / theta ** 2 * wx + (theta - sn) / theta ** 3 * wx * wx
        return V

    def oldExp(self):
        # TODO outdated?
        t = self.w[0:2, 0]
        w = so3(self.w[3:5, 0])
        V = self.V()
        Vt = V * t
        T = np.eye(4, 4)
        T[0:2, 0:2] = w.exp()
        T[0:2, 3] = Vt
        return SE3(T)

    def oldX(self):
        # TODO outdated?
        t = self.w[0:2, 0]
        w = so3(self.w[3:5, 0])
        wx = w.matrix()
        A = np.zeros(4, 4)
        A[0:2, 0:2] = wx
        A[0:2, 3] = t
        return A
