'''
Created on Mar 2, 2016

@author: Francisco Dominguez
Lie group SO(2) a 2x2matrix
Lie algebra so(2) a number
'''
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
    def __init__(self, matrix):
        self.M = matrix

    def __mul__(self, other):
        return SO2(self.M.dot(other.matrix()))

    def matrix(self):
        return self.M

    def log(self):
        return so2(theta=np.arctan2(self.M[1, 0], self.M[0, 0]))

    def J(self):
        # TODO realmente tendria que estar aqui?
        return so2_G * self.M


so2_G = np.matrix([[0, -1],
                   [1, 0]])


class so2(Algebra):
    def __init__(self, **kwargs):
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

    def __add__(self, op):
        R1 = self.exp()
        R2 = op.exp()
        R = R1 * R2
        return R.log()

    def exp(self):
        theta = self.angle
        cs = np.cos(theta)
        sn = np.sin(theta)
        R = np.matrix([[cs, -sn],
                       [sn, cs]])
        return SO2(matrix=R)

    def matrix(self):
        return self.angle * so2_G

    def vector(self):
        return np.array([self.angle])

    def theta(self):
        return self.angle


class SE2(Group):
    def __init__(self, matrix):
        self.M = matrix

    def matrix(self):
        return self.M

    def log(self):
        # TODO testear caso theta = 0
        w = SO2(self.M[0:2, 0:2])
        t = self.M[0:2, 2]
        theta = w.log().angle()
        cs = np.cos(theta)
        sn = np.sin(theta)
        A = sn / theta if theta != 0 else 1
        B = (1 - cs) / theta if theta != 0 else 0
        V1 = 1 / (A ** 2 + B ** 2) * np.matrix([[A, B],
                                                [-B, A]])
        Vt = V1 * t
        return se2(np.array([theta, Vt[0], Vt[1]]))

    def __mul__(self, other):
        return SE2(self.M.dot(other.matrix()))


class se2(Algebra):
    # TODO extract algebra base from class to avoid overusing memory
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
        if "matrix" in kwargs:
            self.w = np.array([0, 0, 0])
            self.w[0] = kwargs["matrix"][1, 0]
            self.w[1:3] = kwargs["matrix"][0:2, 2]
        elif "vector" in kwargs:
            self.w = kwargs["vector"]
        else:
            raise TypeError("Argument must be matrix or vector")

    def __add__(self, op):
        R1 = self.exp()
        R2 = op.exp()
        R = R1 * R2
        return R.log()

    def exp(self):
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
        return self.w[0] * se2.G1 + self.w[1] * se2.G2 + self.w[2] * se2.G3

    def vector(self):
        return self.w


class SO3(Group):
    def __init__(self, matrix):
        self.M = matrix

    def log(self):
        if np.array_equal(np.eye(3), self.M):
            # case theta == 0
            return so3(np.array([0, 0, 0]))
        elif (np.trace(self.M) == -1):
            # case theta == pi
            w = 1 / np.sqrt(2 * (1 + self.M[2, 2])) * np.array([self.M[0, 2], self.M[1, 2], 1 + self.M[2, 2]])
            return so3(w)
        else:
            cs = (self.M.trace() - 1) / 2
            theta = np.acos(cs)
            sn = np.sin(theta)
            logR = theta / (2 * sn) * (self.M - self.M.T)
            return so3(vector=np.array([logR[2, 1], logR[0, 2], logR[1, 0]]))

    def matrix(self):
        return self.M

    def __mul__(self, other):
        return SO3(self.M.dot(other.matrix()))


class so3(Algebra):
    # TODO extract algebra base from class to avoid overusing memory
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

    def __add__(self, op):
        R1 = self.exp()
        R2 = op.exp()
        R = R1 * R2
        return R.log()

    def magnitude(self):
        return np.linalg.norm(self.w)

    def exp(self):
        wx = self.matrix()
        theta = np.linalg.norm(self.w)
        cs = np.cos(theta)
        sn = np.sin(theta)
        I = np.eye(3)
        R = I + (sn / theta) * wx + ((1 - cs) / theta ** 2) * wx * wx if theta != 0 else I
        return SO3(R)

    def vector(self):
        return self.w

    def matrix(self):
        return so3.G1 * self.w[0] + so3.G2 * self.w[1] + so3.G3 * self.w[2]


class SE3(Group):
    def __init__(self, matrix):
        self.M = matrix

    def __mul__(self, other):
        return SE3(self.M.dot(other.matrix()))

    def matrix(self):
        return self.M

    def log(self):
        R = self.M[0:3, 0:3]
        w = SO3(R).log()
        t = self.M[0:3, 3]

        I = np.eye(3)
        theta = w.magnitude()
        wx = w.matrix()
        wx2 = wx.dot(wx)
        A = np.sin(theta) / theta if theta != 0 else 1
        B = (1 - np.cos(theta)) / (theta ** 2) if theta != 0 else 1 / 2
        V = I - 1 / 2 * wx + 1 / (theta ** 2) * (1 - A / (2 * B)) * wx2 if theta != 0 else I

        v = V.dot(t)
        return se3(vector=np.append(w, v))


class se3(Algebra):
    # TODO extract algebra base from class to avoid overusing memory
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

    def __add__(self, op):
        R1 = self.exp()
        R2 = op.exp()
        R = R1 * R2
        return R.log()

    def vector(self):
        return self.w

    def matrix(self):
        x = np.zeros((4, 4))
        for i in range(6):
            x += self.w[i] * eval("G" + str(i + 1))
        return x

    def exp(self):
        w = so3(vector=self.w[0:3])
        R = w.exp().matrix()
        t = self.w[3:6]

        theta = w.magnitude()
        I = np.eye(3)
        wx = w.matrix()
        wx2 = wx.dot(wx)
        A = (1 - np.cos(theta)) / (theta ** 2) if theta != 0 else 1 / 2
        B = (theta - np.sin(theta)) / (theta ** 3) if theta != 0 else 1 / 6
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
        cs = math.cos(theta)
        sn = math.sin(theta)
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
