import numpy as np


def twoOperations(elements, op1, op2):
    passed = True
    for e in elements:
        a = e.matrix()
        b = eval("e." + op1 + "()." + op2 + "().matrix()")
        ok = equalWithError(a, b)
        if not ok:
            print("Error with " + str(e.vector()) + " in " + op1 + "," + op2)
            passed = False
    return passed


def test_expLog(elements):
    return twoOperations(elements, "exp", "log")


def test_logExp(elements):
    return twoOperations(elements, "log", "exp")


def equalWithError(a, b, error=0.00000001):
    # used np.linalg.norm to be able to compare even matrices, not only numbers
    return np.linalg.norm(a - b) < error
