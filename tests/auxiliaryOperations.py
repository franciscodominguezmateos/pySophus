import numpy as np

testvalues = np.linspace(-2, 2, num=1 + 4 * 8) * np.pi
tens = np.power(10, np.linspace(0, 10, 11))
smallvalues = list(map(lambda a: 1 / a, tens))
testvalues = np.append(testvalues, smallvalues)


def twoOperations(elements, op1, op2):
    passed = True
    for e in elements:
        a = e.matrix()
        b = eval("e." + op1 + "()." + op2 + "().matrix()")
        ok = equalWithError(a, b)
        if not ok:
            print("Error with " + str(e.matrix()) + " in " + op1 + "," + op2)
            passed = False
    return passed


def test_expLog(elements):
    return twoOperations(elements, "exp", "log")


def test_logExp(elements):
    return twoOperations(elements, "log", "exp")


def equalWithError(a, b, error=0.00000001):
    # used np.linalg.norm to be able to compare even matrices, not only numbers
    return np.linalg.norm(a - b) < error
