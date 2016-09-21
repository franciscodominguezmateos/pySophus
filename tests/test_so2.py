from unittest import TestCase
from Lie import *
import numpy as np
from tests.auxiliaryOperations import equalWithError
from auxiliaryOperations import test_expLog


class TestSo2(TestCase):
    def setUp(self):
        thetas = np.linspace(-2, 2, num=1 + 4 * 8) * np.pi
        self.testElements = []
        for theta in thetas:
            self.testElements.append(so2(theta=theta))

        tens = np.power(10, np.linspace(0, 10, 11))
        for num in tens:
            smallNum = 1 / num
            self.testElements.append(so2(theta=smallNum))

    def test_expLog(self):
        ok = test_expLog(self.testElements)
        self.assertTrue(ok, "Error with expLog")

    def test_add(self):
        testing = []
        for a in self.testElements:
            for b in self.testElements:
                testing.append({"first": a, "second": b, "expected": so2(theta=(a.theta() + b.theta()))})
        passed = True
        for element in testing:
            testAdd = element["first"] + element["second"]
            ok = equalWithError(testAdd.exp().matrix(), element["expected"].exp().matrix())
            # We test the matrix of the group instead of the theta values of the algebra to avoid some weird behavior
            # with the limit values pi/-pi
            if not ok:
                print("Error adding " + str(element["first"].theta()) + " and " + str(
                    element["second"].theta()) + " which add " + str(
                    element["first"].theta() + element["second"].theta()))
                passed = False
        self.assertTrue(passed, "Error with + operator")
