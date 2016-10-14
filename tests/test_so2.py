from unittest import TestCase
from Lie import *
import numpy as np
from tests.auxiliaryOperations import equalWithError
from auxiliaryOperations import test_expLog, testvalues


class TestSo2(TestCase):
    def setUp(self):
        self.testElements = []
        for value in testvalues:
            self.testElements.append(so2(theta=value))

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
