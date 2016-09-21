from unittest import TestCase
from Lie import *
import numpy as np
from test import test_logExp


class TestSO2(TestCase):
    def setUp(self):

        thetas = np.linspace(-2, 2, num=1 + 4 * 8) * np.pi
        self.testElements = []
        for theta in thetas:
            cs = np.cos(theta)
            sn = np.sin(theta)
            self.testElements.append(SO2(np.array([[cs, -sn],
                                                   [sn, cs]])))

        tens = np.power(10, np.linspace(0, 10, 11))
        for num in tens:
            smallNum = 1 / num
            cs = np.cos(smallNum)
            sn = np.sin(smallNum)
            self.testElements.append(SO2(np.array([[cs, -sn],
                                                   [sn, cs]])))

    def test_logExp(self):
        ok = test_logExp(self.testElements)
        self.assertTrue(ok, "Error with logExp")

    def test_multiplication(self):
        print("Multiplication test included in adding test in so2 algebra")
        self.assertTrue(True)

    def test_J(self):
        self.fail("Not implemented")
