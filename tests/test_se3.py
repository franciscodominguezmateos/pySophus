from unittest import TestCase
from auxiliaryOperations import equalWithError
from Lie import *


class TestSe3(TestCase):
    def setUp(self):
        testvalues = np.linspace(-0.1, 0.1, num=1 + 16) * np.pi
        tens = np.power(10, np.linspace(10, 21, 5))
        smallvalues = list(map(lambda a: 1 / a, tens))
        testvalues = np.append(testvalues, smallvalues)

        self.testElements = []
        for wx in testvalues:
            for wy in testvalues:
                for wz in testvalues:
                    rotationVector=np.array([wx, wy, wz])
                    for x in smallvalues:
                        for y in smallvalues:
                            for z in smallvalues:
                                translationVector = np.array([x,y,z])
                                self.testElements.append(se3(vector=np.append(rotationVector,translationVector)))

    def test_rotations(self):
        passed = True
        for e in self.testElements:
            a = e.exp().matrix()
            b = e.exp().log().exp().matrix()
            # error value passed because with the default it would fail even if the matrices are almost the same
            ok = equalWithError(a, b, error=0.0001)
            if not ok:
                passed = False
                a = e.exp().matrix()
                b = e.exp().log().exp().matrix()
                print("Error with ", e.vector())
        self.assertTrue(passed, "Error in rotations")

    def test_add(self):
        self.fail("Not implemented")
