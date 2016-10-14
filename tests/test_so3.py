from unittest import TestCase
from Lie import *
from auxiliaryOperations import equalWithError


class TestSo3(TestCase):
    def setUp(self):
        testvalues = np.linspace(-1, 1, num=1 + 4 * 8) * np.pi
        tens = np.power(10, np.linspace(0, 10, 11))
        smallvalues = list(map(lambda a: 1 / a, tens))
        testvalues = np.append(testvalues, smallvalues)

        self.testElements = []
        for wx in testvalues:
            for wy in testvalues:
                for wz in testvalues:
                    newelement = so3(vector=np.array([wx, wy, wz]))
                    self.testElements.append(newelement)

    def test_rotations(self):
        passed = True
        for e in self.testElements:
            a = e.exp().matrix()
            b = e.exp().log().exp().matrix()
            # error value passed because with the default it would fail even if the matrices are almost the same
            ok = equalWithError(a, b, error=0.0001)
            if not ok:
                passed = False
                print("Error with ", e.vector())
        self.assertTrue(passed, "Error in rotations")

    def test_add(self):
        self.fail("Not implemented")
