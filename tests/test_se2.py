from unittest import TestCase
from Lie import *
from auxiliaryOperations import *


class TestSe2(TestCase):
    def setUp(self):
        self.testElements = []
        for theta in testvalues:
            for x in smallvalues:
                for y in smallvalues:
                    self.testElements.append(se2(vector=np.array([theta, x, y])))

    def test_expLog(self):
        ok = test_expLog(self.testElements)
        self.assertTrue(ok, "Error with expLog")

    def test_add(self):
        self.fail("Still not decided which order to use in adding, as it means composing functions")
