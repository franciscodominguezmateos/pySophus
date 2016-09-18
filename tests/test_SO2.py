from unittest import TestCase
from Lie import *
import numpy as np


class TestSO2(TestCase):
    def setUp(self):
        self.I = SO2(np.eye(2))
        self.half = SO2(np.array([[-1, 0],
                                  [0, -1]]))
        cs = np.cos(np.pi / 2)
        sn = np.sin(np.pi / 2)
        self.quarter = SO2(np.array([[cs, -sn], [sn, cs]]))
        cs = np.cos(np.pi / 4)
        sn = np.sin(np.pi / 4)
        self.eighth = SO2(np.array([[cs, -sn], [sn, cs]]))

    def tearDown(self):
        pass

    def test_matrix(self):
        ok = np.array_equal(self.I.matrix(), np.eye(2))
        self.assertTrue(ok, "Error with Identity matrix")

        ok = np.array_equal(self.half.matrix(), np.array([[-1, 0],
                                                          [0, -1]]))
        self.assertTrue(ok, "Error with PI rotation")

        cs = np.cos(np.pi / 2)
        sn = np.sin(np.pi / 2)
        ok = np.array_equal(self.quarter.matrix(), np.array([[cs, -sn], [sn, cs]]))
        self.assertTrue(ok, "Error with PI/2 rotation")

        cs = np.cos(np.pi / 4)
        sn = np.sin(np.pi / 4)
        ok = np.array_equal(self.eighth.matrix(), np.array([[cs, -sn], [sn, cs]]))
        self.assertTrue(ok, "Error with PI/4 rotation")

    def test_log(self):
        self.fail("Not implemented")

    def test_J(self):
        self.fail("Not implemented")
