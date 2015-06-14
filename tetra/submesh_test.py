import unittest
import numpy as np
from tetra.submesh import OptimizeGs, Get_k_Orig

class TestOptimizeGs(unittest.TestCase):
    def test_cubic(self):
        Rcub = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
        G_order, G_neg = OptimizeGs(Rcub)
        _vec_equal(self, G_order, (0, 1, 2))
        _vec_equal(self, G_neg, (1, 1, 1))

    def test_hex(self):
        Rhex = np.array([[1.0, 0.0, 0.0], [-0.5, np.sqrt(3)/2.0, 0.0], [0.0, 0.0, 1.0]])
        G_order, G_neg = OptimizeGs(Rhex)
        _vec_equal(self, G_order, (0, 1, 2))
        _vec_equal(self, G_neg, (1, -1, 1))

class TestGet_k_Orig(unittest.TestCase):
    def test_k_no_change(self):
        G_order, G_neg = (0, 1, 2), (1, 1, 1)
        k_opt = (1.0, 2.0, 3.0)
        k_orig = Get_k_Orig(k_opt, G_order, G_neg)
        _vec_equal(self, k_opt, k_orig)

    def test_k_permute(self):
        G_order, G_neg = (1, 0, 2), (1, 1, 1)
        k_opt = (1.0, 2.0, 3.0)
        k_orig = Get_k_Orig(k_opt, G_order, G_neg)
        k_orig_expect = (2.0, 1.0, 3.0)
        _vec_equal(self, k_orig_expect, k_orig)

    def test_k_permute_and_sign(self):
        G_order, G_neg = (1, 0, 2), (1, -1, -1)
        k_opt = (1.0, 2.0, 3.0)
        k_orig = Get_k_Orig(k_opt, G_order, G_neg)
        k_orig_expect = (2.0, -1.0, -3.0)
        _vec_equal(self, k_orig_expect, k_orig)

def _vec_equal(testcase, v, u):
    testcase.assertEqual(len(v), len(u))
    for i in range(len(v)):
        testcase.assertEqual(v[i], u[i])

if __name__ == "__main__":
    unittest.main()
