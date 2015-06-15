import unittest
import numpy as np
from tetra.sum import SumEnergy

def _cubicR(a):
    Da = (a, 0.0, 0.0)
    Db = (0.0, a, 0.0)
    Dc = (0.0, 0.0, a)
    D = np.array((Da, Db, Dc)).T
    R = 2*np.pi*np.linalg.inv(D)
    return R

def _simplebands_Efn(k, num_bands, t, E0, deltaE):
    if deltaE < 0.0:
        raise ValueError("Should have deltaE >= 0.0")
    En0s = []
    for n in range(num_bands):
        En0s.append(E0 + n*deltaE)
    Eks = []
    tk = -2*t*(np.cos(2.0*np.pi*k[0]) + np.cos(2.0*np.pi*k[1]) + np.cos(2.0*np.pi*k[2]))
    for n in range(num_bands):
        Eks.append(En0s[n] + tk)
    return Eks

class TestSumEnergy(unittest.TestCase):
    def test_energy_simpleinsulator(self):
        n = 8
        a = 1.0
        R = _cubicR(a)
        num_bands = 2
        num_electrons = 1
        t = 1.0
        E0 = 6.0
        deltaE = 14.0   # need deltaE > bandwidth = 12t for insulator
        def Efn(k):
            return _simplebands_Efn(k, num_bands, t, E0, deltaE)
        tolerance = 1e-6
        result, n, ws = SumEnergy(n, Efn, R, num_electrons, tolerance)

        expected = E0
        _assert_within(self, result, expected, tolerance*1.1)

def _assert_within(testcase, result, expected, eps):
    err = abs(result - expected)
    testcase.assertTrue(err < eps)

if __name__ == "__main__":
    unittest.main()
