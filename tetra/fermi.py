from scipy.optimize import bisect
from tetra.numstates import NumStates

def FindFermi(num_electrons, tetras, Eks):
    '''Returns the Fermi energy E_F, at which the integrated number of
    states n(E_F) = num_electrons.
    Assumes that Eks contains all the electronic states of the system;
    e.g. for a collinear spin-polarized system, Eks must not contain only
    up-spins or only down-spins but instead contain all states of both types.

    tetras = a list of tuples of the form (kN1, kN2, kN3, kN4) denoting the
    vertices of tetrahedra to include in the summation, where the kN's are
    indices of submesh (i.e. submesh[kN1] = k1, etc.). The tetrahedra must be
    constructed as described in BJA94 Section III.

    Eks = a list in which each element is a sorted list of eigenstate energies
    E_n(k), with k being the k-point at the corresponding element of submesh
    (i.e. Eks[kN][band_index] = E_n(k)).
    '''
    def statecount_error(E):
        count = NumStates(E, tetras, Eks)
        return count - num_electrons
    emin = _minimum_E(Eks)
    emax = _maximum_E(Eks)
    E_Fermi = bisect(statecount_error, emin, emax)
    return E_Fermi

def _minimum_E(Eks):
    minval = None
    for Es_at_k in Eks:
        for energy in E_at_k:
            if minval == None or energy < minval:
                minval = energy
    return minval

def _maximum_E(Eks):
    maxval = None
    for Es_at_k in Eks:
        for energy in E_at_k:
            if maxval == None or energy > maxval:
                maxval = energy
    return maxval
