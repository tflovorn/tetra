from scipy.optimize import bisect
from tetra.numstates import NumStates
from tetra.ksample import OptimizeGs, MakeEks
from tetra.submesh import MakeSubmesh, MakeTetra

def FindFermiToTol(n0, Efn, R, num_electrons, tol=None, tetras0=None, Eks0=None):
    '''Returns the Fermi energy E_F, at which the integrated number of
    states n(E_F) = num_electrons.
    Assumes that Efn contains all the electronic states of the system;
    e.g. for a collinear spin-polarized system, Eks must not contain only
    up-spins or only down-spins but instead contain all states of both types.

    n0 = initial value for Brillouin zone submesh density (the total number of
    k-points sampled is (n+1)**3).

    Efn = a function E(k) which returns a list of the band energies at the
    Brillouin zone point k, with the returned list sorted in ascending order;
    k is expressed in the reciprocal lattice basis.

    R = a numpy matrix with rows given by the reciprocal lattice vectors.

    tol = if given, accuracy to determine the Fermi energy to; quit if
            |E_F(n) - E_F(2*n)| < tol.

    tetras0 = pre-determined tetrahedron list at n=n0.

    Eks0 = pre-determined band energy sampled over submesh with n=n0.
    '''
    # Get optimal reciprocal lattice orientation.
    G_order, G_neg = OptimizeGs(R)

    val, old_val = None, None
    if tetras0 != None and Eks0 != None:
        val = FindFermi(num_electrons, tetras0, Eks0)
    else:
        # Generate submesh and tetrahedra.
        submesh = MakeSubmesh(n0)
        tetras = MakeTetra(n0)
        # Sample E(k).
        Eks = MakeEks(Efn, submesh, G_order, G_neg)
        # Get E_F.
        val = FindFermi(num_electrons, tetras, Eks)
    print("Got E_Fermi = {} at n = {}".format(val, n0))

    if tol == None:
        return val

    n = n0
    while old_val == None or abs(val - old_val) > tol:
        old_val = val
        n *= 2
        # Generate submesh and tetrahedra.
        submesh = MakeSubmesh(n)
        tetras = MakeTetra(n)
        # Sample E(k).
        Eks = MakeEks(Efn, submesh, G_order, G_neg)
        # Get E_F.
        val = FindFermi(num_electrons, tetras, Eks)
        print("Got E_Fermi = {} at n = {}".format(val, n))

    return val

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
        print("At E = ", E, " count = ", count)
        return count - num_electrons
    emin = _minimum_E(Eks)
    emax = _maximum_E(Eks)
    E_Fermi = bisect(statecount_error, emin, emax)
    return E_Fermi

def _minimum_E(Eks):
    minval = None
    for Es_at_k in Eks:
        for energy in Es_at_k:
            if minval == None or energy < minval:
                minval = energy
    return minval

def _maximum_E(Eks):
    maxval = None
    for Es_at_k in Eks:
        for energy in Es_at_k:
            if maxval == None or energy > maxval:
                maxval = energy
    return maxval
