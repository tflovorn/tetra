import numpy as np
from tetra.ksample import OptimizeGs, MakeEks
from tetra.submesh import MakeSubmesh, MakeTetra
from tetra.numstates import _tetra_Es

def DosValues(Emin, Emax, num_Es, n, Efn, R):
    '''Return a list of D(E) values giving the density of states at energy E
    summed over all tetrahedra and band indices; E ranges over num_Es equally
    spaced values from Emin to Emax.

    Also returns the list of energy values used and lists of tetrahedra and
    E(k) values so that these can be used in further calculations (such as
    FindFermi()).

    The calculation of D(E) is implemented as described in BJA94 Appendix C.

    n = Brillouin zone submesh density (the total number of k-points sampled
    is (n+1)**3).

    Efn = a function E(k) which returns a list of the band energies at the
    Brillouin zone point k, with the returned list sorted in ascending order;
    k is expressed in the reciprocal lattice basis.

    R = a numpy matrix with rows given by the reciprocal lattice vectors.
    '''
    tetras, Eks = _dos_setup(n, Efn, R)
    # Get D(E) values.
    E_vals = np.linspace(Emin, Emax, num_Es)
    dos_vals = []
    for E in E_vals:
        dos_vals.append(Dos(E, tetras, Eks))
    return dos_vals, E_vals, tetras, Eks

def _dos_setup(n, Efn, R):
    # Get optimal reciprocal lattice orientation.
    G_order, G_neg = OptimizeGs(R)
    # Generate submesh and tetrahedra.
    submesh = MakeSubmesh(n)
    tetras = MakeTetra(n)
    # Sample E(k).
    Eks = MakeEks(Efn, submesh, G_order, G_neg)
    return tetras, Eks

def DosValuesPerBand(Emin, Emax, num_Es, n, Efn, R):
    '''Return a list of D_i(E) values giving the density of states at energy E
    summed over all tetrahedra separated by band index i; E ranges over num_Es
    equally spaced values from Emin to Emax.

    Also returns the list of energy values used and lists of tetrahedra and
    E(k) values so that these can be used in further calculations (such as
    FindFermi()).

    The calculation of D(E) is implemented as described in BJA94 Appendix C.

    n = Brillouin zone submesh density (the total number of k-points sampled
    is (n+1)**3).

    Efn = a function E(k) which returns a list of the band energies at the
    Brillouin zone point k, with the returned list sorted in ascending order;
    k is expressed in the reciprocal lattice basis.

    R = a numpy matrix with rows given by the reciprocal lattice vectors.
    '''
    tetras, Eks = _dos_setup(n, Efn, R)
    # Get D(E) values.
    E_vals = np.linspace(Emin, Emax, num_Es)
    dos_vals = []
    for E in E_vals:
        dos_vals.append(DosPerBand(E, tetras, Eks))
    return dos_vals, E_vals, tetras, Eks

def DosPerBand(E, tetras, Eks):
    '''Return a list with elements D_i(E), the density of states at energy E
    summed over all tetrahedra separated by band index i.
    The calculation of D_T(E) is implemented as described in BJA94 Appendix C.

    tetras = a list of tuples of the form (kN1, kN2, kN3, kN4) denoting the
    vertices of tetrahedra to include in the summation, where the kN's are
    indices of submesh (i.e. submesh[kN1] = k1, etc.). The tetrahedra must be
    constructed as described in BJA94 Section III.

    Eks = a list in which each element is a sorted list of eigenstate energies
    E_n(k), with k being the k-point at the corresponding element of submesh
    (i.e. Eks[kN][band_index] = E_n(k)).
    '''
    num_bands = len(Eks[0])
    num_tetra = len(tetras)
    Dis = [0.0]*num_bands
    for tet in tetras:
        for band_index in range(num_bands):
            Dis[band_index] += DosContrib(E, tet, num_tetra, Eks, band_index)
    return Dis

def Dos(E, tetras, Eks):
    '''Return D(E), the density of states at energy E summed over all
    tetrahedra and band indices.
    The calculation of D_T(E) is implemented as described in BJA94 Appendix C.

    tetras = a list of tuples of the form (kN1, kN2, kN3, kN4) denoting the
    vertices of tetrahedra to include in the summation, where the kN's are
    indices of submesh (i.e. submesh[kN1] = k1, etc.). The tetrahedra must be
    constructed as described in BJA94 Section III.

    Eks = a list in which each element is a sorted list of eigenstate energies
    E_n(k), with k being the k-point at the corresponding element of submesh
    (i.e. Eks[kN][band_index] = E_n(k)).
    '''
    num_bands = len(Eks[0])
    num_tetra = len(tetras)
    D = 0.0
    for tet in tetras:
        for band_index in range(num_bands):
            D += DosContrib(E, tet, num_tetra, Eks, band_index)
    return D

def DosContrib(E, tetra, num_tetra, Eks, band_index):
    '''Return the contribution to the density of states at energy E from the
    specified band index and tetrahedron.
    The calculation of D_T(E) is implemented as described in BJA94 Appendix C.

    tetra = a tuple of the form (kN1, kN2, kN3, kN4) denoting the
    vertices of a tetrahedron to include in the Brillouin zone summation,
    where the kN's are indices of submesh (i.e. submesh[kN1] = k1, etc.).
    The tetrahedra must be constructed as described in BJA94 Section III.

    num_tetra = total number of tetrahedra in the full Brillouin zone.
    Equal to (volume of tetrahedron) / (volume of full BZ).

    Eks = a list in which each element is a sorted list of eigenstate energies
    E_n(k), with k being the k-point at the corresponding element of submesh
    (i.e. Eks[kN][band_index] = E_n(k)).

    band_index = the band index n to consider, corresponding to n in E_n(k).
    '''
    E1, E2, E3, E4 = _tetra_Es(tetra, band_index, Eks)
    if E <= E1:
        return 0.0
    elif E1 <= E <= E2:
        if E1 == E2:
            return 0.0
        return (1/num_tetra) * 3*(E - E1)**2/((E2 - E1)*(E3 - E1)*(E4 - E1))
    elif E2 <= E <= E3:
        if E2 == E3:
            return (1.0 / num_tetra) * 3.0 * (E2 - E1) / ((E3 - E1) * (E4 - E1))
        fac = (1/num_tetra) / ((E3 - E1)*(E4 - E1))
        elin = 3*(E2 - E1) + 6*(E - E2)
        esq = -3*(((E3 - E1) + (E4 - E2))/((E3 - E2)*(E4 - E2))) * (E - E2)**2
        return fac * (elin + esq)
    elif E3 <= E <= E4:
        if E3 == E4:
            return 0.0
        return (1/num_tetra) * 3*(E4 - E)**2/((E4 - E1)*(E4 - E2)*(E4 - E3))
    else:
        # E >= E4
        return 0.0
