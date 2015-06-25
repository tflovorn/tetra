from tetra.numstates import _tetra_Es

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
