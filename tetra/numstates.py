def NumStatesBands(E, submesh, tetras, Eks):
    '''Return a list of the contributions from each band index to the number
    of states with energy less than or equal to E, summed over all tetrahedra.
    The calculation of n(E) is implemented as described in BJA94 Appendix A.

    submesh = a list of k-points covering the Brillouin zone, determined as
    described in BJA94 Section III.

    tetras = a list of tuples of the form (kN1, kN2, kN3, kN4) denoting the
    vertices of tetrahedra to include in the summation, where the kN's are
    indices of submesh (i.e. submesh[kN1] = k1, etc.). The tetrahedra must be
    constructed as described in BJA94 Section III.

    Eks = a list in which each element is a sorted list of eigenstate energies
    E_n(k), with k being the k-point at the corresponding element of submesh
    (i.e. Eks[kN][band_index] = E_n(k)).
    '''
    pass

def NumStatesContrib(E, tetra, num_tetra, Eks, band_index):
    '''Return the contribution to the number of states with energy less than
    or equal to E (i.e. the integrated density of states n(E)) from the
    specified tetrahedron and band index.
    The calculation of n(E) is implemented as described in BJA94 Appendix A.

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
    elif E1 < E < E2:
        return num_tetra * (E - E1)**3 / ((E2 - E1)*(E3 - E1)*(E4 - E1))
    elif E2 < E < E3:
        fac = num_tetra / ((E3 - E1)*(E4 - E1))
        esq = (E2 - E1)**2 + 3*(E2 - E1)*(E - E2) + 3*(E - E2)**2
        ecub = -(((E3 - E1) + (E4 - E2))/((E3 - E2)*(E4 - E2))) * (E - E2)**3
        return fac * (esq + ecub)
    elif E3 < E < E4:
        return num_tetra * (1 - (E4 - E)**3/((E4 - E1)*(E4 - E2)*(E4 - E3)))
    else:
        # E >= E4
        return num_tetra

def _tetra_Es(tetra, band_index, Eks):
    '''Return a sorted list of the eigenstate energies at the vertices of the
    specified tetrahedron.
    '''
    Es = []
    for kN in tetra:
        Es.append(Eks[kN][band_index])
    return sorted(Es)
