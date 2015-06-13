def DosTs(E, submesh, tetras, Eks):
    '''Return a list of the contribution to the density of states at energy E
    from each tetrahedron. The Nth element of the returned list corresponds to
    the Nth tetrahedron in tetras.
    The calculation of D_T(E) is implemented as described in BJA94 Appendix C.

    submesh = a list of k-points covering the Brillouin zone, determined as
    described in BJA94 Section III.

    tetras = a list of tuples of the form (kN1, kN2, kN3, kN4) denoting the
    vertices of tetrahedra to include in the summation, where the kN's are
    indices of submesh (i.e. submesh[kN1] = k1, etc.). The tetrahedra must be
    constructed as described in BJA94 Section III.

    Eks = a list in which each element is a sorted list of eigenstate energies
    E_n(k), with k being the k-point at the corresponding element of submesh.
    '''
    pass
