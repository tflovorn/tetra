def Weights(E_Fermi, submesh, tetras, Eks):
    '''Return a list in which each element is a list of integration weights.
    The first index for the returned list specifies a band index, and the
    second list specifies a k-point index; i.e. the returned list
    w[n][j] = w_{nj}.
    The calculation of w_{nj} is implemented as described in BJA94 Appendix B
    and Section V.

    E_Fermi = Fermi energy of the system.

    submesh = a list of k-points covering the Brillouin zone, determined as
    described in BJA94 Section III.

    tetras = a list of tuples of the form (kN1, kN2, kN3, kN4) denoting the
    vertices of tetrahedra to include in the summation, where the kN's are
    indices of submesh (i.e. submesh[kN1] = k1, etc.). The tetrahedra must be
    constructed as described in BJA94 Section III.

    Eks = a list in which each element is a sorted list of eigenstate energies
    E_n(k), with k being the k-point at the corresponding element of submesh.

    bandIndex = the band index n to consider, corresponding to n in E_n(k).
    '''
    pass

def WeightContrib(E_Fermi, submesh, tetra, Eks, kN, bandIndex):
    '''Return the specified tetrahedron's contribution to the integration
    weight w_{bandIndex, kN}.
    The calculation of w_{nj} is implemented as described in BJA94 Appendix B
    and Section V.

    E_Fermi = Fermi energy of the system.

    submesh = a list of k-points covering the Brillouin zone, determined as
    described in BJA94 Section III.

    tetra = a tuple of the form (kN1, kN2, kN3, kN4) denoting the
    vertices of a tetrahedron to include in the Brillouin zone summation,
    where the kN's are indices of submesh (i.e. submesh[kN1] = k1, etc.).
    The tetrahedra must be constructed as described in BJA94 Section III.

    Eks = a list in which each element is a sorted list of eigenstate energies
    E_n(k), with k being the k-point at the corresponding element of submesh.

    bandIndex = the band index n to consider, corresponding to n in E_n(k).
    '''
    if kN not in tetra:
        return 0.0
    pass
