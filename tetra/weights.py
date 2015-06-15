from dos import DosContrib

def Weights(E_Fermi, tetras, Eks):
    '''Return a list in which each element is a list of integration weights.
    The first index for the returned list specifies a band index, and the
    second list specifies a k-point index; i.e. the returned list
    w[n][j] = w_{nj}.
    The calculation of w_{nj} is implemented as described in BJA94 Appendix B
    and Section V.

    E_Fermi = Fermi energy of the system.

    tetras = a list of tuples of the form (kN1, kN2, kN3, kN4) denoting the
    vertices of tetrahedra to include in the summation, where the kN's are
    indices of submesh (i.e. submesh[kN1] = k1, etc.). The tetrahedra must be
    constructed as described in BJA94 Section III.

    Eks = a list in which each element is a sorted list of eigenstate energies
    E_n(k), with k being the k-point at the corresponding element of submesh
    (i.e. Eks[kN][band_index] = E_n(k)).
    '''
    num_bands = len(Eks[0])
    num_ks = len(Eks)
    num_tetra = len(tetras)
    # Initialize weights to 0.
    ws = []
    for band_index in range(num_bands):
        ws.append([])
        for k_index in range(num_ks):
            ws[band_index].append(0.0)
    # Calculate weights.
    for tet in tetras:
        for band_index in range(num_bands):
            # Get contributions to ws from this tetrahedron and band.
            tb_ws = WeightContrib(E_Fermi, tet, num_tetra, Eks, band_index)
            # Add this tetrahedron's contributions to the associated k-points.
            for i, k_index in enumerate(tet):
                ws[band_index][k_index] += tb_ws[i]

    return ws

def WeightContrib(E_Fermi, tetra, num_tetra, Eks, band_index):
    '''Return the specified tetrahedron's contribution to the integration
    weights at the k-points of the tetrahedron's vertices; i.e. return
    a list with elements w_{bandIndex, kN, tetra}, where the elements of the
    returned list range over kN values in the order specified by tetra.
    The calculation of the tetrahedron contribution to w_{nj} is implemented
    as described in BJA94 Appendix B and Section V.

    E_Fermi = Fermi energy of the system.

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
    (E1, E2, E3, E4), i_vals = _tetra_Es_ks(tetra, band_index, Eks)
    ws = [None]*4
    if E_Fermi <= E1:
        for i in range(4):
            ws[i] = 0.0
    elif E1 < E_Fermi < E2:
        C = (1 / (4*num_tetra)) * (E_Fermi - E1)**3 / ((E2 - E1)*(E3 - E1)*(E4 - E1))
        ws[0] = C * (4 - (E_Fermi - E1)*(1/(E2 - E1) + 1/(E3 - E1) + 1/(E4 - E1)))
        ws[1] = C * (E_Fermi - E1) / (E2 - E1)
        ws[2] = C * (E_Fermi - E1) / (E3 - E1)
        ws[3] = C * (E_Fermi - E1) / (E4 - E1)
    elif E2 < E_Fermi < E3:
        C1, C2, C3 = _Cs_23(E_Fermi, num_tetra, E1, E2, E3, E4)
        ws[0] = C1 + (C1 + C2)*(E3 - E_Fermi)/(E3 - E1) + (C1 + C2 + C3)*(E4 - E_Fermi)/(E4 - E1)
        ws[1] = C1 + C2 + C3 + (C2 + C3)*(E3 - E_Fermi)/(E3 - E2) + C3*(E4 - E_Fermi)/(E4 - E2)
        ws[2] = (C1 + C2)*(E_Fermi - E1)/(E3 - E1) + (C2 + C3)*(E_Fermi - E2)/(E3 - E2)
        ws[3] = (C1 + C2 + C3)*(E_Fermi - E1)/(E4 - E1) + C3*(E_Fermi - E2)/(E4 - E2)
    elif E3 < E_Fermi < E4:
        C = (1 / (4*num_tetra)) * (E4 - E_Fermi)**3 / ((E4 - E1)*(E4 - E2)*(E4 - E3))
        ws[0] = (1 / (4*num_tetra)) - C*(E4 - E_Fermi)/(E4 - E1)
        ws[1] = (1 / (4*num_tetra)) - C*(E4 - E_Fermi)/(E4 - E2)
        ws[2] = (1 / (4*num_tetra)) - C*(E4 - E_Fermi)/(E4 - E3)
        ws[3] = (1 / (4*num_tetra)) - C*(4 - (1/(E4 - E1) + 1/(E4 - E2) + 1/(E4 - E3))*(E4 - E_Fermi))
    else:
        # E_Fermi >= E4
        for i in range(4):
            ws[i] = 1 / (4*num_tetra)

    dws = _CurvatureCorrection(E_Fermi, tetra, num_tetra, Eks, band_index, (E1, E2, E3, E4))
    for i in range(4):
        ws[i] += dws[i]

    tet_ws = []
    for i in i_vals:
        tet_ws.append(ws[i])
    return tet_ws

def _Cs_23(E_Fermi, num_tetra, E1, E2, E3, E4):
    '''Return coefficients C1, C2, C3 for E2 < E_Fermi < E3.
    '''
    C1 = (1 / (4*num_tetra)) * (E_Fermi - E1)**2 / ((E4 - E1)*(E3 - E1))

    C2_num = (1 / (4*num_tetra)) * (E_Fermi - E1)*(E_Fermi - E2)*(E3 - E_Fermi)
    C2_denom = (E4 - E1)*(E3 - E2)*(E3 - E1)
    C2 = C2_num / C2_denom

    C3_num = (1 / (4*num_tetra)) * (E_Fermi - E2)**2 * (E4 - E_Fermi)
    C3_denom = (E4 - E2)*(E3 - E2)*(E4 - E1)
    C3 = C3_num / C3_denom

    return C1, C2, C3

def _CurvatureCorrection(E_Fermi, tetra, num_tetra, Eks, band_index, Es):
    '''Return a list of the curvature corrections to the k-point weight
    contributions from the specified tetrahedron. The band energies at the
    vertices of the tetrahedron are given in sorted order by Es; the returned
    corrections are given in the same order.
    '''
    D_T = DosContrib(E_Fermi, tetra, num_tetra, Eks, band_index)
    dws = [None]*4
    for i in range(4):
        dws[i] = (D_T / 40) * (sum(Es) - 4*Es[i])
    return dws

def _tetra_Es_ks(tetra, band_index, Eks):
    '''Return a sorted list of the eigenstate energies at the vertices of the
    specified tetrahedron and a list of the indices in tetra corresponding to
    the sorted energy values.
    '''
    Eis = []
    for i, kN in enumerate(tetra):
        Eis.append((Eks[kN][band_index], i))
    sorted_vals = sorted(Eis, key=lambda Ei: Ei[0])
    E_vals, i_vals = [], []
    for E, i in sorted_vals:
        E_vals.append(E)
        i_vals.append(i)
    return E_vals, i_vals
