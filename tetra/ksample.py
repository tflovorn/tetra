import numpy as np

def OptimizeGs(R):
    '''Rearrange the reciprocal lattice vectors such that the Cartesian
    distance from submesh cell point 3 to point 6 (as depicted in BJA94
    Fig. 5) is minimized, in order to minimize interpolation error.
    The reciprocal lattice vectors may be permuted and/or negated.

    R = a numpy matrix with rows given by the reciprocal lattice vectors.

    Returns two vectors G_order = (o0, o1, o2) and G_neg = (n0, n1, n2).
    The value of G_order is drawn from the permutations of (0, 1, 2) and
    indicates the optimal permutation of the reciprocal lattice vectors,
    when combined with G_neg which indicates whether the corresponding
    reciprocal lattice vector is negated.
    Concretely, R_opt[o0, :] = G_neg[0]*R[0, :] and similarly for R[1, :]
    and R[2, :].
    '''
    permutations = ((0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0))
    signs = ((1, 1, 1), (1, 1, -1), (1, -1, 1), (1, -1, -1),
             (-1, 1, 1), (-1, 1, -1), (-1, -1, 1), (-1, -1, -1))
    G_order = None
    G_neg = None
    opt_36 = None
    for perm in permutations:
        for sign in signs:
            this_R = np.zeros((3, 3), dtype=np.float64)
            for i in range(3):
                this_R[perm[i], :] = sign[i] * R[i, :]
            k3_to_k6 = np.linalg.norm(-this_R[0, :] + this_R[1, :] - this_R[2, :])
            if opt_36 == None or k3_to_k6 < opt_36:
                opt_36 = k3_to_k6
                G_order = perm
                G_neg = sign
    return G_order, G_neg

def Get_k_Orig(k_opt, G_order, G_neg):
    '''Convert k_opt (a k-point in the represented in the optimal permutation
    of the reciprocal lattice as determined by OptimizeGs) into the
    corresponding k-point in the original reciprocal lattice vectors.

    The value of G_order is drawn from the permutations of (0, 1, 2) and
    indicates the optimal permutation of the reciprocal lattice vectors,
    when combined with G_neg which indicates whether the corresponding
    reciprocal lattice vector is negated.
    Concretely, R_opt[o0, :] = G_neg[0]*R[0, :] and similarly for R[1, :]
    and R[2, :].
    '''
    k = [0.0]*3
    for i in range(3):
        k[i] = k_opt[G_order[i]] * G_neg[i]
    return k

def GetRopt(R, G_order, G_neg):
    '''Return the value of Ropt corresponding to G_order and G_neg.

    The value of G_order is drawn from the permutations of (0, 1, 2) and
    indicates the optimal permutation of the reciprocal lattice vectors,
    when combined with G_neg which indicates whether the corresponding
    reciprocal lattice vector is negated.
    Concretely, R_opt[o0, :] = G_neg[0]*R[0, :] and similarly for R[1, :]
    and R[2, :].
    '''
    R_opt = np.zeros((3, 3), dtype=np.float64)
    for i in range(3):
        R_opt[G_order[i], :] = G_neg[i]*R[i, :]
    return R_opt

def MakeEks(Efn, submesh, G_order=None, G_neg=None):
    '''Generate Eks, a list in which each element is a sorted list of
    eigenstate energies E_n(k), with k being the k-point at the corresponding
    element of submesh (i.e. Eks[kN][band_index] = E_n(k)).

    Efn = a function E(k) which returns a list of the band energies at the
    Brillouin zone point k, with the returned list sorted in ascending order;
    k is expressed in the reciprocal lattice basis.

    submesh = a list of k-points covering the Brillouin zone, determined as
    described in BJA94 Section III.

    The value of G_order is drawn from the permutations of (0, 1, 2) and
    indicates the optimal permutation of the reciprocal lattice vectors,
    when combined with G_neg which indicates whether the corresponding
    reciprocal lattice vector is negated.
    Concretely, R_opt[o0, :] = G_neg[0]*R[0, :] and similarly for R[1, :]
    and R[2, :].
    '''
    Eks = []
    for k in submesh:
        k_orig = k
        if G_order != None and G_neg != None:
            k_orig = Get_k_orig(k, G_order, G_neg)
        Es = Efn(k_orig)
        _check_sort(Es)
        Eks.append(Es)
    return Eks

def _check_sort(Es):
    for i, E in enumerate(Es):
        if i > 0 and Es[i] < Es[i-1]:
            raise ValueError("Energies returned by Efn not sorted.")

def MakeXks(Xfn, submesh, G_order=None, G_neg=None):
    '''Generate Xks, a list in which each element is a list of matrix elements
    X_n(k), with k being the k-point at the corresponding element of submesh
    (i.e. Eks[kN][band_index] = E_n(k)).

    Xfn = a function X(k) which returns a list of the values of the matrix
    elements of the operator X at k; the returned values are ordered in the
    same way as the band energies and k is expressed in the reciprocal
    lattice basis.

    submesh = a list of k-points covering the Brillouin zone, determined as
    described in BJA94 Section III.

    The value of G_order is drawn from the permutations of (0, 1, 2) and
    indicates the optimal permutation of the reciprocal lattice vectors,
    when combined with G_neg which indicates whether the corresponding
    reciprocal lattice vector is negated.
    Concretely, R_opt[o0, :] = G_neg[0]*R[0, :] and similarly for R[1, :]
    and R[2, :].
    '''
    Xks = []
    for k in submesh:
        k_orig = k
        if G_order != None and G_neg != None:
            k_orig = Get_k_orig(k, G_order, G_neg)
        Xs = Xfn(k_orig)
        Xks.append(Xs)
    return Xks
