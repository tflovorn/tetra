import numpy as np

def MakeSubmesh(n):
    '''Return a list containing the submesh of k-points in the reciprocal
    lattice basis covering the full Brillouin zone. The number of k-points
    in all 3 dimensions is given by n+1 (i.e. we require that n1 = n2 = n3);
    the total number of k-points in the submesh is (n+1)**3.

    The submesh generation is implemented as described in BJA94 Section III.
    As noted there, the submesh should have the symmetry of the reciprocal
    lattice; this restricts the allowed values of n1, n2, and n3. This
    restriction can always be satisfied by setting n1 = n2 = n3.
    We assume the submesh shift (i0, j0, k0) is zero.
    '''
    if n <= 0:
        raise ValueError("Must have n > 0 in MakeSubmesh.")
    submesh = []
    step = 1/n
    for k in range(n+1):
        k3 = k*step
        for j in range(n+1):
            k2 = j*step
            for i in range(n+1):
                k1 = i*step
                submesh.append((k1, k2, k3))
    return submesh

def MakeTetra(n):
    '''Return a list containing the tetrahedra dividing the full Brillouin
    zone. The number of k-points in all 3 dimensions is given by n+1 (i.e. we
    require that n1 = n2 = n3); the total number of k-points in the submesh
    is (n+1)**3.

    The submesh is defined as that returned by MakeSubmesh(n).

    Tetrahedra generation is implemented as described in BJA94 Section III.
    '''
    tetra = []
    subcell_tetras = [(1, 2, 3, 6), (1, 3, 5, 6), (3, 5, 6, 7), (3, 6, 7, 8),
                      (3, 4, 6, 8), (2, 3, 4, 6)]
    # Range over n instead of n+1 here to ensure that the cell is always in
    # the first Brillouin zone.
    # The number of submesh cells is given by n**3.
    for k in range(n):
        for j in range(n):
            for i in range(n):
                points = [(i, j, k), (i+1, j, k), (i, j+1, k), (i+1, j+1, k),
                          (i, j, k+1), (i+1, j, k+1), (i, j+1, k+1), (i+1, j+1, k+1)]
                for sc_t in subcell_tetras:
                    vertices = []
                    for point_index in sc_t:
                        this_i, this_j, this_k = points[point_index-1]
                        vertices.append(_submesh_index(n, this_i, this_j, this_k))
                    tetra.append(vertices)
    return tetra

def _submesh_index(n, i, j, k):
    return i + j*(n+1) + k*(n+1)**2

def _submesh_ijk(n, index):
    i = index % (n+1)
    j = ((index % ((n+1)**2)) - i) / (n+1)
    k = (index - i - j*(n+1)) / (n+1)**2
    return i, j, k
