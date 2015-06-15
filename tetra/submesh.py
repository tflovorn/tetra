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
                submesh.append(k1, k2, k3)
    return submesh
