An implementation of the tetrahedron method for Brillouin zone summation, as
described in: [Blöchl, Jepsen, and Andersen, PRB 49, 16223 (1994)](http://journals.aps.org/prb/abstract/10.1103/PhysRevB.49.16223).
This paper is referred to in the documentation of this package as 'BJA94'.

# Installation

Requires scipy and numpy:

    sudo apt-get install python3-scipy

Get setuptools and install using setup.py:

    sudo apt-get install python3-setuptools
    python3 setup.py install --user

To have changes to the source reflected immediately:

    python3 setup.py develop --user

# Acknowledgement

The implementation of the tetrahedron method in [Quantum ESPRESSO](http://www.quantum-espresso.org/)
was consulted during the development of this implementation in order
to clarify understanding of the method and to help verify correctness.
