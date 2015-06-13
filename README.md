An implementation of the tetrahedron method for Brillouin zone summation, as
described in: [Blöchl, Jepsen, and Andersen, PRB 49, 16223 (1994)](http://journals.aps.org/prb/abstract/10.1103/PhysRevB.49.16223).
This paper is referred to in the documentation of this package as 'BJA94'.

# Installation

Get setuptools and install using setup.py:

    sudo apt-get install python3-setuptools
    sudo python3 setup.py install

To have changes to the source reflected immediately:

    sudo python3 setup.py develop

# Acknowledgement

The implementation of the tetrahedron method in [Quantum ESPRESSO](http://www.quantum-espresso.org/)
was consulted during the development of this implementation in order
to clarify understanding of the method and to help verify correctness.
