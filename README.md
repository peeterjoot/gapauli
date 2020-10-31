# gapauli
A geometric algebra Mathematica module for Euclidean 2D and 3D spaces using the Pauli basis representation as a backend, and a Mathematica module for STA (relativistic space with a +,-,-,- signature) using 

Three modules are provided:

* GA20
* GA30
* GA13

Which implement GA(2,0), GA(3,0), and GA(1,3) algebras respectively.  The first is implemented using only real Pauli matrixes
(\sigma_1, and \sigma_3), the second uses all three Pauli matrices, and the last uses Dirac matrices.

Three generic test notebooks are provided, each of which also provides some documentation

* testGA20.nb   Test cases and documentation for GA20.m (online (w/o Mathematica) viewable version saved as testGA20.pdf)
* testGA30.nb   Test cases and documentation for GA30.m (online (w/o Mathematica) viewable version saved as testGA30.pdf)
* testGA13.nb   Test cases and documentation for GA13.m (online (w/o Mathematica) viewable version saved as testGA13.pdf)

Some other ad-hoc demonstrations are also available:

* bivectorCommutatorGA13.nb -- symbolic verification that <AB>_{0,4} = (1/2(AB + BA), and <AB>_{2} = (1/2(AB - BA) for STA bivectors A,B.
* lorentzForce.nb -- The Force and power components of the covariant (STA) Lorentz force equation are expanded symbolically.
* shortCurrentFilament.nb -- don't remember.
* ellipticParameterization.nb -- don't remember.

