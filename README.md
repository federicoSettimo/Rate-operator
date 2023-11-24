# Rate Operator Quantum Jumps
This project runs the simulations using the ROQJ method.

The files `roqj.h` and `roqj.cpp` contain the library performing the unravelings. 
When used, it requires the definition in the `.cpp` file of the Hamiltonian `H`, the jump term of the master equation `J`, the anti-commutator part `Gamma`, the transformation `Phi` defining the RO, and the `observable` calculated. For the linear algebra, it requires [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page).
