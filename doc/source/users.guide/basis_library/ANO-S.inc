Small ANO basis sets -- ANO-S
..............................

The smallest of the Atomic Natural Orbital (ANO) basis sets are available for
the atoms H--Kr. They have been constructed as eigenfunctions of a density
matrix averaged over several electronic configurations. The ground state
of the atom was included for all atoms, and dependent on the particular atom
one or more of the following states were included:
valence excited states, ground state for the anion and ground state for
the cation. The density matrices were obtained by the SCF, SDCI or MCPF
methods for 1 electron, 2 electron and many electron cases respectively.
The emphasis have been on obtaining good structural properties such as
bond-lengths and -strengths with as small contracted sets as possible.
The quality for electric properties such as polarizabilities have been
sacrificed for the benefit of the properties mentioned above.
See :cite:`anoIV` for further discussions. These basis sets are recommended for
large molecules where the more extended ANO-L basis sets require to much
computational times. One should, however, remember that for a given contraction
it is only the time needed to generate the integrals (or Cholesky vectors) that
is affected and it is usually preferred to use the more accurate ANO-L (or
ANO-RCC) basis sets.

For information about the primitive basis set we refer to the library.
The maximum number of ANO's given in the library is:

* 4s3p for H--He.
* 6s4p3d for Li--Be.
* 7s6p3d for B--Ne.
* 7s5p3d for Na--Mg.
* 7s7p4d for Al--Ar.
* 7s7p4d for K--Ca.
* 8s7p7d4f for Sc--Zn.
* 9s9p5d for Ga--Kr.

However, such contractions are unnecessarily large. Almost converged
results (compared to the primitive sets) are obtained with the basis
sets:

* 3s2p for H--He.
* 4s3p2d for Li--Ne.
* 5s4p3d for Na--Ar.
* 6s5p4d for K--Ca.
* 7s5p4d3f for Sc--Zn.
* 6s5p4d for Ga--Kr.

The results become more approximate below the DZP size:

* 2s1p for H--He.
* 3s2p1d for Li--Ne.
* 4s3p2d for Na--Ar.
* 5s4p3d for K--Ca.
* 6s4p3d2f for Sc--Zn.
* 5s4p3d for Ga--Kr.
