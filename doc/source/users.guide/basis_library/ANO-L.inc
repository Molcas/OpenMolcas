Large ANO basis sets -- ANO-L
.............................

The large ANO basis sets for atoms H--Zn, excluding K and Ca,
have been constructed by averaging the
corresponding density matrix over several atomic states, positive and
negative ions and the atom in an external
electric field :cite:`anoI,anoII,anoIII`.
The different density matrices have been obtained from correlated
atomic wave functions. Usually the SDCI method has been used. The
exponents of the primitive basis have in some cases been optimized.
The contracted basis sets give virtually identical results as the
corresponding uncontracted basis sets for the atomic properties, which
they have been optimized to reproduce. The design objective has been
to describe the ionization potential, the electron affinity, and the
polarizability as accurately as possible. The result is a well
balanced basis set for molecular calculations.

For information about the primitive basis set we refer to the library.
The maximum number of ANO's given in the library is:

* 6s4p3d for H.
* 7s4p3d for He.
* 7s6p4d3f for Li--Be.
* 7s7p4d3f for B--Ne.
* 7s7p5d4f for Na--Ar.
* 8s7p6d5f4g for Sc--Zn.

However, such contractions are unnecessarily large. Almost converged
results (compared to the primitive sets) are obtained with the VQZP basis
sets:

* 3s2p1d for H--He.
* 5s4d3d2f for Li--Ne.
* 6s5p4d3f for Na--Ar.
* 7s6p5d4f3g for Sc--Zn.

The results become more approximate below the size:

* 3s2p for H--He.
* 4s3p2d for Li--Ne.
* 5s4p2d for Na--Ar.
* 6s5p4d3f for Sc--Zn.

It is recommended to use at least two polarization
(3d/4f) functions, since one of them is used for polarization and the
second for correlation. If only one 3d/4f-type function is used one has
to decide for which purpose and adjust the exponents and the
contraction correspondingly. Here both effects are described jointly
by the two first 3d/4f-type ANO's (The same is true for the hydrogen
2p-type ANO's). For further discussions regarding the use of these
basis sets we refer to the
literature :cite:`anoI,anoII,anoIII`.
