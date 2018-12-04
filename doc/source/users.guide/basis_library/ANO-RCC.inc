Relativistic ANO basis sets --- ANO-RCC
.......................................

Extended relativistic ANO-type basis sets are available for the atoms H--Cm.
These basis sets have been generated using the same principles as described
above for the ANO-L basis sets with the difference that the density matrices
have been computed using the CASSCF/CASPT2 method. The basis have been
contracted using the Douglas--Kroll Hamiltonian and should
therefore only be used in calculations where scalar relativistic effects are
included. Seward will automatically recognize this and turn on the DK option
when these basis sets are used :cite:`Roos:03c,Roos:03g,Roos:05a,Roos:05b`. The
basis sets contain functions for correlation of the semi-core electrons. The
new basis sets are called ANO-RCC. More details about the construction and
performance is given in the header for each basis set in the ANO-RCC library.
Basis sets are available for all atoms up to Cm.

Scalar relativistic effect become important already in the second row of the
periodic systems. It is therefore recommended to use these basis sets instead of
ANO-L in all calculations.

For information about the primitive basis set we refer to the library.
The maximum number of ANOs given in the library is:

* 6s4p3d1f for H.
* 7s4p3d2f for He.
* 8s7p4d2f1g for Li--Be.
* 8s7p4d3f2g for Be--Ne.
* 17s12p5d4f for Na.
* 9s8p5d4f for Mg--Al.
* 8s7p5d4f2g for Si--Ar.
* 10s9p5d3f for K.
* 10s9p6d2f for Ca.
* 10s10p8d6f4g2h for Sc--Zn.
* 9s8p6d4f2g for Ga--Kr.
* 10s10p5d4f for Rb--Sr.
* 10s9p8d5f3g for In--Xe.
* 12s10p8d4f for Cs--Ba.
* 11s10p8d5f3g for La.
* 12s11p8d7f4g2h for Ce--Lu.
* 11s10p9d8f4g2h for Hf--Hg.
* 11s10p9d6f4g for Tl--Rn.
* 12s11p8d5f for Fr--Ra.
* 13s11p10d8f6g3h for Ac--Pa.
* 12s10p9d7f5g3h for U--Cm.

However, such contractions are unnecessarily large. Almost converged
results (compared to the primitive sets) are usually obtained with basis
sets of QZP quality. You can get a feeling for the convergence from the test
results presented in the header of each basis set in the library. One should
also remember that larger basis sets are needed for the correlation of semi-core
electrons.

Below is a list of the core electrons correlated for each atom.

======= ===================
Li--B:  1s
C--Ne:  No core correlation
Na:     2s,2p
Mg--Al: 2p
Si--Ar: No core correlation
K:      3s,3p
Ca--Zn: 3p
Ga--Ge: 3d
As--Kr: No core correlation
Rb--Sr: 4p
In--Xe: 4d
Cs--Ba: 5p
La--Lu: 5s,5p
Hf--Re: 4f,5s,5p
Os--Hg: 5s,5p
Tl--Rn: 5d
Fr--Ra: 6p
Ac--Cm: 6s,6p
======= ===================

Basis set label in input:

The general label is given as for the other ANO basis sets:

  *Atom.ANO-RCC...contracted set.* (Note the last dot!).

A short hand notation is also possible:

  *Atom.ANO-RCC-label*, where *label* is one of MB, VDZ, VDZP, VTZP, or VQZP.

A translation between the two possibilities can be found in file:
:file:`$MOLCAS/basis_library/basis.tbl`
