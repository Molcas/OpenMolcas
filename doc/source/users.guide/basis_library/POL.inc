Polarized basis sets
....................

The so-called polarized basis sets are purpose oriented, relatively
small GTO/CGTO sets devised for the purpose of accurate calculations
of dipole electric properties of polyatomic molecules
:cite:`polI,polII,polIII,polIV,polV`.
For each row of the periodic table the performance of the basis sets
has been carefully examined in calculations of dipole moments and
dipole polarizabilities of simple hydrides at both the SCF and
correlated levels of approximation :cite:`polI,polII,polIII,polIV,polV`.
The corresponding results match within a few percent the best
available experimental data. Also the calculated molecular quadrupole
moments turn out to be fairly close to those computed with much larger
basis sets. According to the present documentation the polarized basis
GTO/CGTO sets can be used for safe accurate predictions of molecular
dipole moments, dipole polarizabilities, and also molecular quadrupole
moments by using high-level correlated computational methods. The use
of the polarized basis sets has also been investigated in calculations
of weak intermolecular interactions. The interaction energies,
corrected for the basis set superposition effect (BSSE), which is
rather large for these basis sets, turn out to be close to the best
available data. In calculations for molecules involving the 4th row
atoms, the property data need to be corrected for the relativistic
contribution. The corresponding finite perturbation facility is
available :cite:`polVI,polVII`.

It is recommended to use these basis sets with the contraction given
in the library. It is of course possible to truncate them further,
for example by deleting some polarization functions, but this will
lead to a deterioration of the computed properties.

.. Below is a list of the basis sets available in this class. ::

     H.Pol.Sadlej.6s4p.3s2p.
     Li.Pol.Sadlej.10s6p4d.5s3p2d.
     Be.Pol.Sadlej.10s6p4d.5s3p2d.
     C.Pol.Sadlej.10s6p4d.5s3p2d.
     N.Pol.Sadlej.10s6p4d.5s3p2d.
     O.Pol.Sadlej.10s6p4d.5s3p2d.
     F.Pol.Sadlej.10s6p4d.5s3p2d.
     Na.Pol.Sadlej.14s10p4d.7s5p2d.
     Mg.Pol.Sadlej.14s10p4d.7s5p2d.
     Si.Pol.Sadlej.14s10p4d.7s5p2d.
     P.Pol.Sadlej.14s10p4d.7s5p2d.
     S.Pol.Sadlej.14s10p4d.7s5p2d.
     Cl.Pol.Sadlej.14s10p4d.7s5p2d.
     K.Pol.Sadlej.15s13p4d.9s7p2d.
     Ca.Pol.Sadlej.15s13p4d.9s7p2d.
     Ge.Pol.Sadlej.15s12p9d.9s7p4d.
     As.Pol.Sadlej.15s12p9d.9s7p4d.
     Se.Pol.Sadlej.15s12p9d.9s7p4d.
     Br.Pol.Sadlej.15s12p9d.9s7p4d.
     Rb.Pol.Sadlej.18s15p10d.11s9p4d.
     Sr.Pol.Sadlej.18s15p10d.11s9p4d.
     Sn.Pol.Sadlej.19s15p12d4f.11s9p6d2f.
     Sb.Pol.Sadlej.19s15p12d4f.11s9p6d2f.
     Te.Pol.Sadlej.19s15p12d4f.11s9p6d2f.
     I.Pol.Sadlej.19s15p12d4f.11s9p6d2f.
