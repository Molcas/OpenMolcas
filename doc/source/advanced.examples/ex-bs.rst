Extra information about basis sets and integrals
================================================

.. only:: html

  .. contents::
     :local:
     :backlinks: none

Basis set format
----------------

.. compound::

  The **Inline** option for a basis set will read the basis set
  as defined by the following pseudo code. ::

    Read Charge, lAng
    Do iAng = 0, lAng
       Read nPrim, nContr
       Read (Exp(iPrim),iPrim=1,nPrim)
       Do iPrim=1,nPrim
          Read (Coeff(iPrim,iContr),iContr=1,nContr)
       End Do
    End Do

  where ``Charge`` is the nuclear charge, ``lAng`` is the highest angular
  momentum quantum number, ``nPrim`` is the number of primitive functions
  (exponents) for a given shell, and ``nContr`` is the number of contracted
  functions for a given shell.

The following is an example of a DPZ basis set for carbon.
Normally, however, the basis set will be read from a library file following
the specified label (like, e.g., C.DZP...4s2p1d.),
and not be inserted inline at the input file. ::

  Basis set                                -- Start defining a basis set
  C.DZP.Someone.9s5p1d.4s2p1d. / inline    -- Definition in input stream
         6.0              2                -- charge, max l-quantum no.
      9    4                               -- no. of prim. and contr. s-functions
  4232.61                                  -- s-exponents
  634.882
  146.097
  42.4974
  14.1892
  1.9666
  5.1477
  0.4962
  0.1533
    .002029   .0       .0       .0         -- s-contraction matrix
    .015535   .0       .0       .0
    .075411   .0       .0       .0
    .257121   .0       .0       .0
    .596555   .0       .0       .0
    .242517   .0       .0       .0
    .0       1.0       .0       .0
    .0        .0      1.0       .0
    .0        .0       .0      1.0
      5    2                               -- no. of prim. and contr. p-functions
  18.1557                                  -- p-exponents
  3.98640
  1.14290
  0.3594
  0.1146
   .018534   .0                            -- p-contraction matrix
   .115442   .0
   .386206   .0
   .640089   .0
   .0       1.0
      1    1                               -- no. of prim. and contr. d-functions
     .75                                   -- d-exponents
    1.0                                    -- d-contraction matrix
  C1 0.00000 0.00000 0.00000               -- atom-label, Cartesian coordinates
  C2 1.00000 0.00000 0.00000               -- atom-label, Cartesian coordinates
  End Of Basis                             -- end of basis set definition

The basis set label and the ECP libraries
-----------------------------------------

The label within the :file:`ECP` library
is given as input in the line following the keyword :kword:`BASIS SET`.
The label defines either the valence basis set and core potential
which is assigned to a
frozen-core atom
or the embedding potential
which is assigned to an environmental froze-ion.
Here, all the comments made about this label in the section
**The basis set label and the basis set library**
for all-electron basis sets
stand, except for the following changes:

#. The identifier ``type`` must be ``ECP`` or ``PP``.

#. The identifier ``aux`` specifies the kind of the potential.
   It is used, for instance, to choose between non-relativistic, Cowan--Griffin, or no-pair
   Douglas--Kroll relativistic core potentials
   (i.e. ``Pt.ECP.Casarrubios.13s10p9d5f.1s2p2d1f.16e-NR-AIMP.``
   or ``Pt.ECP.Casarrubios.13s10p9d5f.1s2p2d1f.16e-CG-AIMP.``
   or ``Pt.ECP.Rakowitz.13s10p9d6f.5s4p4d2f.18e-NP-AIMP.``)
   and to pick up one among all the embedding potentials available
   for a given ion
   (i.e. ``F.ECP.Lopez-Moraza.0s.0s.0e-AIMP-KMgF3.``
   or ``F.ECP.Lopez-Moraza.0s.0s.0e-AIMP-CsCaF3.``).

#. The identifier ``contracted`` is used here
   in order to produce the actual basis set
   out of the basis set included in the :file:`ECP` library,
   which is a minimal basis set (in general contraction form) augmented
   with some polarization, diffuse, ... function.
   It indicates the number of s, p, ... contracted functions
   in the actual basis set,
   the result being always a many-primitive contracted function
   followed by a number of primitives.
   As an example,
   ``At.ECP.Barandiaran.13s12p8d5f.3s4p3d2f.17e-CG-AIMP.``
   will generate a (13,1,1/12,1,1,1/8,1,1/5,1) formal contraction pattern
   which is in this case a (13,1,1/12,1,1,1/7,1,1/5,1) real pattern.
   Other contraction patters should be input "Inline".

#. The user is suggested to read carefully :numref:`TUT:sec:ecp`
   of the tutorials and examples manual before using the ECP utilities.

.. _UG\:sec\:one-electron_integral_labels:

One-Electron Integral Labels
----------------------------

.. compound::

  The storage of one-electron integrals on disk is facilitated by the
  one-electron integral I/O facility. The internal structure of the
  one-electron file and the management is something which the user normally
  do not need to worry about. However, for the general input section of the
  :program:`FFPT`, the user need to know the name and structure of the internal
  labels which the one-electron integral I/O facility associates with each type
  of one-electron integral. The labels are listed and explained here below for reference.
  The component index is also utilized by the one-electron integral I/O facility to
  discriminate the various components of the one-electron integrals of a certain type,
  for example, the dipole moment integrals have three components (1=x-component,
  2=y-component, 3=z-component). The component index is enumerated as a canonical
  index over the powers of the Cartesian components of the operator (e.g. multipole
  moment, velocity, electric field, etc.). The order is defined by following pseudo
  code, ::

    Do ix = nOrder, 0, -1
       Do iy = nOrder-ix, 0, -1
          iz = nOrder-ix-iy
       End Do
    End Do,

  where ``nOrder`` is the total order of the operator, for example, ``nOrder=2`` for
  the electric field gradient and the quadrupole moment operator.

.. _tab\:bs:

============ =========================================================================
Label        Explanation
============ =========================================================================
``Mltpl nn`` the ``nn``\ th order Cartesian multipole moments.
``MltplS``   the overlap matrix used in the semi-empirical NDDO method.
``Kinetic``  the kinetic energy integrals.
``Attract``  the electron attraction integrals.
``AttractS`` the electron attraction integrals used in the semi-empirical NDDO method.
``PrjInt``   the projection integrals used in ECP calculations.
``M1Int``    the M1 integrals used in ECP calculations.
``M2Int``    the M2 integrals used in ECP calculations.
``SROInt``   the spectrally resolved operator integrals used in ECP calculations.
``XFdInt``   the external electric field integrals.
``MassVel``  the mass-velocity integrals.
``Darwin``   the Darwin one-electron contact integrals.
``Velocity`` the velocity integrals.
``EF0nnnnn`` the electric potential at center ``nnnnn``.
``EF1nnnnn`` the electric field at center ``nnnnn``.
``EF2nnnnn`` the electric field gradient at center ``nnnnn``.
``AngMom``   the angular momentum integrals.
``DMS``      the diamagnetic shielding integrals.
``Wellnnnn`` the ``nnnn``\ th set of spherical well integrals.
``OneHam``   the one-electron Hamiltonian.
``AMProd``   the hermitized product of angular momentum integrals.
``AMFI``     the atomic mean field integrals.
============ =========================================================================
