.. index::
   single: Tool; MOLDEN
   single: MOLDEN

.. _TUT\:sec\:Molden:

Writing MOLDEN input
====================

By default the :program:`GUESSORB`, :program:`SCF`, :program:`MBPT2`, :program:`RASSCF`,
:program:`SLAPAF`, :program:`LOCALISATION`, and :program:`MCLR` modules
generate input in Molden format. The :program:`SCF`, :program:`MBPT2`, :program:`RASSCF`,
and :program:`LOCALISATION` modules generate input for molecular orbital
analysis, :program:`SLAPAF` for geometry optimization analysis, minimum energy paths,
Saddle optimization paths and IRC TS analysis,
and the :program:`MCLR` module generates input for
analysis of harmonic frequencies. Molden files can be visualized by :program:`GV`
or by :program:`Molden` (http://www.cmbi.ru.nl/molden/).

The generic name of the input file and the actual
name are different for the nodes as a reflection on the data generated
by each module. Hence, the actual names (generic name) for the Molden files in each module are

* :program:`GUESSORB` module:
  :file:`$Project.guessorb.molden` (:file:`MD_GSS`)
* :program:`SCF` module:
  :file:`$Project.scf.molden` (:file:`MD_SCF`)
* :program:`MBPT2` module:
  :file:`$Project.mp2.molden` (:file:`MD_MP2`)
* :program:`RASSCF` module:
  :file:`$Project.rasscf.molden` (:file:`MD_CAS`) for the state-averaged natural orbitals, and
  :file:`$Project.rasscf.x.molden` (:file:`MD_CAS.x`) for the state-specific natural spin orbitals,
  where :file:`x` is the index of a CI root.
* :program:`SLAPAF` module:
  :file:`$Project.geo.molden` (:file:`MD_GEO`) for geometry optimizations,
  :file:`$Project.mep.molden` (:file:`MD_MEP`) for minimum energy paths,
  :file:`$Project.irc.molden` (:file:`MD_IRC`) for IRC analysis of a TS, and
  :file:`$Project.saddle.molden` (:file:`MD_SADDLE`) for Saddle method TS optimizations.
* :program:`LOCALISATION` module:
  :file:`$Project.local.molden` (:file:`MD_LOC`)
* :program:`MCLR` module:
  :file:`$Project.freq.molden` (:file:`MD_FREQ`)
