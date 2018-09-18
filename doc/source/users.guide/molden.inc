.. index::
   single: Tool; LUSCUS
   single: LUSCUS

.. _UG\:sec\:Molden:

Writing LUSCUS/MOLDEN input
---------------------------

By default the :program:`GUESSORB`, :program:`SCF`, :program:`MBPT2`, :program:`RASSCF`,
:program:`SLAPAF`, :program:`LOCALISATION`, and :program:`MCLR` modules
generate input in Molden format. The :program:`SCF`, :program:`MBPT2`, :program:`RASSCF`,
and :program:`LOCALISATION` modules generate input for molecular orbital
analysis, :program:`SLAPAF` for geometry optimization analysis, minimum energy paths,
Saddle optimization paths and IRC TS analysis,
and the :program:`MCLR` module generates input for
analysis of harmonic frequencies.

The generic name of the input file and the actual
name are different for the nodes as a reflection on the data generated
by each module. Hence, the actual names (generic name) for the Molden files in each module are

* GUESSORB module:
  $Project.guessorb.molden (MD_GSS)
* SCF module:
  $Project.scf.molden (MD_SCF)
* MBPT2 module:
  $Project.mp2.molden (MD_MP2)
* RASSCF module:
  $Project.rasscf.molden (MD_CAS) for the state-averaged natural orbitals, and
  $Project.rasscf.x.molden (MD_CAS.x) for the state-specific natural spin orbitals,
  where x is the index of a CI root.
* SLAPAF module:
  $Project.geo.molden (MD_GEO) for geometry optimizations,
  $Project.mep.molden (MD_MEP) for minimum energy paths,
  $Project.irc.molden (MD_IRC) for IRC analysis of a TS, and
  $Project.saddle.molden (MD_SADDLE) for Saddle method TS optimizations.
* LOCALISATION module:
  $Project.local.molden (MD_LOC)
* MCLR module:
  $Project.freq.molden (MD_FREQ)

Files in MOLDEN format can be visualized by Molden program,
or by |molcas| GUI program :program:`LUSCUS`.

Note, the MOLDEN code cannot handle cases in which some shells, with the
same angular momentum,
are Cartesian and some are real spherical harmonics. All
shells for a give angular momentum have to be either Cartesian or
real spherical harmonics, otherwise no MOLDEN file for orbitals will be generated.

For further details with respect to Molden consult
http://www.cmbi.ru.nl/molden/.

For further details with respect to LUSCUS consult
http://luscus.sourceforge.net/.
