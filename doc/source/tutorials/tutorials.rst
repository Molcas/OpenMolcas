.. _TUT\:sec\:pg-based-tut:

Program Based Tutorials
=======================

The |molcas| |molcasversion| suite of Quantum Chemical programs is modular in
design. The desired calculation is achieved by executing a list of
|molcas| program modules in succession, while potentially manipulating
the program information files. If the information files from a previous
calculation are saved, then a subsequent calculation need not recompute
them. This is dependent on the correct information being preserved in
the information files for subsequent calculations.
Each module has keywords to specify the
functions to be carried out, and many modules rely on the
specification of keywords in previous modules.

The following sections describe the use of the |molcas| modules and
their inter-relationships. Each module is introduced in the
approximate order for performing a typical calculation.
A complete flowchart for the |molcas| |molcasversion| suite of programs follows.

.. toctree::

   flowchart
   tut_emil

   tut_gateway
   tut_seward
   tut_scf
   tut_mbpt2
   tut_rasscf
   tut_caspt2
   tut_rassi
   tut_casvb
   tut_motra
   tut_guga
   tut_mrci
   tut_cpf
   tut_ccsdt
   tut_optimiza
   tut_mckinley
   tut_mclr
   tut_genano
   tut_ffpt
   tut_vibrot
   tut_single_aniso
   tut_poly_aniso
   tut_grid_it
   tut_molden
   tut_expbas

   tut_errors
   tut_hints

.. tut_alaska
