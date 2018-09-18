.. index::
   single: Program; Seward
   single: Seward
   single: Integrals

.. _TUT\:sec\:seward:

:program:`SEWARD` --- An Integral Generation Program
====================================================

An *ab initio* calculation always requires integrals. In the
|molcas| suite of programs, this function is supplied by the :program:`SEWARD`
module. :program:`SEWARD` computes the one- and two-electron integrals for the
molecule and basis set specified in the input to the program :program:`GATEWAY`,
which should be run before :program:`SEWARD`. :program:`SEWARD` can also be used
to perform some property expectation calculations on the isolated molecule.
The module is also used as an input parser for the reaction field and
numerical quadrature parameters.

We commence our tutorial by calculating the integrals for a water molecule. The
input is given in :numref:`block:seward_input`. Each |molcas| module
identifies input from a file by the name of the module. In the case of
:program:`SEWARD`, the program starts with the label
``&SEWARD``, which is the first statement in the file shown below.

In normal cases no input is required for :program:`SEWARD`, so the following
input is optional. The first keyword used is :kword:`TITLe`. Only the first
line of the title is printed in the output. The first title line is also saved
in the integral file and appears in any subsequent programs that use the
integrals calculated by :program:`SEWARD`.

.. code-block:: none
   :caption: Sample input requesting the :program:`SEWARD` module
              to calculate the integrals for water in :math:`C_{2v}` symmetry.
   :name: block:seward_input

   &SEWARD
   Title
   Water - A Tutorial. The integrals of water are calculated using C2v symmetry

.. Sample input requesting the :program:`SEWARD` module
   to calculate the integrals for water in :math:`C_{2v}` symmetry.

In more complicated cases more input may be needed, to specify certain types of
integrals, that use of Cholesky decomposition techniques (:kword:`CHOLesky` keyword), etc. We refer to the
specific sections of the User's Guide for more information.
The output from a :program:`SEWARD` calculation is small and contains in principle
only a list of the different types of integrals that are computed.

The integrals produced by the :program:`SEWARD` module are stored in
two files in the working directory. They are ascribed the :program:`FORTRAN`
names :file:`ONEINT` and :file:`ORDINT` which are
automatically symbolically linked by the |molcas| script to the file
names :variable:`$Project`:file:`.OneInt` and
:variable:`$Project`:file:`.OrdInt`, respectively
or more specifically, in our case, :file:`water.OneInt` and
:file:`water.OrdInt`, respectively. The default name for each
symbolical name is contained in the corresponding program files of the
directory :file:`$MOLCAS/shell`.
The :file:`ONEINT` file contains the one-electron integrals.
The :file:`ORDINT` contains the ordered and packed two-electron integrals.
Both files are used by later |molcas| program modules.

.. :program:`SEWARD` --- Basic and Most Common Keywords
   ----------------------------------------------------

   .. class:: keywordlist

   :kword:`Cholesky`
     Use Cholesky decomposition

   :kword:`AMFI`
     Atomic mean-field integrals for relativistic calculations.
     Required for spin-coupling. Automatic for ANO-RCC basis sets
