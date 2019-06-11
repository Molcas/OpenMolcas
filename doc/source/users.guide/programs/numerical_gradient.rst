.. index::
   single: Program; Numerical_Gradient
   single: Numerical_Gradient

.. _UG\:sec\:numerical_gradient:

:program:`numerical_gradient`
=============================

.. only:: html

  .. contents::
     :local:
     :backlinks: none

.. xmldoc:: %%Description:
            The Numerical_Gradient module is a program which numerically evaluates the
            gradient of the energy with respect to nuclear perturbations.

The :program:`Numerical_Gradient` module is a program which numerically evaluates the gradient
of the energy with respect to nuclear perturbations.

Note that this module is automatically invoked by the :program:`Alaska` module if the wave function
method is MBPT2, CCSDT, CASPT2, MS-CASPT2, or a calculation using the Cholesky decomposition.
The user should normally never request the execution of this module; instead it is advised to use the
:kword:`NUMErical` keyword in Alaska, if it is necessary to force the use of numerical gradients rather than
analytical ones.

The module is parallelized over the displacements, which in case of large jobs gives a linear
speed up compared to a serial execution, although in order to obtain this it is important to
choose the number of nodes such that the number of contributing perturbations is a multiple of
the number of nodes. For a given molecule the number of perturbations equals the number of atoms
times 6 (a perturbation with plus and minus delta for each of the three axes). Symmetry can of
course reduce this number. If the request of execution originates from the :program:`Slapaf`
module further reduction in perturbations is achieved due to the utilization of rotational and
translational invariance.

.. index::
   pair: Dependencies; Numerical_Gradient

.. _UG\:sec\:numerical_gradient_dependencies:

Dependencies
------------

The dependencies of the :program:`Numerical_Gradient` module is the union
of the dependencies of the :program:`SEWARD`, :program:`SCF`, :program:`RASSCF`,
:program:`MBPT2`, :program:`MOTRA`, :program:`CCSDT`, and
:program:`CASPT2` modules.

.. index::
   pair: Files; Numerical_Gradient

.. _UG\:sec\:numerical_gradient_files:

Files
-----

The files of the :program:`Numerical_Gradient` module is the union
of the files of the :program:`SEWARD`, :program:`SCF`, :program:`RASSCF`,
:program:`MBPT2`, :program:`MOTRA`, :program:`CCSDT`, and
:program:`CASPT2` modules.
