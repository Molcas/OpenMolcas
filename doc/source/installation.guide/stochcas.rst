.. _sec\:StochCAS_installation:

Installation of |molcas| for Stochastic-CASSCF calculations
===========================================================

.. only:: html

  .. contents::
     :local:
     :backlinks: none

The Stochastic-CASSCF method is based on the interface of the RASSCF program of |molcas|,
responsible for the orbital rotations via Super-CI, and the :program:`NECI` program,
responsible for the FCIQMC dynamics, replacing the deterministic Direct-CI based algorithm for large active space selection.
In principle, two installation protocols can be adopted that are referred to as embedded and uncoupled
forms. In the embedded form, the :program:`NECI` program is treated as a dependent subroutine of the
RASSCF program. This form effectively leads to an automatized version of the
Stochastic-CASSCF within the |openmolcas| software.
In the uncoupled form of Stochastic-CASSCF, :program:`NECI` is installed as a stand-alone program
and the |molcas|-:program:`NECI` interface is controlled manually by the user.
In this guide the uncoupled form will be discussed. It is the form preferred by
the developers of the method due to the non-black-box nature of the approach.

Uncoupled Form
++++++++++++++

The necessary routines in |molcas| are installed automatically,
but for improved communication between |molcas| and :program:`NECI` it is
recommended to compile with HDF5.

The :program:`NECI` code is available at https://github.com/ghb24/NECI_STABLE.

The :program:`NECI` code requires some external software and libraries:

* MPI: For builds intended to be run in parallel. OpenMPI, MPICH2 and its derivatives (IBM MPI, Cray MPI, and Intel MPI) have been tested.
* Linear algebra: ACML, MKL, BLAS/LAPACK combination.
* HDF5: To make use of the structured HDF5 format for reading/writing POPSFILES
   (files storing the population of walkers, and other information, to restart calculations).
   This library should be built with MPI and fortran support.

For configuring and compiling :program:`NECI` cmake is recommended::

  cmake -DENABLE_BUILD_HDF5=ON -DENABLE_HDF5=ON -DCMAKE_BUILD_TYPE=Cluster ~/neci/
  make -j hdf5
  make -j neci dneci

Cmake flag ``-DENABLE_BUILD_HDF5=ON`` builds the HDF5 library from source, and use that instead of one provided by the system.
Cmake flag ``-DENABLE_HDF5=ON`` makes use of HDF5 for popsfiles (default=on).

Two executable files will be generated: :file:`neci.exe` and :file:`dneci.exe`. The latter is compulsory for sampling one- and two-body
density matrices necessary for performing the orbital optimization. For a more detailed description of the :program:`NECI` configuration
the users are invited to read the available :program:`NECI` documentation.

More details about configuration/installation of the :program:`NECI` code can be found in the :program:`NECI` documentation.

There are currently no default verification tests for the Stochastic-CASSCF method. However, after installation of |molcas| one test is possible
to verify that MO integrals are correctly dumped into the FCIDUMP file. Simply use: ::

  molcas verify limannig


Embedded form
+++++++++++++

For the embedded form the :program:`NECI` source code has to be downloaded into the
|molcas| source directory.
Just execute in the |molcas| repository::

   git submodule update --init External/NECI

Then compile |molcas| with the ``-DNECI=ON`` cmake flag.
