.. _sec\:StochCAS_installation:

Installation of |molcas| for Stochastic-CASSCF calculations
===========================================================

.. only:: html

  .. contents::
     :local:
     :backlinks: none

The Stochastic-CASSCF method is based on the interface of the RASSCF program of |molcas|,
responsible for the orbital rotations via Super-CI, and the NECI program,
responsible for the FCIQMC dynamics, replacing the deterministic Direct-CI based algorighm for large active space selection.
In principle, two installation protocols can be adopted that are referred to as embedded and uncoupled
forms. In the embedded form, the NECI program is treated as a dependent subroutine of the
RASSCF program. This form effectively leads to an automatized version of the
Stochastic-CASSCF within the OpenMolcas software.
In the uncoupled form of Stochastic-CASSCF, NECI is installed as a stand-alone program
and the |molcas|-NECI interface is controlled manually by the user. In this guide only
the installation of the uncloupled form will be discussed as it is the form preferred by
the developers of the method due to the non-black-box nature of the approach.
In order to configure the uncoupled form of the Stochastic-CASSCF method
simply use the cmake flag -DNECI=ON::

  cmake -DNECI=ON ~/molcas/

The NECI code is available at https://github.com/ghb24/NECI_STABLE. 

The code requires some external software and libraries:

* MPI: For builds intended to be run in parallel. OpenMPI, MPICH2 and its derivatives (IBM MPI, Cray MPI, and Intel MPI) have been tested.
* Linear algebra: ACML, MKL, BLAS/LAPACK combination.
* HDF5: To make use of the structured HDF5 format for reading/writing POPSFILES (files storing the population of walkers, and other information, to restart calculations). This library should be built with MPI and fortran support.

For configuring and compiling NECI cmake is recommended::

  cmake -DENABLE_BUILD_HDF5=ON -DENABLE_HDF5=ON -DCMAKE_BUILD_TYPE=Cluster ~/neci/
  make -j hdf5
  make -j neci dneci

Cmake flag -DENABLE_BUILD_HDF5=ON builds the hdf5 library from source, and use that instead of one provided by the system.
Cmake flag -DENABLE_HDF5=ON makes use of hdf5 for popsfiles (default=on).

Two executable files will be generated: neci.exe and dneci.exe. The latter is compulsory for sampling one- and two-body
density matrices necessary for performing the orbital optimization. For a more detailed description of the NECI configuration
the users are invited to read the available NECI documentation.

More details about configuration/installation of the NECI code can be found in the NECI documentation.

There are currently no verification tests for the Stochastic-CASSCF method. However, after installation of Molcas one test is possible
to verify that MO integrals are correctly dumped into the FCIDUMP file. Simply use
::

  molcas verify limannig

