.. _sec\:chemps2_installation:

Installation of CheMPS2--|molcas| interface for DMRG calculations
=================================================================

.. only:: html

  .. contents::
     :local:
     :backlinks: none

The CheMPS2--|molcas| interface requires the following components:
The CheMPS2--|molcas| interface :cite:`Phung2016,Wouters2016`,
based on the Block--|molcas| interface :cite:`Nakatani2017`,
can support DMRG-SS-CASPT2 and DMRG-SA-CASPT2 calculations.

It requires the CheMPS2 binary. For installation of CheMPS2, consult
http://sebwouters.github.io/CheMPS2/index.html if it is not already
available in your OS.

Note that only the version with the Open Multi-Processing (OpenMP) is supported,
thus build CheMPS2 with:

::

  -D WITH_MPI=OFF

In order to efficiently run the CheMPS2--|molcas| interface,
it is advisible to compile either serial or parallel |molcas| with MPI.
An example:

::

  ./configure -compiler intel -parallel -64 -mpiroot /path/to/mpi/root \
              -mpirun /path/to/mpi/bin/mpirun -blas MKL -blas_lib -mkl=sequential \
              -hdf5_inc /path/to/hdf5/include \
              -hdf5_lib /path/to/hdf5/lib \
              -chemps2 /path/to/chemps2/binary

The CheMPS2--|molcas| interface can also be activated with CMake:

::

  -D CHEMPS2=ON -D CHEMPS2_DIR=/path/to/chemps2/binary

Before testing the CheMPS2--|molcas| interface, make sure to increase stack size,
export :variable:`OMP_NUM_THREADS`, the CheMPS2 binary, and all the required libraries for CheMPS2.

::

  ulimit -s unlimited
  [export OMP_NUM_THREADS=...]
  export PATH=/path/to/chemps2/binary:$PATH

Verify the installation:

::

  molcas verify extra:850,851
  molcas verify benchmark:970
