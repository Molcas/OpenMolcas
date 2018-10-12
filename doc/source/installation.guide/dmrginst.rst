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
It requires the following components:

* HDF5: https://support.hdfgroup.org/HDF5
* CheMPS2 binary: https://github.com/SebWouters/CheMPS2/archive/v1.8.4.tar.gz

HDF5 must be configured with Fortran, Fortran 2003, and C++ enabled:

::

  ./configure --enable-fortran --enable-fortran2003 --enable-cxx

Make sure that both :file:`libhdf5_fortran.so` and :file:`libhdf5.so` are generated.

For CheMPS2 binary installation, consult
http://sebwouters.github.io/CheMPS2/sourcecode.html

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

  molcas verify 850 851
  molcas verify benchmark:970
