This Molcas utility adds HDF5 functionality to Molcas.  In short, it
provides a runfile-like functionality, but for any number of files
(several HDF5 files can be active at the same time) and without the
need to pre-define the data that is in there.

This readme is intended to give a quick overview of the functionality
and reasons to use this library, but for the full information on HDF5
you can look here:

  https://www.hdfgroup.org/HDF5/


What exactly is HDF5?

HDF5 stands for Hierarchical Data Format, and is a format used to
store data in binary format in a portable way. The HDF5 library
provides I/O functionality for dealing with HDF5 files, and provides a
set of APIs for multiple programming languages. The library also comes
with several tools for manipulating and viewing the data.

Why choose HDF5?

The properties I was looking for were portability (being able to copy
binary files to any machine and use them with Molcas) and ease of
accessing the data from outside Molcas. The other candidate was
NetCDF, but since they have HDF5 anyway as a dependency, I chose to go
with the least amount of extra dependencies.

HDF5 has a Fortran 90 API, why not use those?

The Fortran 90 API for HDF5 has one big problem: it relies on the use
of Fortran modules. This means that in practice, a user would have to
compile the HDF5 library themselves with the same compiler (and
version) as they use for Molcas. Even on my own laptop, where HDF5
comes preinstalled, the package had been compiled with a slightly
older version of gfortran and resulted in a problem with the modules.

The general advice is to not use Fortran modules for external libraries,
but people at HDF5 apparently didn't know or bother with this. Even though
I found that on the computer clusters I used, HDF5 is compiled for each
available Fortran compiler, the need for maximum portability to user's
computers made me decide to use the C API, and write only Fortran wrappers
for the convenience functions used in Molcas
