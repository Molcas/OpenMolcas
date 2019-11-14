Installation
============

.. only:: html

  .. contents::
     :local:
     :backlinks: none

The present installation guide describes the necessary steps for installing
and tailoring |molcas|. It also describes the steps for applying updates
whenever necessary.

The installation procedure can be reduced to a few simple steps:

#. Extract the contents of the tar
#. Configure the package
#. Build the package
#. Build GUI and documentation (optional)
#. Make the package generally available

Prerequisites
-------------

Prerequisite hardware
.....................

In general, |molcas| can be built on any hardware that runs under a UNIX operating system.
Some of these variants of hardware and software have been tested by us, and you
should not have any problems to install |molcas| on any of these.
For other platforms you will most likely need to put some extra effort into the installation.
In many cases the only effort on your part is setting some compiler flags,
paths to system software etc.
For a list of the platforms where we have
successfully installed |molcas| see our homepage:
|MolcasWWW|.

To load the executables resident, sufficient memory is required.
In addition, the programs are enabled to allocate work space dynamically.
To avoid excessive paging we recommend that your machine should be
equipped with at least 2 GB of memory per running application. Note, that
|molcas| will run faster with more memory.

To build |molcas| you may need up to 2 GB of free disk space depending on
what you choose to install and how (e.g. CMake keeps the object files around).
The actual size of a |molcas| installation is around 300-600 MB.
To run the verification tests of |molcas| you should have a scratch disk
with up to 1 GB of free disk space, depending on the suite you run. For the
"small" set about 400 MB will suffice.
To perform larger calculations, ample amount of scratch disk space is necessary.
The exact amount varies with the type of systems studied, but a general
recommendation is at least 4 GB of disk space, per production run.

Prerequisite software
.....................

If you obtain the source code of |molcas|, then you need to
make certain that the necessary software is available to build |molcas|.
The minimum requirements are:

* A Fortran compiler (with support for Fortran 2003 and later)
* A C compiler
* GNU make
  See URL |GnuWWW| and navigate to the gnumake page or go directly
  to |GnuMakeWWW|
* Perl (version 5.008 or higher)

Also, you can benefit from following optional dependencies:

* CMake (version 2.8.11 or higher, recommendeded for easier configuration)
* an optimized BLAS/LAPACK library (alternative for the slow built-in library based on Netlib)
* MPI-2 (to enable parallelization features of |molcas|)
* Global Arrays (version 5 or higher, an alternative to the built-in DGA library)

.. warning::

   The DGA library is not available in |openmolcas|.

The Graphical User Interface codes in |molcas| require additional software,
including OpenGL and glut library. However, in most of the cases there is no need
to install these libraries, since executables for GUI are included into the
distribution, or they can be downloaded from |molcas| webpage (|MolcasWWW|).

In order to get TaskFarm working, the following packages should be installed
on your system:

* Debian/Ubuntu: :file:`libipc-shareable-perl`, :file:`libparallel-forkmanager-perl`
* RedHat/CentOS: :file:`perl-IPC-Shareable`, :file:`perl-Parallel-ForkManager`
* Other OS: :file:`Shareable.pm`, and :file:`ForkManager.pm` should be available from
  :variable:`$PERL5LIB`

Windows
^^^^^^^

To install |molcas| under MS Windows (98/NT/XP/7/8) one should install Cygwin
(freeware from RedHat Inc., which can be downloaded from |CygwinWWW|).
The minimal installation of Cygwin to run |molcas| includes:

* check that user name (under Windows) does not contain spaces
* select a disk, which has enough space for installation of Cygwin and |molcas|
* install Cygwin to the root of selected disk with all defaults
* run setup again and install the following packages: Devel\ :math:`\rightarrow`\gcc-fortran,
  Devel\ :math:`\rightarrow`\make, Devel\ :math:`\rightarrow`\gcc-gcc, Utils\ :math:`\rightarrow`\time, Perl\ :math:`\rightarrow`\perl
* optionally install editors: Editors\ :math:`\rightarrow`\mc, Editors\ :math:`\rightarrow`\vim
* run cygwin.bat to create Cygwin environment for the user
* copy |molcas| tar file into your home directory in Cygwin, and
  proceed with installation in the same way as under Linux.

MacOS
^^^^^

Installation of |molcas| under MacOS requires installation of the Apple
Developer Tools (Xcode) and a Fortran compiler. These programs could be
downloaded from:

| https://developer.apple.com/xcode/downloads/
| https://opensource.apple.com/
| https://gcc.gnu.org/wiki/GFortranBinaries#MacOS
| http://hpc.sourceforge.net/
| https://www.macports.org

However, if you are looking for an out of the box solution, you can download a Free PGI for Mac OS X
distribution available at https://www.pgroup.com/products/freepgi/index.htm

Preparing the installation
..........................

In order to install |molcas| you need to choose a directory
where the |molcas| driver script is to be installed. The driver
executes scripts and programs form the |molcas| package and must be
located in a directory included into the :variable:`PATH` variable.
Usually this will be :file:`/usr/bin` when installing as root,
and :file:`~/bin` when installing as an unprivileged user.

The driver script :file:`molcas` uses the value of the environment variable
:variable:`MOLCAS` to identify which version to use. The major advantage with this
mechanism is that it is easy to switch between different versions of |molcas|
by simply changing the environment variable :variable:`MOLCAS`.
However if the current directory is a subdirectory (up to 3rd level) of a
|molcas| tree, the latter will be used regardless of the value of the :variable:`MOLCAS` variable.

|molcas| itself can be located in any place on the disk.
The installation can be done by root, or by an unprivileged user.
In the later case you can copy the :file:`molcas` driver script to an appropriate
location, e.g. :file:`/usr/local/bin`, after the installation.

All files are contained in a tar archive file
with the name :file:`molcasXX.tar.gz`, where :file:`XX` depends on the version number, you need to
uncompress the file with the command
:command:`gunzip molcasXX.tar.gz`, and untar the package with
:command:`tar -xvf molcasXX.tar`.

.. _sec\:configure_molcas:

Configuring |molcas|
--------------------

Before you can build |molcas| you have to configure it.
Most common platforms have been setup by the |molcas| developers, so for a serial
installation with default settings for compiler and compiler flags configuration
of |molcas| can be done without specifying any special extra options.

There are two ways to configure |molcas|: using the :command:`configure` script (alternative 1)
or using :command:`cmake` (alternative 2). The latter is more recent and does not support all the
original platforms, but it supports most commonly used environments
and is easier as it is usually able to autodetect the necessary libraries.
The CMake alternative also makes it easier to synchronize different installations
based on the same source. Therefore, we recommend you to use alternative 2.
If you are already familiar with building |molcas| 8.0 or an earlier version,
it might be easier to keep using alternative 1, as you can then port your exisiting
configuration options.
Note that for certain external software (e.g. DMRG), CMake is mandatory,
and thus alternative 2 is needed.

For new users, use of the :command:`setup` script is recommended. It will also allow
you to choose between the two configuration alternatives if CMake is detected.

Simple configuration with the setup script
..........................................

If you are new to building |molcas|, it is recommended that you use the
:file:`setup` script to configure |molcas|. Advanced users and/or users that need
further customization should skip this section and use one of the alternatives
in the next sections to configure |molcas|.

.. compound::

  The :file:`setup` script will prompt you interactively about the most important settings.
  To get started, first unpack the |molcas| source and then run::

    ./setup

  in the main |molcas| directory. Answer the questions and then proceed to
  :numref:`sec:building_molcas` to build |molcas|.

For advanced users that need further customization, one of the alternatives
in the next sections will be needed.

Advanced configuration with the configure script (alternative 1)
................................................................

You can start the configuration by running the :file:`configure` script::

  ./configure

To know which flags can be used, run :command:`{./configure -h`. The above command
(without any flags) will use a default configuration, i.e. a serial |molcas|
installation using built-in BLAS/LAPACK.

When configuration is finished, you should review the log file
:file:`configure.log`
to see if everything is ok.
There is no harm in running the configuration script even if it should fail,
you simply rerun it with correct parameters.

If the configuration step was not successful, you probably are missing some
prerequisite software, or this software is located in an unusual location on the disk.
In the later case you might need to update your :variable:`PATH`, or use flag :command:`-path`
in :command:`configure`.

|molcas| supports out-of-source installation. If for some reason,
you would like to install molcas under a separate tree, you can create
a directory, and call :command:`configure` with appropriate flags, e.g. ::

  mkdir $HOME/molcas
  cd $HOME/molcas
  /sw/molcas_dist/configure -speed safe

The list all the options for :command:`configure`, run ::

  ./configure -help

Advanced configuration with CMake (alternative 2)
.................................................

Start configuration by creating a build directory and then running the :command:`cmake`
program with the location of the |molcas| source. You can create the build directory
as a subdirectory of the |molcas| root directory, but we recommend you create it
outside. For example, suppose the |molcas| root is located at :file:`molcas` in your
home directory and you want to build it at :file:`molcas-build`::

  mkdir ~/molcas-build
  cd ~/molcas-build/
  cmake ~/molcas/

After the first run, CMake keeps a cache of the configuration around
as a file :file:`CMakeCache.txt`.
If you want to change certain options (e.g. use a different compiler),
you need to remove this file first. As an example, if we wish to reconfigure
with the intel compilers instead we do::

  rm CMakeCache.txt
  CC=icc FC=ifort cmake ~/molcas/

Once the cache file is there, subsequent additional options do not
require pointing to the |molcas| source, just to the build directory
itself. Here we add an option to build parallel |molcas|::

  cmake -DMPI=ON .

When using MPI, CMake will pick up the correct compiler and MPI library from the wrappers.
So if you configure with MPI or you later wish to use a different MPI library/compiler wrapper, it is better
to remove the cache file first, then (re)configure pointing to the wrappers::

  rm CMakeCache.txt
  CC=mpicc FC=mpifort cmake -DMPI=ON ~/molcas/

Summary of main options for CMake
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To see all the CMake options, you can run the :command:`ccmake` command to interactively
select the options and their values. The following is a list of some of the most
important options.

* ``BUILD_SHARED_LIBS`` --- (``ON``/``OFF``): Enable building shared libraries (reduced disk space).
* ``CMAKE_BUILD_TYPE`` --- (``Debug``/``Garble``/``RelWithDebInfo``/``Release``/``Fast``): Normally
  use ``Release``. ``RelWithDebInfo`` may be useful for reporting a problem. ``Fast`` in unsupported
  and may give wrong results in some cases.
* ``CMAKE_INSTALL_PREFIX`` --- Specify the directory where |molcas| will be installed when running
  :command:`make install`.
* ``GA`` --- (``ON``/``OFF``): Enable interface with Global Arrays (see :numref:`sec:parallel_installation`).
* ``HDF5`` --- (``ON``/``OFF``): Enable HDF5 files (portable binary format) in some programs.
* ``LINALG`` --- (``Internal``/``Runtime``/``MKL``/``ACML``/``OpenBLAS``): Select the linear algebra library (BLAS + LAPACK)
  against which to link |molcas|. ``Internal`` uses the default (and slow) Netlib version included with |molcas|. ``Runtime``
  (experimental) offers the possibility of choosing the library at run time.
* ``MPI`` --- (``ON``/``OFF``): Enable multi-process parallelization.
* ``OPENMP`` --- (``ON``/``OFF``): Enable multi-thread parallelization (usually restricted to using
  a multi-threaded version of linear algebra libraries.
* ``TOOLS`` --- (``ON``/``OFF``): Compile the tools that have CMake support.

Example 1: GCC C/Fortran compilers with GA/OpenMPI and OpenBLAS
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* 64-bit Linux OS
* MPI preinstalled (install OpenMPI or MPICH with package manager)

::

  # OpenBLAS
  tar zxvf OpenBLAS-v0.2.15.tar.gz
  cd OpenBLAS-0.2.15/
  make USE_OPENMP=1 NO_LAPACK=0 INTERFACE64=1 BINARY=64 DYNAMIC_ARCH=1 libs netlib shared
  [sudo] make PREFIX=/opt/openblas-lapack-ilp64 install
  # GA
  tar zxvf /path/to/ga-5-4b.tgz
  cd ga-5-4b/
  ./configure --enable-i8 --with-blas8 --with-lapack8 --with-scalapack8 --prefix=/opt/ga54b-ilp64.OpenMPI
  make
  [sudo] make install
  # Molcas
  tar zxvf molcas.tgz
  cd molcas
  mkdir build && cd build/
  export GA=/opt/ga54b-ilp64.OpenMPI
  export OPENBLASROOT=/opt/openblas-lapack-ilp64
  CC=mpicc FC=mpifort cmake -DMPI=ON -DGA=ON -DLINALG=OpenBLAS ../
  make -j4

Any other configurations can be easily derived from the above by simply
removing the unneeded optional stuff.

Example 2: Intel C/Fortran compilers with GA/IntelMPI and MKL
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* 64-bit Linux OS
* Intel Compilers/MPI preinstalled (usually on a cluster as a module)
* Infiniband (OpenIB)

::

  # make sure Intel compilers/MPI are loaded
  module load ...
  # check compilers/environment
  type mpiicc
  type mpiifort
  echo $MKLROOT
  # GA
  tar zxvf /path/to/ga-5-4b.tgz
  cd ga-5-4b/
  ./configure --prefix=/opt/ga54b-ilp64.IntelMPI --enable-i8 --with-openib \
              --with-blas8="-L$MKLROOT/lib/intel64 -lmkl_intel_ilp64
                            -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm" \
              --with-scalapack8="-L$MKLROOT/lib/intel64 -lmkl_scalapack_ilp64
                                 -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core
                                 -lmkl_blacs_intelmpi_ilp64 -liomp5 -lpthread -lm"
  make
  [sudo]make install
  # Molcas
  tar zxvf molcas.tgz
  cd molcas
  mkdir build && cd build/
  export GA=/opt/ga54b-ilp64.IntelMPI
  CC=mpiicc FC=mpiifort cmake -DMPI=ON -DGA=ON -DLINALG=MKL ../
  make -j4

.. _sec\:building_molcas:

Building |molcas|
-----------------

When the configuration step (:numref:`sec:configure_molcas`) is completed
successfully, you can build |molcas|.
This is simply done by typing :command:`make` in the |molcas| root directory.
It is recommended that you save the output from :command:`make` in a log file
for tracing of potential problems. ::

  make > make.log 2>&1

.. compound::

  In order to speed up the build process, you can perform a parallel compilation.
  When you configure |molcas| with CMake, a simple ::

    make -jN

  will suffice, while if you configured |molcas| with the :command:`configure` script,
  you need to build in two steps::

    make -jN build
    make install

  In the above commands, ``N`` should be replaced with the number of threads
  you wish to use for building. Typically, a number equal to the number of cores
  will be fine.

When |molcas| is being compiled some compilers give a lot of warnings.
These are not serious in most cases. We are working on eliminating them, but
the job is not yet completely finished.

After |molcas| has been built correctly, you should absolutely run a basic
verification to ensure that the installation is sane. See the :ref:`next section <sec:verify>`
for details on verification.

.. _sec\:verify:

Verifying the |molcas| installation
...................................

After a successful build of |molcas| you should verify that the various
modules run correctly. Directory :file:`test/`
contains various subdirectories with test inputs for |molcas|. Use the command
:command:`molcas verify [parameters]` to start verification.
Running this command without parameters will check
the main modules and features of |molcas| and we recommend this default
for verifying the installation.
You can also specify a keyword as argument that translates into a
sequence of test jobs,
or you can specify a list of test jobs yourself.
The script has extended documentation, just run
:command:`molcas verify --help`.

To generate a report after performance tests you should execute
a command :command:`molcas timing`. The report is then located in the file
:file:`test/timing/user.timing`. The results of benchmark tests
for some machines are collected on our webpage
http://www.molcas.org/benchmark.html
At the completion of the test suite a log of the results is generated
in the file :file:`test/results`. If installation was performed
by another user (e.g. root), you can redefine the location of
output files by adding the flag :command:`-path PATH`.
Each test job is signaled as either
ok of failed. If there are any failed jobs, the outputs are saved in
:file:`Test/Failed_Tests`. Each test job tests for a resulting checksum
for the modules tested. This checksum is typically the energy for a
wavefunction program such as :program:`RASSCF`, whereas other types of
codes use other checksums.

The checksums will not match exactly with our reference values since
different machines use different arithmetics. We have tried to make
the acceptable tolerances as small as possible and at the same time
make all tests pass successfully. It might be the case that your
particular platform will produce one or more results that are
just outside our tolerances, and in such a case the test is most
likely ok.

.. _sec\:building_html_doc:

GUI and documentation
---------------------

Normally, there is no need to build the GUI used in |molcas|, since we
provide executables for most common platforms.
You can download executables for the GUI from the |molcas| webpage (|MolcasWWW|).

In order to build documentation in various formats, use the command :command:`make doc`.
Alternatively the documentation is also available no the web page (|MolcasWWW|).

