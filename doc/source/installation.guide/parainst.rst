.. _sec\:parallel_installation:

Parallel Installation
=====================

.. only:: html

  .. contents::
     :local:
     :backlinks: none

Installation of |molcas| for execution in multi-processor environments can be a
bit more involved than the standard installation, this chapter considers those
particulars not covered previously.

The parallelization of |molcas| is achieved through the use of the Global Arrays
(GA) API and some direct MPI calls. The API can be linked to an external GA library
or to our own DGA library (an internal PGAS framework built upon MPI-2).

.. warning::

   The DGA library is not available in OpenMolcas.

When using DGA (the default), the current list of supported MPI-2.2 implementations is given below:
MPICH2/MPICH3, MVAPICH2, OpenMPI, Intel MPI.

When one wants to use an external GA library, it has to be configured and
compiled separately. In that case, please read the section on using an
external GA installation to properly configure and install GA first!!!

**IMPORTANT**: not all modules support distribution of work and/or
resources through parallel execution, and even if they do it might be that some
functionality is limited to serial performance. This is a list of core modules
which *can* benefit from parallel execution: gateway, seward, scf, rasscf,
caspt2. More detailed information regarding parallel behaviour can be found in
the documentation of the respective module and in the table at the beginning of
the manual about supported parallelism. If no information is available, you
should conclude that there is nothing to be gained from parallel execution.

Supported MPI implementations
-----------------------------

Most probably, you will use a free MPI-2 implementation such as MPICH2/MPICH3,
MVAPICH2, or Open MPI.

* MPICH2: :file:`https://www.mpich.org/`
* MPICH3: :file:`https://www.mpich.org/`
* MVAPICH2: :file:`http://mvapich.cse.ohio-state.edu/`
* Open MPI: :file:`https://www.open-mpi.org/`

**NOTE**: Open MPI versions older than v1.6.5 are not supported. More specifically,
only Open MPI v1.6.5, and v1.8.1 are tested and known to work correctly with |molcas|.

It is a very good idea to verify that the correct compiler environment is
present before configuring |molcas|. You should therefore check that the
backend compiler of the wrappers is correct by running :command:`/path/to/mpif77
-show` (MPICH2/MPICH3 and MVAPICH2) or :command:`/path/to/mpif77 --showme` (Open MPI),
which will list the actual executed command. If the backend compiler seems to
be correct, also try to run it to see if it is properly detected (on some
clusters you will need to load the appropriate module for the compiler). If all
is well, you should be able to configure |molcas| without any problems.

It is highly recommended to use the compiler that was used for the MPI library
to build GA (optional) and |molcas| to avoid compatibility issues. However, if you really
want to use a different compiler than the compiler which was used for building
the MPI library, you can do so by passing the :command:`-fc` and :command:`-cc`
command line arguments (MPICH2/MPICH3 and MVAPICH2) to the wrappers, or setting the
environment variables :variable:`OMPI_F77`/:variable:`OMPI_F90` and :variable:`OMPI_CC` (Open MPI).

Several commercial MPI implementations exist such as HP-MPI, IBM's MPI-F, Intel
MPI, SGI's MPT. Currently we only support Intel MPI. For the others that
are not (yet) supported, it is recommended to either configure |molcas| without parallel
options and change the :file:`Symbols` file after the serial configuration, or rely on
cmake to use the correct options.

Please refer to the documentation of your MPI implementation for details on how
to build programs, i.e. which wrappers to use and if necessary what libraries
you need to link in.

Using an external Global Arrays installation (optional step)
------------------------------------------------------------

If you wish to use an external GA library, it has to be installed before
you build |molcas|. You could e.g. use this if you have trouble with the built-in
DGA solution.
The installation instructions may be found at the Global Arrays home page:
:file:`http://hpc.pnl.gov/globalarrays/`

Note that any problems with installation or other issues specific to GA are
best resolved by contacting the GA authors directly, rather than the
|molcas| group. It is therefore a very good idea to run the GA
testing code as a job on the cluster where you want to use |molcas| to make sure
that it works properly before continuing to install |molcas|.

Global Arrays needs to be installed with 8-byte integer support using
the flag(s) :command:`--enable-i8 --with-blas8[=...] [--with-scalapack8[=...]]`,
and for infiniband clusters you probably need to use the :command:`--with-openib` flag.
When linking to an external library, e.g. the Intel MKL, do not forget to include
the proper :command:`ilp64` library versions.

Please read the documentation of GA for more details about installation.

General overview of the procedure with the configure script (alternative 1)
---------------------------------------------------------------------------

In the simplest case, the parallel version of |molcas| may be installed
simply by specifying the flag :command:`-parallel`
to :file:`configure`. For example: ::

  ./configure -parallel

When using an external GA, pass the location of the installation to |molcas| configure: ::

  ./configure -parallel -ga /opt/ga-5.1

When the locations of the MPI :file:`lib` and :file:`include` directories is set
incorrectly, you can specify them by setting their common root directory with
the :command:`par_root` flag or if they are in different directories you can use the
separate :command:`par_inc` and :command:`par_lib` flags: ::

  ./configure -parallel -par_root /usr/lib/openmpi
  ./configure -parallel -par_inc /usr/lib/openmpi/include -par_lib /usr/lib/openmpi/lib

More likely, some individual tailoring will be required, the following
summarizes the necessary steps:

#. Check that the correct wrapper compilers were detected, as specified in :file:`$MOLCAS/Symbols`.
#. If needed, change the :variable:`F77`/:variable:`F90` and :variable:`CC` variables in the :file:`Symbols` file for any custom modifications you made to the wrappers.
#. Optionally install (and test) the external Global Arrays library.
#. Check the command for executing binaries in parallel, as specified by :variable:`RUNBINARY` in
   :file:`$MOLCAS/molcas.rte`.
#. Install (and test) |molcas|.

Provided that steps 1--4 can be successfully accomplished, the installation
of |molcas| itself is unlikely to present many difficulties.

General overview of the procedure with cmake (alternative 2)
------------------------------------------------------------

CMake accepts two main flags for parallel installation, one to specify the
use of parallelization :command:`-DMPI=ON`, and one to speficy an external GA
library :command:`-DGA=ON` instead of DGA (the default is :command:`-DGA=OFF`, meaning no
external GA is used, so do not confuse the option ``-DGA`` which means "define GA"
with DGA). When using the latter :command:`-DGA=ON` flag, make sure the :command:`GAROOT` environment
variable is exported and contains the path of the GA installation, before running
cmake.

CMake will determine an appropriate MPI libary based on the compiler it finds, so in order
to use a specific MPI library, just make sure the :variable:`CC` and :variable:`FC`
variables point to the correct MPI wrappers!

The whole procedure is summarized below (square brackets showing optional commands): ::

  [export GAROOT=/path/to/external/GA]
  [CC=/path/to/mpicc] [FC=/path/to/mpifort] cmake -DMPI=ON [-DGA=ON] /path/to/molcas
  make [-j4]

Running |molcas| in parallel
----------------------------

A few comments on running on a cluster:

The very old MPICH versions sometimes need a file with a list of the nodes the job at hand is allowed
to use. At default the file is static and located in the MPICH installation
tree. This will not work on a workstation cluster, though, because then all
jobs would use the same nodes.

Instead the queue system sets up a temporary file, which contains a list of the
nodes to be used for the current task. You have to make sure that this filename
is transfered to $mpirun. This is done with the :command:`-machinefile` flag. On a
Beowulf cluster using PBS as queue system the :variable:`RUNBINARY` variable in
:file:`$MOLCAS/molcas.rte` should look something like: ::

  RUNBINARY='/path/to/mpirun -machinefile $PBS_NODEFILE -np $MOLCAS_NPROCS $program'

The newer MPICH2/MPICH3 as well as MVAPICH2, which works through the use of the HYDRA daemons and does not need
this command line argument, as well as Open MPI most likely only need the :command:`-np
$MOLCAS_NPROCS` command line option. They use mpiexec instead of mpirun.

.. compound::

  Parallel execution of |molcas| is achieved by exporting the environment
  variable :command:`MOLCAS_NPROCS`, for example when running on 4 nodes use: ::

    export MOLCAS_NPROCS=4

  and continuing as usual.

In this section, we assume you will be using PBS on a cluster in order to
submit jobs. If you don't use PBS, please ask your system administrator or
consult the cluster documentation for equivalent functionality.

Example of a submit script
..........................

::

  #!/bin/sh
  #PBS -l walltime=10:00:00
  #PBS -l nodes=4
  #PBS -l pmem=3000mb

  ######## Job settings ###########
  export MOLCAS_MEM=800
  export SUBMIT=/home/molcasuser/project/test/
  export Project=test000
  export MOLCAS_NPROCS=4

  ######## modules ###########
  . use_modules
  module load intel/11.1
  module load openmpi/1.4.1/intel/11.1

  ######## molcas settings ###########
  export MOLCAS=/usr/local/molcas80.par/
  export WorkDir=/disk/local/

  ######## run ###########
  cd $SUBMIT
  molcas $Project.input -f

Memory
......

The maximum available memory is set using the PBS option pmem. Typically,
:variable:`MOLCAS_MEM` will then be set to around 75% of the available physical
memory. So for a parallel run, just divide the total physical memory by the
number of processes you will use and take a bit less. For example, for a system
with 2 sockets per node and 64 GB of memory, running 1 process per socket, we
would set pmem to 30000 MB.

I/O
...

The important thing to consider for I/O is to have enough scratch space
available and enough bandwidth to the scratch space. If local disk is large
enough, this is usually preferred over network-attached storage. |molcas|
requires the absolute pathname of the scratch directory to be the same across
nodes.

Pinning
.......

Process pinning is sometimes required to achieve maximum performance. For CASPT2
for example, processes need to be pinned to their socket or NUMA domain.

The pinning configuration can usually be given as an option to the MPI runtime.
With Intel MPI for example, one would set the :variable:`I_MPI_PIN_DOMAIN`
variable to :command:`socket`. Alternatively, you can use a third-party program
to intervene on your behalf, e.g. https://code.google.com/p/likwid/.
Please ask your system administrator how to correctly pin your processes.

GA specific issues
..................

When using GA, several problems can occur when trying to run jobs with a large
amount of memory per process. A few example error messages are given here with
their proposed solution.

::

  (rank:0 hostname:node1011 pid:65317):ARMCI DASSERT fail.
   src/devices/openib/openib.c:armci_pin_contig_hndl():1142
   cond:(memhdl->memhndl!=((void *)0))

The error output in the Molcas errfile (stderr) then says: ::

  Last System Error Message from Task 2:: Cannot allocate memory

Related messages that display a problem with :file:`armci_server_register_region`
instead of :file:`armci_pin_contig_hndl` can also occur, and point to similar problems.

This can have two causes:

* Some parameters of the Mellanox :file:`mlx4_core` kernel module were
  set too low, i.e., :file:`log_num_mtt` and :file:`log_mtts_per_seg`.
  These should be set according to the instructions on
  :file:`https://community.mellanox.com/docs/DOC-1120`. Values of 25 and 0
  respectively, or 24 and 1 should be fine.
* The "max locked memory" process limit was set too low. You can check
  this value by running :command:`ulimit -a` or :command:`ulimit -l`. Make
  sure you check this through an actual job! Easiest is to start an
  interactive job and then execute the command. The value should be set
  to unlimited, or at least to the amount of physical memory available.

::

  0: error ival=4 (rank:0 hostname:node1011 pid:19142):ARMCI DASSERT fail.
   src/devices/openib/openib.c:armci_call_data_server():2193
   cond:(pdscr->status==IBV_WC_SUCCESS)

This error is related to the value of the variable
:variable:`ARMCI_DEFAULT_SHMMAX`, try setting it at least to 2048. If this is
still too low, you should consider patching GA to allow higher values.
