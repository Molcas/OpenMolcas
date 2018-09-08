Maintaining the package
=======================

.. only:: html

  .. contents::
     :local:
     :backlinks: none

Tailoring
---------

|molcas|, as shipped, is configured with some default settings. You can change
some of these easily.
You can change default settings used in |molcas| (like memory usage, default
scratch area, policy in saving files, etc.)
by editing |molcas| resource file: global resource file :file:`$MOLCAS/molcasrc` or
user resource file :file:`$HOME/.Molcas/molcasrc`.

.. _sec\:dynamic_memory:

Dynamic memory
..............

Most modules in |molcas| utilize dynamic memory allocation. The amount of
memory each module allocate is controlled by the environment variable
:variable:`MOLCAS_MEM`. The amount of memory allocated is

* :variable:`MOLCAS_MEM` is undefined --- 1024 MB of memory is allocated
  (on 32 bit installation)
* :variable:`MOLCAS_MEM`\=\ ``nn`` --- ``nn`` MB is allocated.
  If this amount cannot be allocated, the module stops.

.. _sec\:disk_usage:

Disk usage
..........

Today many workstations utilize 64-bit integers
and addressing. However, old UNIX workstations and PC's had 32-bit integers
resulting in a file size limit of 2 GB.
To circumvent these limitations, the I/O routines
of |molcas| support multifile files, where a "file" is in reality a
logical file consisting of several physical files. The size limit of
these physical files is controlled by the environment variable
:variable:`MOLCAS_DISK` according to

* :variable:`MOLCAS_DISK` is undefined --- The modules will use a 2 GB size of the
  physical files. This might be the appropriate setting for machines
  with 32-bit addressing.
* :variable:`MOLCAS_DISK`\=\ ``nn`` --- The modules will use a ``nn`` MB size
  of the physical files.

To use files with a size bigger than 2 GB |molcas| should be compiled as 64-bit
executable.

Improving CPU performance
.........................

|molcas| is shipped with a number of default setup files located
in directory :file:`cfg/`. The defaults in these files are set to
a fairly safe level, but not necessary optimal. What you can change
to improve performance is

* Compiler flags
* Mathematical (blas) libraries

The simplest way to set up optimization level, and/or compile |molcas|
with various BLAS libraries is to use :command:`configure -setup`. This
interactive script helps to make a proper selection of flags for
improvement of |molcas| performance.

If you do decide to try to improve the performance we recommend that you
create a new setup file, for example, :file:`cfg/local.cfg` and
modify this file. It is not unlikely that your attempts to optimize
the codes will lead you to a case where some modules work and others do not.
In such a scenario it can be fruitful to have two copies of
|molcas|, one "safe" where all modules work and one "fast" where
some modules do not function properly.

Changing the compiler flags is the easiest. Using the most
aggressive optimization flags do sometimes lead to problems for
some of the modules. We have tried to choose an optimization level
that yields functioning code, but still reasonable fast.
For some systems there is a predefined set of compiler flags for
aggressive optimization. To compile |molcas| with these flags you
should run :command:`configure` with flag :command:`-speed fast`.
Note that this aggressive optimization level is not supported
by the |molcas| team. In other words, you are using it at your own
risk.

For some platforms you can utilize the vendor blas libraries. This
will certainly yield better performance, but may not work on all
platforms.

.. compound::

  During configuration of |molcas| it is possible to specify
  an external BLAS/LAPACK library. Use a flag :command:`-blas TYPE`
  to specify the type of BLAS libary: lapack (for a standard lapack
  library), Goto (for GotoBLAS), Atlas (for ATLAS), MKL (for Intel MKL).
  You should also specify a flag ::

    -blas_lib -Wl,--start-group -L/path/to/blas -lmy-blas -Wl,--end-group

  specifying the link options.
  For example, to configure |molcas| with Intel MKL library,
  you should issue a command ::

    ./configure -compiler intel -blas MKL -blas_lib -Wl,--start-group /opt/intel/mkl/lib/intel64 -lmkl_gf_ilp64 -lmkl_sequential -lmkl_core -Wl,--end-group

.. compound::

  To compile |molcas| with CUDA BLAS library, first, you have to compile
  the fortran wrapper provided by nVIDIA: ::

    CUDA=/path/to/cuda/
    FLAGS=-m64
    gcc $FLAGS -I$CUDA/include/ -I$CUDA/src/ -c $CUDA/src/fortran_thunking.c -o \
    $MOLCAS/lib/fortran.o
    ./configure -blas CUDA -blas_dir $CUDA/lib

  or, if on a 64bit system: ::

    ./configure -blas CUDA -blas_dir $CUDA/lib64

.. It is also possible to make a manual installation of a vendor
   supplied BLAS library.
   One should issue commands :command:`molcas uninstall blas_util`,
   :command:`molcas uninstall essl_util` and :command:`molcas uninstall %lapack_util` to remove BLAS/LAPACK related directories
   from the |molcas| source code, then export XLIB variable to set the
   location of blas library, e.g. :command:`XLIB="-lblas"; export XLIB`,
   and finally reconfigure and build |molcas|. If the library is in a
   nonstandard location you may have to issue the command
   :command:`XLIB="-Lpath_to_lib -lblas"; export XLIB`. Alternatively,
   define XLIB in the system specific configuration file.

After making changes to the setup files you have to issue the commands
:command:`make veryclean`, :command:`./configure` and :command:`make` in the |molcas| root
directory. It is highly recommended to run the verification suite after
any changes in configuration file.

.. MT:sec:prgm:

Customizing handling of files
.............................

The location and attributes of files used by |molcas| are defined in PRGM files.
The master copy of these files is located at the :file:`data` directory in the |molcas| root.
A user can copy these files and modify them. The highest priority is given
to the files located in the subdirectory :file:`prgm` in the current (:variable:`CurrDir`) directory,
next the :file:`$HOME/.Molcas` and finally the original location at :file:`$MOLCAS/data`.

The simplest way to manipulate prgm files is to use the :command:`molcas prgm` command.
A command :command:`molcas prgm init global` makes a copy of the prgm files in the :file:`.Molcas`
directory. :command:`molcas prgm init local` creates a :file:`prgm` subdirectory, and copies the prgm files
into it.
The editing of the files can be performed by your favourite editor, or by the :command:`molcas prgm` script.

The structure of PRGM files is simple. The field ``(file)`` is followed by the
FORTRAN name (as this file is known for |molcas|), real file name (as this file is known by the Operating system), and finally the attributes.

The attributes are listed here:

====== =====================================================
``ro`` the file is accessed for reading
``rw`` the file is accessed for reading and writing
``s``  the file will be saved (copied) from the scratch area
``m``  the file will be saved (moved) from the scratch area
``g``  the file can be visualized by program "gv/luscus"
``t``  the file is an ASCII text
``*``  multifile
``.``  multifile (for internal use)
``l``  use lustre filesystem in parallel run
``p``  the file will be deleted
``e``  the file will be placed to memory (see FiM)
``f``  use an alternative file location (FastDir)
====== =====================================================

A command :command:`molcas prgm +x ScfOrb` will add the attribute ``x`` to all ScfOrb files. A regexp can be used for a filename.
An opposite command :command:`molcas prgm -x ScfOrb` will remove the attribute ``x``.
A command :command:`molcas prgm comp [global]` shows the list of modifications in the local (:file:`prgm`), or global (:file:`.Molcas`)
directories.

.. _MT\:sec:fim:

Improving I/O performance
.........................

In order to activate this technology for a |molcas| scratch file, one needs to
do three things. First, please edit an external resource :file:`*.prgm` (for example,
:file:`$MOLCAS/data/seward.prgm`) from the :file:`$MOLCAS/data/` directory. If you don't
have access to the root |molcas| directory, then you can simply copy the
needed resource file into your home :file:`$HOME/.Molcas/` directory and edit it there.
The editing of the file consists in adding the "``e``" character to its
attributes: ::

  original: (file) ORDINT "$WorkDir/$Project."OrdInt rw*
  modified: (file) ORDINT "$WorkDir/$Project."OrdInt rw*e

Second, you need to set up the :variable:`MOLCAS_FIM` environment variable to ``1``
i.e.: ::

  export MOLCAS_FIM=1

The third and the final step is to specify the :variable:`MOLCAS_MAXMEM` (:math:`\geq`\ :variable:`MOLCAS_MEM`) parameter such that the
:variable:`MOLCAS_MAXMEM`\ |-|\ :variable:`MOLCAS_MEM` difference (in MW) is sufficient to host an entire
file in RAM. In other words, the :variable:`MOLCAS_MAXMEM`\ |-|\ :variable:`MOLCAS_MEM` difference should
exceed the original filesize.

In general, not all |molcas| files are suitable for placing in RAM. In
particular, it is a bad idea to activate FiM for :file:`RUNFILE`. In order to
identify which |molcas|'s files are proper candidates for FiM, you can simply
inspect the section ``II. I/O Access Patterns`` from a |molcas|'s output.
All files with high ratio of I/O ``random Write/Read calls`` are good candidates for
FiM. In particular case of the :program:`SEWARD` module, the :file:`ORDINT` file is
a very good candidate for FiM: ::

  II. I/O Access Patterns
  - - - - - - - - - - - - - - - - - - - -
  Unit  Name               % of random
                         Write/Read calls
  - - - - - - - - - - - - - - - - - - - -
   1  RUNFILE             28.6/  11.5
   2  ORDINT             100.0/  24.0
   3  DNSMAT               0.0/   0.0
   4  TWOHAM               0.0/   0.0
   5  GRADIENT            88.9/   0.0
   6  DNSMAX               0.0/   0.0
   7  TWOHAX               0.0/   0.0
   8  SODGRAD             85.7/   0.0
   9  SOXVEC              85.7/   0.0
  10  SODELTA             88.9/   0.0
  11  SOYVEC              88.9/   0.0
  12  ONEINT             100.0/  53.3
  - - - - - - - - - - - - - - - - - - - -

Applying patches
----------------

All program systems do contain bugs and |molcas| is certainly no exception.
We prepare patches for all problems as soon as we identify and fix the
problem.

.. You can get these patches from our web server in an easy and automatic way.

For important updates we provide Service Packs. A service pack is a shell
script, which makes a backup of your current |molcas| installation, and installs
updates.

Local modifications
...................

|molcas| is shipped with source code so you can make modifications
yourself. You are, of course, responsible for the correctness of any
such modification.

If you do make changes/additions to the source code that you feel is of
interest to other users, we encourage you to make these available.
Perhaps the best mechanism is to use the bulletin board on out homepage:
|MolcasWWW|.

Check Molcas Programming Guide for a detailed description of development
and distribution of modified code in |molcas|.
