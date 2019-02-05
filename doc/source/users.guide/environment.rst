Commands and environment variables
==================================

This section will describe the usage of |molcas| in an
UNIX environment.

Production jobs using |molcas| in an UNIX environment can be
performed as batch jobs.
This requires the creation of a shell script that
contains a few simple commands. Further you need to create input for
each program module that you intend to use. This section describes the
necessary steps you have to take in order to make a successful job using
|molcas|.
Input examples for a typical |molcas|
run can be found in :file:`doc/samples/problem_based_tutorials/` directory.
Also you can use some input examples
in :file:`Test/input` subdirectory.

Commands
--------

There is a command supplied with the |molcas| package, named
:command:`molcas`, that the user issue to perform a given task.
A sequence of such commands will perform the calculation requested by
the user.

.. class:: commandlist

:command:`molcas`
  This command tells which |molcas| installation will be used, and gives
  some help about usage of |molcas| command

:command:`molcas` :file:`input-file`
  This command executes a command in the |molcas| system.

:command:`molcas help` :program:`prgm`
  This command gives the list of available keywords for program :program:`prgm`.

:command:`molcas help` :program:`prgm keyword`
  This command gives description of a :program:`keyword`.

:command:`molcas help` :program:`environment`
  This command gives a list of |molcas| specific environment variables.

:command:`molcas help` :program:`basis element`
  This command gives a list of basis sets available for an :program:`element`.

The following is an example of running |molcas| by using a single input file: ::

  molcas $Project.input

An alternative way of running |molcas| as a sequence of separate calls: ::

  molcas $Project.seward.input    # Execute seward
  molcas $Project.scf.input       # Execute scf

By default, the output will go directly to the screen. It can be redirected
by using flag :command:`-f`, e.g. :command:`molcas -f water.inp` will store the output
in :file:`water.log` and :file:`water.err` files.

The default behavior of |molcas| execution can be altered by setting environment variables.

Project name and working directory
----------------------------------

When running a project, |molcas| uses the variable
:variable:`Project` giving a project name, and a scratch directory defined by
the variable :variable:`WorkDir`.
This serves the purpose of maintaining structure of the files and
facilitating automatic file mapping.

There are several ways to set up these variables.
By default, the name of the Project constructed from the name of the input file,
by removing the last suffix, e.g. for example for an input name :file:`Water.SCF.input`
the :variable:`Project` name will be :variable:`Water.SCF`.
Alternatively,
user can set environment variable :variable:`Project`, or :variable:`MOLCAS_PROJECT`.

Scratch directory can be set by environment variable :variable:`MOLCAS_WORKDIR`.
If it is set to value ``PWD``, current directory will be used. Otherwise,
it can be set to a directory name. In this case scratch area will be located
in a subdirectory :file:`$MOLCAS_WORKDIR/$Project`. It is also possible to
overwrite the value of scratch area, by setting environment variable
:variable:`WorkDir`.

* :command:`Project=...; export Project`
* :command:`WorkDir=...; export WorkDir`

|molcas| modules communicates between each other via files, located in the :variable:`WorkDir`.
The description of internal filenames and file mapping can be found at Appendix.

Input
-----

.. compound::

  When you have decided which program modules you need to use to perform your
  calculation, you need to construct input for each of these. There is no
  particular structure enforced on the input files, but it is recommended that
  you follow:

  * :file:`$Project.<prgm-name>.input`

  which is the name of the input files assumed in the sample shell script.

Preparing a job
---------------

When you prepare a job for batch processing, you have to create a shell script.
It is recommended that you use the sample shell script supplied with
|molcas| as a starting point when building your own shell script.
The following steps are taken in the shell script:

#. Define and export the |molcas| variables

   * Project (or use :variable:`MOLCAS_PROJECT`)
   * WorkDir (or :variable:`MOLCAS_WORKDIR`)

#. Issue a sequence of |molcas| commands.
#. Remove the scratch directory and all files in it.

The following is an example of a shell script. ::

  Project=HF; export Project                               # Define the project id
  WorkDir=/temp/$LOGNAME/$Project.$RANDOM; export WorkDir  # Define scratch directory
  molcas $Project.input                                    # Run molcas with input file, which
                                                           # contains inputs for several modules
  rm -r $WorkDir                                           # Clean up

The file :file:`$ThisDir/$Project.input` contains the ordered sequence
of |molcas| inputs and the EMIL interpreter will call the appropriate
programs. See section :ref:`UG:sec:EMIL` for an explanation of the
additional tools available in the EMIL interpreter.

The following is an example of a shell script to be submitted for batch
execution. ::

  Project=HF; export Project                               # Define the project id
  WorkDir=/temp/$LOGNAME/$Project.$RANDOM; export WorkDir  # Define scratch directory
  molcas $Project.seward.input                             # Execute seward
  molcas $Project.scf.input                                # Execute scf
  rm -r $WorkDir                                           # Clean up

An alternative way to control the usage of the WorkDir is to use flags in molcas command:

.. class:: commandlist

:command:`-new`
  clean WorkDir before the usage

:command:`-clean`
  clean WorkDir after the usage

Note, that if you configured your working environment by using :command:`setuprc` script,
the only command you have to place into the shell script is: ::

  molcas $Project.input

.. _UG\:sec\:sysvar:

System variables
----------------

|molcas| contains a set of system variables that the user can
set to modify the default behaviour of |molcas|. Two of them
(Project and WorkDir) must be set in order to make |molcas| work at all.
There are defaults for these but you are advised not to use the defaults.

There are several ways of using |molcas| environment variables:

* These variables can be exported in your shell script ::

    export MOLCAS_MEM=512
    molcas input

* These variables can be included into |molcas| input: ::

    * begin of the input file
    >>> export MOLCAS_MEM=512

      . . .

* variables can be included directly into :command:`molcas` command in the form: ::

    molcas MOLCAS_MEM=512 input

The simplest way to set up default environment for |molcas| is
to use script :file:`setuprc`, which can be run as command
:command:`molcas setuprc`. This interactive script creates
a resource file :file:`molcasrc`, located either in :file:`$MOLCAS` or :file:`$HOME/.Molcas`
directory. The priority of these settings is: user defined settings
(e.g. in :command:`molcas` command), user resource file, |molcas| resource file.

Two flags in |molcas| command are related to resource files:

.. class:: variablelist

:variable:`-env`
  Display current |molcas| environment
  e.g. :command:`molcas -env input` will print information about environment
  variables, used during execution of the input file.

:variable:`-ign`
  Ignore resource files
  e.g. :command:`molcas -ign input` will process input file without settings,
  which are stored in :file:`$MOLCAS/molcasrc` and in :file:`$HOME/molcasrc` files.

The most important environment variables, used in |molcas|:

.. xmldoc:: <MODULE NAME="ENVIRONMENT" LEVEL="HIDDEN">
            %%Description:
            List of environment variables

.. class:: variablelist

:variable:`Project`
  This variable can be set in order to overwrite the default name of the
  project you are running. The default (and recommended) value of the project name is the
  name of the input file (without the file extension).

  .. xmldoc:: <KEYWORD MODULE="ENVIRONMENT" NAME="Project" APPEAR="Project" KIND="STRING" LEVEL="ADVANCED">
              %%Keyword: Project <advanced>
              <HELP>
              This variable can be set in order to overwrite the default name of the
              project you are running. The default (and recommended) value of the project name
              is the name of the input file (without the file extension).
              </HELP>
              </KEYWORD>

:variable:`WorkDir`
  This variable can be used to specify directly the directory where all files
  that |molcas| creates are placed. See :kword:`MOLCAS_WORKDIR` for more options.

  .. xmldoc:: <KEYWORD MODULE="ENVIRONMENT" NAME="WorkDir" APPEAR="WorkDir" KIND="STRING" LEVEL="ADVANCED">
              %%Keyword: WorkDir <advanced>
              <HELP>
              This variable can be used to specify directly the directory where all files
              that molcas creates are placed. See MOLCAS_WORKDIR for more options.
              </HELP>
              </KEYWORD>

:variable:`CurrDir`
  This variable corresponds to the location of the input, and it is used as
  a default location for all output files, generated by |molcas| modules.

  .. xmldoc:: <KEYWORD MODULE="ENVIRONMENT" NAME="CurrDir" APPEAR="CurrDir" KIND="STRING" LEVEL="ADVANCED">
              %%Keyword: CurrDir <advanced>
              <HELP>
              This variable corresponds to the location of the input, and it is used as
              a default location for all output files, generated by molcas modules.
              </HELP>
              </KEYWORD>

:variable:`MOLCAS`
  This variable indicates the location of |molcas|. The default version of Molcas
  to be used is specified at file :file:`.Molcas/molcas`, located at user HOME directory.

  .. xmldoc:: <KEYWORD MODULE="ENVIRONMENT" NAME="MOLCAS" APPEAR="MOLCAS" KIND="STRING" LEVEL="ADVANCED">
              %%Keyword: MOLCAS <advanced>
              <HELP>
              This variable indicates the location of molcas. The default version of Molcas
              to be used is specified at file .Molcas/molcas, located at user HOME directory.
              </HELP>
              </KEYWORD>

:variable:`MOLCAS_NPROCS`
  This variable should be used to run |molcas| code in parallel. It defines the
  number of computational units (cores or nodes) which will be used.

  .. xmldoc:: <KEYWORD MODULE="ENVIRONMENT" NAME="MOLCAS_NPROCS" APPEAR="MPI Processes" KIND="STRING" LEVEL="BASIC">
              %%Keyword: MOLCAS_NPROCS <basic>
              <HELP>
              This variable should be used to run molcas code in parallel. It defines the
              number of computational units (cores or nodes) which will be used.
              </HELP>
              </KEYWORD>

:variable:`MOLCAS_MEM`
  .. compound::

    This environment variable controls the size (soft limit) of the
    work array utilized in the programs that offer dynamic memory.
    It is specified in Megabytes, i.e.

    :command:`MOLCAS_MEM=256; export MOLCAS_MEM`

    will assign 256MB for the working arrays.
    It is also possible to use Gb (Tb) to specify memory in Gb or Tb.

  * MOLCAS_MEM is undefined --- The default amount of memory (1024MB),
    will be allocated for the work arrays.
  * MOLCAS_MEM is defined but nonzero --- This amount of memory
    will be allocated.

  See also :kword:`MOLCAS_MAXMEM`.

  .. xmldoc:: <KEYWORD MODULE="ENVIRONMENT" NAME="MOLCAS_MEM" APPEAR="MOLCAS_MEM (Mb)" KIND="INT" LEVEL="BASIC">
              %%Keyword: MOLCAS_MEM <basic>
              <HELP>
              This environment variable controls the size (soft limit) of the
              work array (in Mb) utilized in the programs that offer dynamic memory.
              It is also possible to set up memory in Gb, e.g. 2Gb
              </HELP>
              </KEYWORD>

The complete list of |molcas|-related environment variables:

.. class:: variablelist

:variable:`MOLCAS_COLOR`
  By default molcas uses markup characters in the output. To overwrite, set the key to NO.

  .. xmldoc:: <KEYWORD MODULE="ENVIRONMENT" NAME="MOLCAS_COLOR" APPEAR="Use markup in the output" KIND="CHOICE" LIST="----,NO" LEVEL="ADVANCED">
              %%Keyword: MOLCAS_COLOR <advanced>
              <HELP>
              By default molcas uses markup characters in the output.
              To overwrite, set the key to NO
              </HELP>
              </KEYWORD>

:variable:`MOLCAS_NPROCS`
  See above

:variable:`MOLCAS_DEBUGGER`
  This variable can be set to the name of debugger (or another code) which will be used on top of
  molcas executables. The option is useful for tracing an error in the code

  .. xmldoc:: <KEYWORD MODULE="ENVIRONMENT" NAME="MOLCAS_DEBUGGER" APPEAR="Debugger" KIND="STRING" LEVEL="ADVANCED">
              %%Keyword: MOLCAS_DEBUGGER <advanced>
              <HELP>
              This variable can be set to the name of debugger (or another code) which will be used on top of
              molcas executables. The option is useful for tracing an error in the code
              </HELP>
              </KEYWORD>

:variable:`MOLCAS_DISK`
  The value of this variable is used to split large files into a set of
  smaller datasets, as many as are needed (max. 20 subsets). It is specified
  in Megabytes, for instance, :command:`MOLCAS_DISK=1000; export MOLCAS_DISK`,
  and the following rules apply:

  * MOLCAS_DISK is undefined --- The program modules will ignore this
    option and the file size limit will be defined by your hardware
    (2 GBytes for 32-bit machines).
  * MOLCAS_DISK=0 (zero) --- The programs will assume a file size limit
    of 2 GBytes (200GBytes on 64-bit machines).
  * MOLCAS_DISK is defined but nonzero --- The files will be limited to
    this value (approximately) in size.

  .. xmldoc:: <KEYWORD MODULE="ENVIRONMENT" NAME="MOLCAS_DISK" APPEAR="MOLCAS_DISK" KIND="INT" LEVEL="ADVANCED">
              %%Keyword: MOLCAS_DISK <advanced>
              <HELP>
              The value of this variable is used to split large files into a set of
              smaller datasets, as many as are needed (max. 20 subsets).
              </HELP>
              </KEYWORD>

:variable:`MOLCAS_ECHO_INPUT`
  An environment variable to control echoing of the input.
  To suppress print level, set MOLCAS_ECHO_INPUT to ``NO``.

  .. xmldoc:: <KEYWORD MODULE="ENVIRONMENT" NAME="MOLCAS_ECHO_INPUT" APPEAR="Echo input" KIND="CHOICE" LIST="----,NO" LEVEL="ADVANCED">
              %%Keyword: MOLCAS_ECHO_INPUT <advanced>
              <HELP>
              An environment variable to control echoing of the input.
              To suppress print level, set MOLCAS_ECHO_INPUT to 'NO'.
              </HELP>
              </KEYWORD>

:variable:`MOLCAS_FIM`
  Activates the Files In Memory I/O layer. See section :ref:`MT:sec:fim` for more details.
  *Note that this setting is available only in MOLCAS compiled without Global
  Arrays.*

  .. warning::

     This feature is not available in OpenMolcas.

  .. xmldoc:: <KEYWORD MODULE="ENVIRONMENT" NAME="MOLCAS_FIM" APPEAR="FiM" KIND="CHOICE" LIST="----,YES" LEVEL="ADVANCED">
              %%Keyword: MOLCAS_FIM <advanced>
              <HELP>
              Activate the Files in Memory I/O layer
              </HELP>
              </KEYWORD>

:variable:`MOLCAS_INPORB_VERSION`
  Selects the version used for writing orbital files (`$Project.ScfOrb`, `$Project.RasOrb`, etc.).
  The value should be a version number such as ``1.0`` or ``2.2``.
  If the version is not known, the default (usually latest) version will be used.

  .. xmldoc:: <KEYWORD MODULE="ENVIRONMENT" NAME="MOLCAS_INPORB_VERSION" APPEAR="INPORB version" KIND="REAL" LEVEL="BASIC">
              %%Keyword: MOLCAS_INPORB_VERSION <basic>
              <HELP>
              Selects the version used for writing orbital files.
              The value should be a version number such as 1.0 or 2.2.
              </HELP>
              </KEYWORD>

:variable:`MOLCAS_KEEP_WORKDIR`
  If set to NO |molcas| will remove scratch area after a calculation.
  This setting can be overwritten by running :command:`molcas` with flag :command:`-clean`.

  .. xmldoc:: <KEYWORD MODULE="ENVIRONMENT" NAME="MOLCAS_KEEP_WORKDIR" APPEAR="Keep WorkDir" KIND="CHOICE" LIST="NO" LEVEL="BASIC">
              %%Keyword: MOLCAS_KEEP_WORKDIR <basic>
              <HELP>
              If set to NO molcas will remove scratch area after a calculation.
              This setting can be overwritten by running molcas with flag -clean:
              || molcas -clean input
              </HELP>
              </KEYWORD>

:variable:`MOLCAS_LICENSE`
  An environment which specifies the directory with |molcas| license file :file:`license.dat`.
  The default value of this variable is directory :file:`.Molcas/` in user home directory.

  .. xmldoc:: <KEYWORD MODULE="ENVIRONMENT" NAME="MOLCAS_LICENSE" APPEAR="License Directory" KIND="DIR" LEVEL="ADVANCED">
              %%Keyword: MOLCAS_LICENSE <advanced>
              <HELP>
              An environment which specifies the directory with molcas license file license.dat.
              The default value of this variable is directory .Molcas/ in user home directory.
              </HELP>
              </KEYWORD>

:variable:`MOLCAS_LINK`
  An environment variable to control information about linking of files.
  By default (MOLCAS_LINK is not set) only essential
  information about linking will be printed. To increase/decrease the
  print level, set MOLCAS_LINK to ``Yes``/``No``.

  .. xmldoc:: <KEYWORD MODULE="ENVIRONMENT" NAME="MOLCAS_LINK" APPEAR="Link information" KIND="CHOICE" LIST="YES,NO" LEVEL="ADVANCED">
              %%Keyword: MOLCAS_LINK <advanced>
              <HELP>
              An environment variable to control information about linking of files.
              By default (MOLCAS_LINK is not set) only essential
              information about linking will be printed. To increase/decrease the
              print level, set MOLCAS_LINK to 'Yes'/'No'.
              </HELP>
              </KEYWORD>

:variable:`MOLCAS_MAXITER`
  An environment variable to control maximum number of iterations in DO WHILE loop.

  .. xmldoc:: <KEYWORD MODULE="ENVIRONMENT" NAME="MOLCAS_MAXITER" APPEAR="Max Iter" KIND="INT" LEVEL="ADVANCED">
              %%Keyword: MOLCAS_MAXITER <advanced>
              <HELP>
              An environment variable to control maximum number of iterations in DO WHILE loop
              </HELP>
              </KEYWORD>

:variable:`MOLCAS_MAXMEM`
  An environment variable to set up a hard limit for allocated memory (in Mb).
  If is not specified, then it takes value of MOLCAS_MEM. Otherwise, the (MOLCAS_MAXMEM-MOLCAS_MEM)
  amount of RAM will be primarily used for keeping files in memory (FiM), or allocating Distributed Global Arrays.
  *Note that this setting is available only in MOLCAS compiled without Global Arrays.*

  .. xmldoc:: <KEYWORD MODULE="ENVIRONMENT" NAME="MOLCAS_MAXMEM" APPEAR="Max Memory" KIND="INT" LEVEL="ADVANCED">
              %%Keyword: MOLCAS_MAXMEM <advanced>
              <HELP>
              An environment variable to set up a hard limit for allocated memory (in Mb).
              If is not specified, then it takes value of MOLCAS_MEM. Otherwise, the
              (MOLCAS_MAXMEM-MOLCAS_MEM) amount of RAM will be primarily used for keeping
              files in memory (FiM), or allocating Distributed Global Arrays. Note that
              this setting is available only in MOLCAS compiled without GA.
              </HELP>
              </KEYWORD>

:variable:`MOLCAS_MEM`
  See above.

:variable:`MOLCAS_MOLDEN`
  If MOLCAS_MOLDEN set to ``ON`` a Molden style input file will be generated regardless of the number of orbitals.

  .. xmldoc:: <KEYWORD MODULE="ENVIRONMENT" NAME="MOLCAS_MOLDEN" APPEAR="MOLDEN output" KIND="CHOICE" LIST="ON,OFF" LEVEL="ADVANCED">
              %%Keyword: MOLCAS_MOLDEN <advanced>
              <HELP>
              If MOLCAS_MOLDEN set to 'ON' a Molden style input file will be generated regardless of the number of orbitals.
              </HELP>
              </KEYWORD>

:variable:`MOLCAS_NEW_WORKDIR`
  If set to YES |molcas| will never reuse files in scratch area.
  This setting can be overwritten by running :command:`molcas` with flag :command:`-old`:
  :command:`molcas -old input`

  .. xmldoc:: <KEYWORD MODULE="ENVIRONMENT" NAME="MOLCAS_NEW_WORKDIR" APPEAR="Use new WorkDir" KIND="CHOICE" LIST="YES" LEVEL="BASIC">
              %%Keyword: MOLCAS_NEW_WORKDIR <basic>
              <HELP>
              If set to YES molcas will never reuse files in scratch area.
              This setting can be overwritten by running molcas with flag -old:
              || molcas -old input
              </HELP>
              </KEYWORD>

:variable:`MOLCAS_OUTPUT`
  This variable can alter the default directory for extra output files,
  such as orbitals files, molden files, etc.
  If set, |molcas| will save output files to the specified directory.
  The directory name can be set in the form of absolute PATH, or
  relative PATH (related to the submit directory).
  A special value ``WORKDIR`` will keep all output files in WorkDir.
  A special value ``NAME`` will create a subdirectory with a name of Project.
  If the variable is not set, all output files will be copied or moved
  to the current directory. Default value can be forced by MOLCAS_OUTPUT=PWD.

  .. xmldoc:: <KEYWORD MODULE="ENVIRONMENT" NAME="MOLCAS_OUTPUT" APPEAR="Output Directory" KIND="CHOICE" LIST="????,PWD,?DIR,NAME,WORKDIR" LEVEL="BASIC">
              %%Keyword: MOLCAS_OUTPUT <basic>
              <HELP>
              This variable can alter the default directory for extra output files,
              such as orbitals files, molden files, etc.
              If set, molcas will save output files to the specified directory.
              The directory name can be set in the form of absolute PATH, or
              relative PATH (related to submit directory)
              A special value 'WORKDIR' will keep all output files in WorkDir.
              A special value 'NAME' will create a subdirectory with a name of Project.
              If the variable is not set, all output files will be copied or moved
              to the current directory.
              Default value can be forced by MOLCAS_OUTPUT=PWD
              </HELP>
              </KEYWORD>

:variable:`MOLCAS_PRINT`
  MOLCAS_PRINT variable controls the level of output. The value could be numerical or mnemonic:
  SILENT (0), TERSE (1), NORMAL (2), VERBOSE (3), DEBUG (4) and INSANE (5).

  .. xmldoc:: <KEYWORD MODULE="ENVIRONMENT" NAME="MOLCAS_PRINT" APPEAR="Print level" KIND="CHOICE" LIST="----,0:SILENT,1:TERSE,2:NORMAL,3:VERBOSE,4:DEBUG,5:INSANE" LEVEL="BASIC">
              %%Keyword: MOLCAS_PRINT <basic>
              <HELP>
              MOLCAS_PRINT variable controls the level of output. The value could be numerical or mnemonic:
              SILENT (0), TERSE (1), NORMAL (2), VERBOSE (3), DEBUG (4) and INSANE (5).
              </HELP>
              </KEYWORD>

:variable:`MOLCAS_PROJECT`
  If set to value NAME, |molcas| will use the prefix of the input file
  as a project name. Otherwise, it set a project name for the calculation.
  If set to the value NAMEPID, the Project name still will be constructed
  from the name of input file, however, the name of scratch area will
  be random.

  .. xmldoc:: <KEYWORD MODULE="ENVIRONMENT" NAME="MOLCAS_PROJECT" APPEAR="Project Policy" KIND="CHOICE" LIST="????,NAME,NAMEPID" LEVEL="BASIC">
              %%Keyword: MOLCAS_PROJECT <basic>
              <HELP>
              If set to the value NAME, molcas will use the prefix of the input file
              as a project name,
              If set to the value NAMEPID, the Project name still will be constructed
              from the name of input file, however, the name of scratch area will
              be random
              </HELP>
              </KEYWORD>

:variable:`MOLCAS_PROPERTIES`
  If MOLCAS_PROPERTIES is set to ``LONG`` properties with the individual MO contributions will be listed.

  .. xmldoc:: <KEYWORD MODULE="ENVIRONMENT" NAME="MOLCAS_PROPERTIES" APPEAR="Property output" KIND="CHOICE" LIST="SHORT,LONG" LEVEL="ADVANCED">
              %%Keyword: MOLCAS_PROPERTIES <advanced>
              <HELP>
              If MOLCAS_PROPERTIES is set to 'LONG' properties with the individual MO contributions will be listed.
              </HELP>
              </KEYWORD>

:variable:`MOLCAS_REDUCE_PRT`
  If set to NO, print level in DO WHILE loop is not reduced.

  .. xmldoc:: <KEYWORD MODULE="ENVIRONMENT" NAME="MOLCAS_REDUCE_PRT" APPEAR="Verbose input in do loops" KIND="CHOICE" LIST="----,NO" LEVEL="ADVANCED">
              %%Keyword: MOLCAS_REDUCE_PRT <advanced>
              <HELP>
              If set to NO, print level in DO WHILE loop is not reduced
              </HELP>
              </KEYWORD>

:variable:`MOLCAS_REDUCE_NG_PRT`
  If set to NO, print level in :program:`numerical_gradient` loop is not reduced.

  .. xmldoc:: <KEYWORD MODULE="ENVIRONMENT" NAME="MOLCAS_REDUCE_NG_PRT" APPEAR="Verbose input in numerical_gradient" KIND="CHOICE" LIST="----,NO" LEVEL="ADVANCED">
              %%Keyword: MOLCAS_REDUCE_NG_PRT <advanced>
              <HELP>
              If set to NO, print level in numerical_gradient is not reduced.
              </HELP>
              </KEYWORD>

:variable:`MOLCAS_SAVE`
  This variable can alter the default filenames for output files.
  If not set (default), all files will overwrite old files.
  If set to ``INCR`` all output files will get an incremental
  filenames.
  If set to ``ORIG`` --- an existent file will be copied with
  an extension ``.orig``

  .. xmldoc:: <KEYWORD MODULE="ENVIRONMENT" NAME="MOLCAS_SAVE" APPEAR="Save files as" KIND="CHOICE" LIST="----,INCR,ORIG" LEVEL="BASIC">
              %%Keyword: MOLCAS_SAVE <basic>
              <HELP>
              This variable can alter the default filenames for output files.
              If not set (default), all files will overwrite old files.
              If set to 'INCR' all output files will get an incremental
              filenames.
              If set to 'ORIG' - an existent file will be copied with
              an extension '.orig'
              </HELP>
              </KEYWORD>

:variable:`MOLCAS_TIME`
  If set, switch on timing information for each module

  .. xmldoc:: <KEYWORD MODULE="ENVIRONMENT" NAME="MOLCAS_TIME" APPEAR="Timing" KIND="CHOICE" LIST="----,YES" LEVEL="ADVANCED">
              %%Keyword: MOLCAS_TIME <advanced>
              <HELP>
              If set, switch on timing information for each module
              </HELP>
              </KEYWORD>

:variable:`MOLCAS_TIMELIM`
  Set up a timelimit for each module (in minutes). By default, the maximum
  execution time is set to unlimited. *Note that this setting is available only
  in MOLCAS compiled without Global Arrays.*

  .. xmldoc:: <KEYWORD MODULE="ENVIRONMENT" NAME="MOLCAS_TIMELIM" APPEAR="Time Limit" KIND="INT" LEVEL="ADVANCED">
              %%Keyword: MOLCAS_TIMELIM <advanced>
              <HELP>
              Set up a timelimit for each module (in minutes). By default, the maximum
              execution time is set to unlimited. Note that this setting is available only
              in MOLCAS compiled without Global Arrays.
              </HELP>
              </KEYWORD>

:variable:`MOLCAS_TRAP`
  If MOLCAS_TRAP set to ``OFF`` |molcas| modules will continue to be executed,
  even if a non-zero return code was produced.

  .. xmldoc:: <KEYWORD MODULE="ENVIRONMENT" NAME="MOLCAS_TRAP" APPEAR="Trap on Error" KIND="CHOICE" LIST="----,OFF" LEVEL="ADVANCED">
              %%Keyword: MOLCAS_TRAP <advanced>
              <HELP>
              If MOLCAS_TRAP set to 'OFF' molcas modules will continue to be executed,
              even if a non-zero return code was produced.
              </HELP>
              </KEYWORD>

:variable:`MOLCAS_WORKDIR`
  A parent directory for all scratch areas. It can be set to an
  absolute PATH (recommended), to a relative PATH, or to a special value PWD
  (to use current directory for scratch files)

  .. xmldoc:: <KEYWORD MODULE="ENVIRONMENT" NAME="MOLCAS_WORKDIR" APPEAR="Top WORKDIR" KIND="CHOICE" LIST="?DIR" LEVEL="BASIC">
              %%Keyword: MOLCAS_WORKDIR <basic>
              <HELP>
              A parent directory for all scratch areas. It can be set to an
              absolute PATH (recommended), to a relative PATH, or to a special value PWD
              (to use current directory for scratch files)
              </HELP>
              </KEYWORD>

User can customize his installation by adding MOLCAS environment variable into :file:`molcasrc` file.

Another way of customizing |molcas| is to use prologue and epilogue scripts. If user created a file
:file:`prologue` in :file:`$HOME/.Molcas` directory it will be executed (as :command:`./prologue`) before |molcas| calculation
starts. :file:`epilogue` in :file:`$HOME/.Molcas` directory will be executed at the end of calculation.
Files :file:`module.prologue` and :file:`module.epilogue` contains commands executing before and after
each executable molcas module. These files may use internal |molcas| variables, such as
:variable:`$Project`, :variable:`$WorkDir`, :variable:`$MOLCAS_MEM`, etc. Note that prologue/epilogue scripts should be executable.
For debug purposes, the location of prologue and epilogue files can be set by :variable:`$MOLCAS_LOGUE_DIR` variable.

Example:

:file:`prologue`: ::

  echo Calculation of $Project input will start at `date`

:file:`module.prologue`: ::

  echo Running module $MOLCAS_CURRENT_PROGRAM at $WorkDir

.. xmldoc:: </MODULE>

.. xmldoc:: <MODULE NAME="COMMENT" LEVEL="HIDDEN">
               <KEYWORD MODULE="COMMENT" NAME="UNDEFINED" APPEAR="Unrecognized Content" KIND="STRINGS" LEVEL="BASIC" />
            </MODULE>
