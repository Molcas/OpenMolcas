.. _UG\:sec\:EMIL:

General input structure. EMIL commands
======================================

.. xmldoc:: %%Description:
            EMIL commands

This is a general guide to the input structure of the programs in the
|molcas| program system. All programs conform to the same conventions
except where explicitly stated otherwise.

.. TODO: Note nested markup is currently not supported
   (|molcas| is not substituted if included in the italic text)

The programs are driven by keywords, which are either used without
further information, or followed by additional specifications on the
line(s) following the keyword, and is normally numeric in nature.
*All numerical inputs are read in free format, note that in general*
|molcas| *will not be able to process lines longer than 120 characters.*
The keywords can be given in mixed case (both upper and lower case are
allowed). In the input stream you can
insert comment lines anywhere, except between a keyword and the
following additional specifications, with a comment line identified by
an asterisk (*) in the first position on the line.

Most codes look at the first 4 characters of the keyword and ignores the
rest.
The entries in the lists of keywords below follow the
standard that the significant characters are in upper case and larger
than the nonsignificant characters.
This do not imply that the keywords have to be typed in upper case;
they can be typed freely in mixed case.

.. compound::

  All inputs begin with a name of the program followed by the keywords: ::

    &PROGRAM
    * here follows the keywords

  where PROGRAM is the name of the |molcas| module. The input listing is finished
  when a new program name, precede by the symbol &, is found (or the end of
  file or an EMIL command).

The following is an example of a list of keywords common to most of the
programs:

.. class:: keywordlist

:kword:`TITLe`
  This keyword starts the reading of title line.
  The following line is treated as title line.

The programs only decode the first four characters of a keyword
(except otherwise specifically indicated). For clarity it is however
recommended to write the full keyword name. The keywords can be typed freely
in upper, lower or mixed case.

An example for an input file used to run the :program:`SCF` program follows: ::

  &SCF
  Title
   Water molecule. Experimental equilibrium geometry
  * The symmetries are: a1, b2, b1 and a2.
  Occupied
  3 1 1 0
  * The ivo keyword prepares virtual orbitals for MCSCF.
  Ivo

Interpretation of |molcas| input is performed by :file:`molcas.exe`.
The internal language used by :file:`molcas.exe` is EMIL (Extended Molcas Input
Language).
It includes two different types of input commands:

* Sections with |molcas| input.
* EMIL commands (a line started with > character)

|molcas| input
--------------

EMIL allows to write |molcas| input in a more compact way:
user can omit :variable:`&END`, as well as a compulsory (in previous versions of |molcas|) keyword :kword:`End of input`.
As soon as a new module (or EMIL command) is requested in
a user input, the input for the module is terminated.

Also, it is possible to separate lines by ``;`` sign, or by ``=`` sign (to create
a pair ``keyword=value``). In some rare occasions signs ``;`` and ``=`` are used in
the input for a |molcas| module. In order to keep these symbols unchanged, user can
mark a part of an input, containing these symbols, by EMIL commands ``>> verbatim`` and
``>> end verbatim``.

.. compound::

  It means that the input: ::

    &SCF &END
    CHARGE
     1
    End of input
    &ALASKA &END
    End of input
    &SLAPAF &END
    End of input

  could be written as: ::

    &SCF; CHARGE=1
    &ALASKA; &SLAPAF

User can comment parts of input, by using ``*`` at the beginning of line,
or use C-style comments (``/* ... */``) to comment several lines.

In a rare occasion user might want to execute a UNIX command from
the input. It is important to understand that not all UNIX commands
can be understood and interpreted by EMIL. Also, EMIL should know
where to execute a command -- only at the master node, or for all
parallel tasks. In the past, EMIL supports the usage of commands
started from an exclamation mark, or with command :command:`UNIX`.
To avoid confusions, the serial execution of a command is now related to :command:`SHELL`,
and the parallel execution to :command:`EXEC`.

.. _UG\:sec\:emil_commands:

EMIL commands
-------------

.. compound::

  EMIL commands can be written in a short form: ::

    > KEY [VALUE]

  or in a nice form: ::

    >>>>>>>>>>  KEY  [VALUE]  <<<<<<<<<

EMIL commands are not case sensitive, but the variables used in commands must be written in upper case.
Also, it is important to place spaces in between elements (words) in the commands.

Here is a list of EMIL commands:

.. xmldoc:: <EMIL>

.. class:: commandlist

:command:`>> EXPORT A=B`
  a command to set environment variable A to value B

  .. xmldoc:: <COMMAND NAME="EXPORT" APPEAR="Export" FORMAT="EXPORT %s = %s" SHOWVALUE="TRUE">
              %%Keyword: EXPORT <basic>
              <HELP>
              A command to export environment variable in a form A=B
              </HELP>
              </COMMAND>

:command:`>> EXIT`
  a command to terminate execution.
  An optional value for this command is the return code (default value is 0)

  .. xmldoc:: <COMMAND NAME="EXIT" APPEAR="Exit" FORMAT="EXIT">
              %%Keyword: EXIT <basic>
              <HELP>
              A command to terminate execution.
              An optional value for this command is the return code (default value is 0)
              </HELP>
              </COMMAND>

:command:`>> INCLUDE file`
  a command to include a file into the input
  A compulsory value for this command is the filename.

  .. xmldoc:: <COMMAND NAME="INCLUDE" APPEAR="Include file" FORMAT="INCLUDE %s" LEVEL="HIDDEN">
              %%Keyword: INCLUDE <basic>
              <HELP>
              A command to include file

                >> INCLUDE filename
              </HELP>
              </COMMAND>

:command:`>> FILE file`
  A compulsory value for this command is the filename. A command to inline a
  file in the input file. The file will be extracted into WorkDir before the
  start of the calculation. The end of file should be marked as :command:`EOF` command.
  Not that the file is only created in the master process WorkDir, if the slaves
  need access to it, you'll need to use the :command:`COPY` command (see below).
  All files specified with :command:`FILE` are created at the beginning of the calculation,
  regardless of their placement in the input.

  .. xmldoc:: <COMMAND NAME="FILE" APPEAR="Inline file" FORMAT="FILE %s" LINK_ANCHOR="EOF" CONTENT="ANY">
              %%Keyword: FILE <basic>
              <HELP>
              A command to inline a file in the input file. The file will be extracted into
              WorkDir before the start of the calculation

                >> FILE filename
                ...
                >> EOF
              </HELP>
              </COMMAND>

:command:`>> EOF`
  A command to close inlined file.

  .. xmldoc:: <COMMAND NAME="EOF" APPEAR="EOF" FORMAT="EOF">
              %%Keyword: EOF <basic>
              <HELP>
              A command to close inlined file.
              </HELP>
              </COMMAND>

:command:`>> SHELL`
  a command to execute a unix command in serial.

  .. xmldoc:: <COMMAND NAME="SHELL" APPEAR="Serial shell" FORMAT="SHELL %s" SHOWVALUE="TRUE">
              %%Keyword: SHELL <basic>
              <HELP>
              A command to define a unix command to be executed in serial.
              </HELP>
              </COMMAND>

:command:`>> EXEC`
  a command to execute a unix command in parallel.

  .. xmldoc:: <COMMAND NAME="EXEC" APPEAR="Parallel shell" FORMAT="EXEC %s" SHOWVALUE="TRUE">
              %%Keyword: EXEC <basic>
              <HELP>
              A command to define a unix command to be executed in parallel.
              Note that any special characters are ignored.
              </HELP>
              </COMMAND>

:command:`>> LINK`
  a command to make a link between two files, located in WorkDir. The command is similar to
  ``!ln -s FILE1 FILE2`` but in parallel environment it is executed in all WorkDirs. The command assumes that
  FILE1 does exist, and FILE2 does not at the moment. >>LINK -FORCE allows
  to link a file which does not exist. User should avoid the usage of LINK commands in the input.

  .. xmldoc:: <COMMAND NAME="LINK" APPEAR="Link" FORMAT="LINK %s %s" FILE_INDEX="0" VALUES="????|????.OR.ITER" SHOWVALUE="TRUE">
              %%Keyword: LINK <basic>
              <HELP>
              A command to link two files located in WorkDir.
              The command is similar to '!ln -s FILE1 FILE2' but in parallel environment it is
              executed in all WorkDirs. The command assumes that FILE1 does exist, and FILE2 does not.
              >>LINK -FORCE allows to link a file FILE1 which does not exist at the moment. User should
              avoid the usage of LINK commands in the input.
              </HELP>
              </COMMAND>

:command:`>> COPY`
  a command to make a copy. The command is similar to ``!cp -f
  /path/to/FILE1 FILE2`` but can be used also in a parallel environment, in which case it
  will take the source file and distribute to the work directories of all
  processes. The destination must be located in the work directory. Note that
  EMIL command does not allow to use masks in the command. If FILE1 does not
  exist, the command returns an error code.

  .. xmldoc:: <COMMAND NAME="COPY" APPEAR="Copy to all" FORMAT="COPY %s %s" FILE_INDEX="0" SHOWVALUE="TRUE">
              %%Keyword: COPY <basic>
              <HELP>
              A command to copy one file to another. The command is similar to '!cp -f
              /path/to/FILE1 FILE2' but can be used also in a parallel environment, in which case it
              will take the source file and distribute to the work directories of all
              processes. The destination must be located in the work directory. Note that
              EMIL command does not allow to use masks in the command. If FILE1 does not
              exist, the command returns an error code.
              </HELP>
              </COMMAND>

:command:`>> CLONE`
  a command to make a clone copy of a file, doing a local copy on
  all slaves if parallel. It is mostly used internally, e.g. to distribute an input
  file to all WorkDirs.

  .. xmldoc:: <COMMAND NAME="CLONE" APPEAR="Copy on all" FORMAT="CLONE %s %s" FILE_INDEX="0" SHOWVALUE="TRUE">
              %%Keyword: CLONE <basic>
              <HELP>
              A command to make a clone copy of a file, doing a local copy on
              all slaves if parallel. It is mostly used internally, e.g. to distribute an input
              file to all WorkDirs.
              </HELP>
              </COMMAND>

:command:`>> COLLECT`
  A command to copy one file to another, collecting files on
  slaves and put them on the master if parallel. It is mostly used internally, e.g.
  to collect output files.

  .. xmldoc:: <COMMAND NAME="COLLECT" APPEAR="Collect from slaves" FORMAT="COLLECT %s %s" FILE_INDEX="0" SHOWVALUE="TRUE">
              %%Keyword: COLLECT <basic>
              <HELP>
              A command to copy one file to another, collecting files on slaves and put them on the master if parallel
              </HELP>
              </COMMAND>

:command:`>> SAVE`
  A command to copy one file to another, only on the master if parallel

  .. xmldoc:: <COMMAND NAME="SAVE" APPEAR="Copy on master" FORMAT="SAVE %s %s" FILE_INDEX="0" SHOWVALUE="TRUE">
              %%Keyword: SAVE <basic>
              <HELP>
              A command to copy one file to another, only on the master if parallel
              </HELP>
              </COMMAND>

:command:`>> RM`
  a command to delete a file. The command is similar to ``!rm
  FILE`` but can be used also in parallel environment. Note that EMIL command
  does not allow to use masks in the command. An attempt to remove non existent
  file leads to an error. It is possible to use -FORCE flag to allow deleting
  of non-existent file.

  .. xmldoc:: <COMMAND NAME="RM" APPEAR="Delete" FORMAT="RM %s" FILE_INDEX="0" SHOWVALUE="TRUE">
              %%Keyword: RM <basic>
              <HELP>
              A command to delete a file. The command is similar to '!rm
              FILE' but can be used also in parallel environment. Note that EMIL command
              does not allow to use masks in the command. An attempt to remove non existent
              file leads to an error. It is possible to use -FORCE flag to allow deleting
              of non-existent file.
              </HELP>
              </COMMAND>

:command:`>> EVAL A=B`
  evaluate a numerical value

  .. xmldoc:: <COMMAND NAME="EVAL" APPEAR="Evaluate" FORMAT="EVAL %s = %s">
              %%Keyword: EVAL <basic>
              <HELP>
              A command to evaluate numerical expression, e.g. eval A=$A+1
              </HELP>
              </COMMAND>

Keywords to organize loops in input, and execute modules conditionally:

.. class:: commandlist

:command:`>> DO WHILE`
  a command to start a loop. The loop should be terminated by SLAPAF or LOOP module,
  followed by ENDDO command

  .. xmldoc:: <COMMAND NAME="DO" APPEAR="Do loop" FORMAT="DO WHILE" LINK_ANCHOR="ENDDO">
              %%Keyword: DO WHILE <basic>
              <HELP>
              A command to start a loop. The loop should be terminated by SLAPAF or LOOP module,
              followed by ENDDO command
              </HELP>
              </COMMAND>

:command:`>> DO GEO`
  a command to start a special loop for geometry optimization
  with constrained internal coordinates. The loop should be terminated by
  ENDDO command. (See documentation for :program:`GEO` for more details.)

  .. xmldoc:: <COMMAND NAME="DOGEO" APPEAR="Geo loop" FORMAT="DO GEO" LINK_ANCHOR="ENDDO">
              %%Keyword: DO GEO <basic>
              <HELP>
              A command to start a constrained geometry optimization loop. The loop should be
              terminated by an ENDDO command
              </HELP>
              </COMMAND>

:command:`>> FOREACH A in (B, C, D)`
  a command to loop when the value of A is in the comma or space separated list.
  The list also can be written in the format ``From .. To``. Note that variable in the loop must be uppercased.

  .. xmldoc:: <COMMAND NAME="FOREACH" APPEAR="Foreach loop" FORMAT="FOREACH %s IN ( %s )" LINK_ANCHOR="ENDFOREACH" SHOWVALUE="TRUE">
              %%Keyword: FOREACH <basic>
              <HELP>
              A command to loop when the value of A is in the comma or space separated list.

                >> foreach A in (B,C,D)

              The list also can be written in the format "From .. To".
              Note that variable in the loop must be uppercased.
              </HELP>
              </COMMAND>

  .. xmldoc:: <COMMAND NAME="ENDFOREACH" FORMAT="END FOREACH" />

:command:`>> ENDDO`
  a command to finish the loop. If last module (before ENDDO command) returns
  1 --- the loop will be executed again (if number of iterations is less than MAXITER).
  If the return code is equal to 0 the loop will be terminated.

  .. xmldoc:: <COMMAND NAME="ENDDO" FORMAT="END DO">
              %%Keyword: ENDDO <basic>
              <HELP>
              A command to finish the loop.
              </HELP>
              </COMMAND>

:command:`>> IF ( ITER = N )`
  a command to make conditional execution of modules/commands on iteration N (N possibly could be a space separated list)

:command:`>> IF ( ITER NE N )`
  a command to skip execution of modules/commands on iteration N

:command:`>> IF ( ITER != N )`
  same as above

  .. xmldoc:: <COMMAND NAME="IF" APPEAR="Condition" FORMAT="IF ( %s %s %d )" VALUES="ITER.OR.????|=.OR.!=" LINK_ANCHOR="ENDIF" SHOWVALUE="TRUE">
              %%Keyword: IF ITER <basic>
              <HELP>
              A command to make conditional execution of modules/commands on iteration N
              </HELP>
              </COMMAND>

:command:`>> IF ( $VAR = N )`
  a command to make conditional execution if $VAR value equals to N (if statement terminated by ENDIF command)

:command:`>> IF ( $VAR = N ) GOTO JUMP`
  a command to make conditional goto to a label JUMP

:command:`>> IF ( -FILE file )`
  test for existence of a file

  .. xmldoc:: <COMMAND NAME="IFGOTO" APPEAR="Conditional jump" FORMAT="IF ( %s %s %d ) GOTO %s" VALUES="ITER.OR.????|=.OR.!=" LINK_ANCHOR="LABEL" LINK_VALUE_INDEX="3" SHOWVALUE="TRUE">
              %%Keyword: IF  <advanced>
              <HELP>
              A command to make conditional execution.
              Allowed syntax:

              * IF ( $VAR = 7 )           (IF statement terminated by ENDIF)
              * IF ( $VAR = 7 ) GOTO JUMP (jump to label JUMP)
              * IF ( -FILE file )         (test for existence of a file, terminated by ENDIF)
              </HELP>
              </COMMAND>

:command:`>> LABEL JUMP`
  a command to define a label. Note! Only forward jumps are allowed.

  .. xmldoc:: <COMMAND NAME="LABEL" APPEAR="Label" FORMAT="LABEL %s" LINK_TARGET="0" SHOWVALUE="TRUE">
              %%Keyword: LABEL <advanced>
              <HELP>
              A command to define a label. Note! Only forward jumps are allowed.
              </HELP>
              </COMMAND>

:command:`>> ENDIF`
  terminate :command:`IF` block. Note nested if's are not allowed.

  .. xmldoc:: <COMMAND NAME="ENDIF" FORMAT="END IF">
              %%Keyword: ENDIF  <basic>
              <HELP>
              Terminate IF block. Note nested if's are not allowed.
              </HELP>
              </COMMAND>

EMIL interpreter automatically stops calculation if a module returns a returncode
higher than 0 or 1. To force the interpretor to continue calculation even if a
returncode equal to 16 (which is a return code for non-convergent calculation) one
should set environment variable MOLCAS_TRAP=\ ``OFF``.

:program:`SLAPAF` returns a special return code in the case of converged (non converged) geometry.
So, to organize a structure calculation one should place the call to
:program:`SLAPAF` as a last statement of loop block. The summary of geometry optimization
convergence located in a file :file:`$Project.structure`.
The programs following a geometry optimization will automatically
assume the optimized geometry and wave function. Any new :program:`SEWARD`
calculation after an optimization (minimum or transition state) will
disregard the input coordinates and will take the geometry optimized at previous step.

It is also possible to use a special dummy program :program:`LOOP` to organize
infinite loops, or loops terminated by the counter (set by MOLCAS_MAXITER)

Keyword :kword:`SET` is obsolete and should be changed to :kword:`EXPORT`.

Verbatim input.

If an input for a module must contain special symbols, such as ``;`` or ``=``, user can
mark a corresponding part of the input by EMIL command :command:`VERBATIM`

.. class:: commandlist

:command:`>> VERBATIM <<`
  start verbatim input

  .. xmldoc:: <COMMAND NAME="VERBATIM" APPEAR="Verbatim" FORMAT="VERBATIM" LEVEL="HIDDEN" LINK_ANCHOR="ENDVERBATIM" CONTENT="ANY" />
              %%Keyword: VERBATIM <advanced>
              Start verbatim input

:command:`>> END VERBATIM <<`
  finish verbatim input

  .. xmldoc:: <COMMAND NAME="ENDVERBATIM" FORMAT="END VERBATIM" LEVEL="HIDDEN" />
              %%Keyword: END VERBATIM <advanced>
              Finish verbatim input

.. xmldoc:: </EMIL>

.. .. xmldoc:: <MODULE NAME="COMMENT" LEVEL="HIDDEN">
                 <KEYWORD MODULE="COMMENT" NAME="UNDEFINED" APPEAR="Unrecognized Content" KIND="STRINGS" LEVEL="BASIC" />
               </MODULE>

   .. xmldoc:: <MODULE NAME="DEMO">
                 <KEYWORD MODULE="DEMO" NAME="CHECKBOX" APPEAR="Check Appear" KIND="SINGLE" />
                 <KEYWORD MODULE="DEMO" NAME="INTTEXT" KIND="INT" REQUIRE="CHECKBOX" DEFAULT_VALUE="1" />
                 <KEYWORD MODULE="DEMO" NAME="REALTEXT" KIND="REAL" MIN_VALUE="0" MAX_VALUE="100.0" />
                 <KEYWORD MODULE="DEMO" NAME="STRINGTEXT" KIND="STRING" />
                 <KEYWORD MODULE="DEMO" NAME="INTTABLE" KIND="INTS" SIZE="3" />
                 <KEYWORD MODULE="DEMO" NAME="REALTABLE" KIND="REALS" SIZE="3" />
                 <KEYWORD MODULE="DEMO" NAME="MULTILINESTRING" KIND="STRINGS" SIZE="3" />
                 <KEYWORD MODULE="DEMO" NAME="MANYLINESTRING" KIND="STRINGS" />
                 <KEYWORD MODULE="DEMO" NAME="COMBOBOX" KIND="CHOICE" LIST="CHO1,CHO2,CHO3" />
                 <KEYWORD MODULE="DEMO" NAME="INTTABLE_COMPUTED" KIND="INTS_COMPUTED" SIZE="2" />
                 <KEYWORD MODULE="DEMO" NAME="COMPUTED_REALTABLE" KIND="REALS_COMPUTED" SIZE="3" />
                 <KEYWORD MODULE="DEMO" NAME="LOOKUP_INTTABLE" KIND="INTS_LOOKUP" SIZE="NSYM" />
                 <KEYWORD MODULE="DEMO" NAME="FILELOAD" KIND="FILE" />
                 <SELECT MODULE="DEMO" NAME="SELECTION" CONTAINS="SEL1,SEL2,SEL3">
                   <KEYWORD MODULE="DEMO" NAME="SEL1" KIND="REALS" SIZE="2" EXCLUSIVE="SEL2,SEL3" />
                   <KEYWORD MODULE="DEMO" NAME="SEL2" KIND="INTS" SIZE="3" EXCLUSIVE="SEL1,SEL3" />
                   <KEYWORD MODULE="DEMO" NAME="SEL3" KIND="SINGLE" EXCLUSIVE="SEL1,SEL2" />
                 </SELECT>
                 <GROUP MODULE="DEMO" NAME="BOXGROUP" KIND="BOX" WINDOW="INPLACE">
                   <KEYWORD MODULE="DEMO" NAME="CHECKBOX1" KIND="SINGLE" EXCLUSIVE="CHECKBOX2,CHECKBOX3" />
                   <KEYWORD MODULE="DEMO" NAME="CHECKBOX2" KIND="SINGLE" EXCLUSIVE="CHECKBOX1,CHECKBOX3" />
                   <KEYWORD MODULE="DEMO" NAME="CHECKBOX3" KIND="SINGLE" EXCLUSIVE="CHECKBOX1,CHECKBOX2" />
                   <KEYWORD MODULE="DEMO" NAME="REALTEXT1" KIND="REAL" MIN_VALUE="0" MAX_VALUE="100.0" />
                 </GROUP>
                 <GROUP MODULE="DEMO" NAME="BLOCKGROUP" KIND="BLOCK" WINDOW="POPUP">
                   <KEYWORD MODULE="DEMO" NAME="CHECKBOX4" KIND="SINGLE" EXCLUSIVE="CHECKBOX5,CHECKBOX6" />
                   <KEYWORD MODULE="DEMO" NAME="CHECKBOX5" KIND="SINGLE" EXCLUSIVE="CHECKBOX4,CHECKBOX6" />
                   <KEYWORD MODULE="DEMO" NAME="CHECKBOX6" KIND="SINGLE" EXCLUSIVE="CHECKBOX4,CHECKBOX5" />
                   <KEYWORD MODULE="DEMO" NAME="REALTEXT2" KIND="REAL" MIN_VALUE="0" MAX_VALUE="100.0" />
                 </GROUP>
                 <GROUP MODULE="DEMO" NAME="RADIOGROUP" KIND="RADIO" WINDOW="TAB">
                   <KEYWORD MODULE="DEMO" NAME="CHECKBOX7" KIND="SINGLE" />
                   <KEYWORD MODULE="DEMO" NAME="CHECKBOX8" KIND="SINGLE" />
                   <KEYWORD MODULE="DEMO" NAME="REALTEXT3" KIND="REAL" MIN_VALUE="0" MAX_VALUE="100.0" />
                 </GROUP>
               </MODULE>

Below are different input examples.

The first example shows the procedure to perform first a :program:`CASSCF`
geometry optimization of the water molecule, then a numerical hessian
calculation on the optimized geometry, and later to make a :program:`CASPT2`
calculation on the optimized geometry and wave function. Observe that
the position of the :program:`SLAPAF` inputs controls the data required
for the optimizations.

.. extractfile:: ug/EMIL.loop.input

  *
  *    Start Structure calculation
  *
  >>EXPORT MOLCAS_MAXITER=50
   &GATEWAY
  coord
  $MOLCAS/Coord/Water.xyz
  BASIS = ANO-S

  >>>>>>>>>>>>> Do while <<<<<<<<<<<<
   &SEWARD
  >>>>>>>> IF ( ITER = 1 ) <<<<<<<<<<
   &SCF
  >>>>>>> ENDIF <<<<<<<<<<<<<<<<<<<<

   &RASSCF
  Title
   H2O ANO(321/21).
  Nactel   = 6  0  0
  Spin     = 1
  Inactive = 1  0  0  0
  Ras2     = 3  1  0  2

   &ALASKA; &SLAPAF

  >>>>>>>>>>>>> ENDDO <<<<<<<<<<<<<<

   &CASPT2
  Maxit = 20
  Lroot = 1
   &GRID_IT

Another example demonstrate a possibility to use loops. SCF module
will be called twice --- first time with BLYP functional,
second time with B3LYP functional.

.. extractfile:: ug/EMIL.foreach.input

  *------------------------------------------------------
   &GATEWAY
  coord
  $MOLCAS/Coord/C2H6.xyz
  basis
  ANO-S-VDZ
  group
  y xz
  *------------------------------------------------------
   &SEWARD
  Title
  Ethane DFT test job
  *------------------------------------------------------
  >>foreach DFT in (BLYP, B3LYP )
   &SCF ; KSDFT = $DFT
  >>enddo
  *------------------------------------------------------

The next examples calculates HF energy for the several structures:

.. extractfile:: ug/EMIL.foreach.coord.sample

  * modify coordinates in place
  >>foreach DIST in (1.0, 2.0, 20.0)
   &GATEWAY
  Coord
  2
  hydrogen molecule
  H 0 0 0
  H $DIST 0 0
  BASIS= ANO-S-MB
  GROUP= C1
   &SEWARD
   &SCF
  UHF
  SCRAMBLE=0.3
  >>enddo

  * incremental change of coordinates
  >>export DIST=1.0
  >>foreach L in ( 1 .. 3 )
  >>eval DIST=$DIST+0.1
   &GATEWAY
  Coord
  2
  hydrogen molecule
  H 0 0 0
  H $DIST 0 0
  BASIS= ANO-S-MB
  GROUP= C1
   &SEWARD
   &SCF
  >>enddo

  * different coordinate files
  >> FILE H2001.xyz
  2

  H  0.300000000  0.000000000  0.000000000
  H -0.300000000  0.000000000  0.000000000

  >> FILE H2002.xyz
  2

  H  0.350000000  0.000000000  0.000000000
  H -0.350000000  0.000000000  0.000000000

  >> FILE H2003.xyz
  2

  H  0.400000000  0.000000000  0.000000000
  H -0.400000000  0.000000000  0.000000000

  >>foreach COO in ( 000, 001, 002)
   &GATEWAY
  Coord = H2$COO.xyz
  BASIS= ANO-S-MB
  GROUP= C1
   &SEWARD
   &SCF
  >>enddo

Use of shell parameters in input
--------------------------------

The |molcas| package allows the user to specify parts or variables in the the
input file with shell variables, which subsequently are dynamically defined
during execution time. **Note:** the shell variable names must be in upper
case. Find below a simple example where a part of the :math:`\ce{H2}` potential curve
is computed. First, the script used to run the calculation: ::

  #! /bin/sh
  #
  Home=`pwd` ;                     export Home
  Project=H2 ;                     export Project
  WorkDir=/tmp/$Project ;          export WorkDir
  #
  # Create workdir and cd to it
  #
  rm -fr $WorkDir
  mkdir $WorkDir
  #
  # Loop over distances
  #
  for R in 0.5 0.6 0.7 0.8 0.9 1.0
  do
     export R
     molcas $Home/$Project.input > $Home/$Project-$R-log 2> $Home/$Project-$R-err
  done
  #
  # Cleanup WorkDir
  #
  rm -fr $WorkDir

In this sh shell script we have arranged the call to the |molcas| package inside a loop over
the various values of the distances. This value is held by the variable $R which is
exported every iterations. Below is the input file used, note that the third cartesian
coordinate is the variable $R. ::

  &SEWARD
  Symmetry
   x y z
  Basis set
  H.sto-3g....
  H   0.000   0.000   $R
  End of basis
  End of input

  &SCF

Customization of molcas input
-----------------------------

EMIL interpretor supports templates (aliases) for a group of program calls
or/and keywords. The definition of these templates can be located in file :file:`alias`
located at |molcas| root directory, or at :file:`.Molcas/` directory.
The definition should be written in the following format: ``@name { sequence of EMIL commands }``.
In order to use the alias, the input should contain ``@name``.

.. compound::

  For example, user can define ::

    @DFTgeometry {
    >> DO WHILE
     &SEWARD
     &SCF; KSDFT=B3LYP;
     &SLAPAF
    >>ENDDO
    }

  and so, an input for geometry optimization can be written in the following form: ::

    &GATEWAY; Coord=Water.xyz; Basis = ANO-L-MB;
    @DFTgeometry

It is also possible to use parameters. In the alias file, possible parameters
have names: ``$1``, ``$2``, etc. up to 5 parameters.
In the user input an alias should be followed by parenthesis with comma separated list
of values.

Modifying the previous example: ::

  @DFTgeometry {
  >> DO WHILE
   &SEWARD
   &SCF; CHARGE=$1; KSDFT=$2;
   &SLAPAF
  >>ENDDO
  }

Input file now looks like: ::

  @DFTgeometry(0,B3LYP)
