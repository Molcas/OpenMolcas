.. index::
   single: Program; ESPF
   single: ESPF

.. _UG\:sec\:espf:

:program:`espf` (+ QM/MM interface)
===================================

.. only:: html

  .. contents::
     :local:
     :backlinks: none

.. _UG\:sec\:espf_description:

Description
-----------

.. xmldoc:: <MODULE NAME="ESPF">
            %%Description:
            <HELP>
            The ElectroStatic Potential Fitted (ESPF) method adds contributions
            to the one-electron hamiltonian to compute the interaction between
            the charge distribution and any external electrostatic potential,
            field, and field derivatives. This module may be used for QM/MM
            computations. Use of symmetry is forbidden.
            </HELP>

.. compound::

  The ElectroStatic Potential Fitted (ESPF) method adds contributions to the one-electron Hamiltonian for computing the interaction between the charge distribution in |molcas| and *any* external electrostatic potential, field, and field derivatives. The approximate interaction energy is expressed as:

  .. math:: \Delta E^{\text{ESPF}} = \left ( \sum_a \braopket {\Psi} {Q^a} {\Psi} + Z_a \right ) V^a

  with :math:`Q^a` a multipole-like operator whose matrix elements are fitted to the electron potential integrals (determined on a grid surrounding the QM atoms) and :math:`V^a` the *external* electrostatic potential (and derivatives) at nucleus :math:`a`. Both energy and gradient computations are available. A call to :program:`ESPF` right after :program:`SEWARD` is required to carry out such calculations.

*NOTE:* Always run :program:`SEWARD` + :program:`ESPF`. If not, very strange results may happen due to interactions counted twice or more!

*NOTE:* Symmetry is ignored since the external potential usually breaks the one given in :program:`GATEWAY`.

If no external potential is given, the :program:`ESPF` module can be used to compute atomic point charges fitted to the electrostatic potential produced by the nuclei and electrons.

ESPF and QM/MM
--------------

Whereas the ESPF method can be used standalone, it has been developed for hybrid quantum mechanics/molecular mechanics (QM/MM) computations, in which an extended molecular system is divided into two subsystems: the "active" center described by any QM method and its surroundings qualitatively treated with an empirical atomic forcefield. The current implementation can be used with either a modified version of the :program:`TINKER` program or with the :program:`GROMACS` program as MM code.

Using |molcas| together with :program:`Tinker`
..............................................

In order to obtain the modified :program:`TINKER` code, you must run the ":command:`molcas get_tinker`" command.

The current patched version of :program:`TINKER`\ [#fn1]_ is **6.3.3**.

.. [#fn1] https://dasher.wustl.edu/tinker/

*IMPORTANT:* The environment variable :variable:`TINKER` must point to the directory in which the :program:`TINKER` executable binaries are stored (usually in :file:`$MOLCAS/tinker/bin`).

The most convenient way to define (i) the QM and MM subsystems and (ii) which atoms are to be known by |molcas| (all the QM ones and some MM ones, see below) requires to simply add the keyword :kword:`TINKER` in :program:`GATEWAY`. This way, :program:`GATEWAY` will ask :program:`TINKER` to pass it all information needed.

Alternatively, the normal coordinate input in :program:`GATEWAY` can be used. For MM atoms that are to be known by |molcas|, if the atomic symbol is Xx, specify ``Xx...... / MM`` or ``Xx_MM`` in native or XYZ format, respectively. In this case you should make sure there is no mismatch between the |molcas| and :program:`TINKER` coordinates.

A final option is specifying :kword:`COORD`, with the atom labels to be used |molcas| and some dummy coordinates; and then :kword:`TINKER`, which will pick up the :program:`TINKER` coordinates, but keep the |molcas| labels. In order to use this combination, the :kword:`Expert` keyword must be specified before :kword:`TINKER` too.

Using |molcas| together with :program:`Gromacs`
...............................................

The interface to :program:`GROMACS` differs from the :program:`TINKER` interface in that the MM code is not run as a separate program but included in |molcas| as a library. In this way, the communication between the QM and MM codes is handled by simple function calls instead of using data files. The interface is automatically installed along with |molcas| provided that the :program:`GROMACS` library (currently a development version\ [#fn2]_) is available at configuration time\ [#fn3]_. Instructions how to install the :program:`GROMACS` library can be found at the official web site\ [#fn4]_. Make sure that the installation is done in double precision since this is the precision used by |molcas|. Also make sure to source the :program:`GROMACS` GMXR script in your shell startup file. otherwise the |molcas| configuration procedure will not be able to detect the relevant library path.

.. [#fn2] http://repo.or.cz/w/gromacs.git/shortlog/refs/heads/qmmm
.. [#fn3] Configuration with CMake requires the flag ``-D GROMACS=ON``
.. [#fn4] http://www.gromacs.org/

The recommended (and the only verified) approach of using the |molcas|/:program:`GROMACS` interface is to define the full QM+MM system in the :program:`GROMACS` input. The system definition can then be imported into |molcas| by adding the keyword :kword:`GROMACS` in :program:`GATEWAY` (see :numref:`UG:sec:gateway` for details). For efficiency reasons, the |molcas| part of the interface separates the MM subsystem into two different atom types: *inner* MM atoms and *outer* MM atoms. These are completely equivalent as far as interactions are concerned. However, whereas the coordinates of the inner MM atoms are stored and updated using |molcas| standard mechanism, the outer MM atoms are handled using a mechanism specifically designed with large systems in mind. The division into inner and outer MM atoms can be controlled with options to the :kword:`GROMACS` keyword in :program:`GATEWAY` (see :numref:`UG:sec:gateway`).

Please note that the |molcas|/:program:`GROMACS` interface is still under development and is currently provided for evaluation purposes only.

The QM/MM method
................

.. compound::

  The Hamiltonian of the full QM/MM system is divided into three terms:

  .. math:: H=H_{\text{QM}}+H_{\text{MM}}+H_{\text{QM/MM}}

  The first one describes the QM part as it would be *in vacuo*, the second one describes the surroundings using a classical MM forcefield and the last one deals with the interactions between the QM and the MM subsystems. In its usual formulation, the last term is (for :math:`q` point charges interacting with :math:`N` nuclei and :math:`n` electrons):

  .. math:: H_{\text{QM/MM}}=\sum_{a=1}^{q}\sum_{b=1}^{N}\frac{Q_{a}Z_{b}}{R_{ab}}-
          \sum_{a=1}^{q}\sum_{i=1}^{n}\frac{Q_{a}}{r_{ai}}+\sum_{a=1}^{q}\sum_{b=1}^{N}E_{ab}^{\text{vdw}}+
          E^{\text{bonded}}

The first two terms deal with the electrostatic interactions between the QM charge distribution and the MM electrostatic potential. In |molcas| the ESPF method is used for this purpose.
A short-range van der Waals term is added (van der Waals parameters are assigned to all the atoms --- both QM and MM). If the frontier between the two subsystems involves a bond, some empirical bonded terms may also be used. For the sake of simplicity, the standard MM parameters are kept unchanged for the MM atoms but should be modified (or calculated) for the QM atoms (e.g. it may be necessary to fit the QM van der Waals parameters).

The usual forcefields use the "1--4 condition" to separate the bonded interactions (stretching, bending, torsion) from the non-bonded ones (electrostatic and vdw). This means than the non-bonded potentials are applied only if atoms are separated by 3 bonds or more. As for the QM/MM interactions, this procedure is kept with the exception that all the QM atoms experience the electrostatic potential generated by *all* the MM point charges (the QM/MM frontier case is considered later).

*NOTE:* Starting with |molcasviii|, all MM point charges interact with the QM charge distribution using the ESPF method (at variance with previous |molcas| versions in which the few MM atoms defined in :program:`GATEWAY` were interacting directly with the QM electrons and nuclei).

Link atoms
..........

When no bonds are involved between the QM and the MM parts, the QM/MM frontier definition is obvious and only the electrostatic and vdw interactions are taken into account. However, if one or several chemical bonds exist, the definition of a smooth but realistic frontier is needed. Several schemes, more or less sophisticated, have been proposed. In the current implementation, only the most basic one, the link atom (LA) approach is included. In the LA approach, each QM/MM bond that should be cut is saturated with a monovalent atom --- most often a hydrogen atom --- on the QM side. The position of a link atom is often restrained: frozen distance from the corresponding QM frontier atom and always on the segment defined by the two frontier atoms (Morokuma's method, selected by the :kword:`LAMOROKUMA` keyword).

From the macromolecular point of view, link atoms do not exist, i.e. they should not interact with the MM part. However, this leads to severe overpolarization of the frontier, due to unbalanced interactions. Hence interactions between the link atoms and the MM potential is kept. To remove problems that may arise from too strong interactions between a link atom and the closest MM point charges, these point charges may be spread in the MM neighborhood. For instance, in a protein, this procedure is mainly justified if the MM frontier atom is an :math:`\alpha` carbon (Amber- or Charmm-typed forcefields usually set these point charges close to zero).

Geometry optimization --- microiterations
.........................................

In a QM/MM geometry optimization job, a |molcas| step costs as hundreds of :program:`TINKER` or :program:`GROMACS` steps. Thus it is very convenient to use the microiteration technique, that is, converging the MM subsystem geometry every |molcas| step. In the case of :program:`TINKER`, this is requested in the :program:`TINKER` keyword file, whereas if :program:`GROMACS` is used, it is requested directly in :program:`ESPF`. In order to improve the optimization convergence, an improved QM/MM Hessian can be built in :program:`SLAPAF` using the :kword:`RHIDDEN` keyword (note that adding the keyword :kword:`CARTESIAN` may help too).

.. TEMPORARILY REMOVED.
   The :program:`TINKER` package incorporates several polarizable forcefields, eg Amoeba, based on the induced dipoles approach. They can be selected for a QM/MM calculation. In that case, both the QM wavefunction and the MM induced dipoles are converged self-consistently during the SCF procedure, similarly to reaction-field calculations. In case of a State-Average CASSCF calculation, you must use the :kword:`RFRoot` keyword (see :program:`RASSCF`) to select the electronic state which is electrostatically coupled to the polarizable surroundings.

.. _UG\:sec\:espf_dependencies:

Dependencies
------------

The :program:`ESPF` program depends on :program:`SEWARD` for modifying the core Hamiltonian matrix and on :program:`ALASKA` for computing the extra contributions to the gradient.

.. index::
   pair: Files; ESPF

.. _UG\:sec\:espf_files:

Files
-----

:program:`ESPF` will use the following input
files: :file:`RYSRW`, :file:`ABDATA`, :file:`RUNFILE`, :file:`ONEINT` (for more information see :numref:`UG:sec:files_list`).
In addition, :program:`ESPF` uses :file:`ESPFINP` (the ESPF input file) and :file:`SEWARINP` (the Seward input file).

Please note that the external potential can be given within a file, separated from the :program:`ESPF` input file.

In calculations using the |molcas|/:program:`GROMACS` interface, :program:`ESPF` will additionally need access to the :program:`GROMACS` tpr file.

Intermediate files
..................

All the intermediate files are related to the use of :program:`ESPF` together :program:`TINKER`. The files allow for communication between the :program:`ESPF` program and the MM code. |molcas| uses one file to pass the QM atoms coordinates and ESPF-derived point charges to :program:`TINKER`. :program:`TINKER` uses the same file to pass the external potential, the MM-only energy and gradient components to |molcas|.

.. class:: filelist

:file:`TINKER.LOG`
  The log file of the Tinker run.

:file:`$Project.xyz`
  The coordinate file for :program:`TINKER`.

:file:`$Project.key`
  The keyword file for :program:`TINKER`.

:file:`$Project.qmmm`
  The communication file between |molcas| and :program:`TINKER`.

Output files
............

.. class:: filelist

:file:`ONEINT`
  One-electron integral file generated by the :program:`SEWARD` program.

:file:`RUNFILE`
  Communication file for subsequent programs.

:file:`ESPF.DATA`
  Ascii file containing some specific informations needed for subsequent calls to the :program:`ESPF` module.

:file:`GMX.LOG`
  Logfile for the :program:`GROMACS` library routines.

.. _UG\:sec\:espf_input:

Input
-----

Below follows a description of the input to :program:`ESPF`.

In addition to the keywords and the comment lines the input may contain blank lines. The input for each module is preceded by its name like: ::

  &ESPF

Compulsory keywords

.. class:: keywordlist

:kword:`EXTErnal`
  Specify how the external potential is given. This keyword is compulsory in the first run of :program:`ESPF`. On the next line, one integer or a text string must be given:

  * One integer :math:`n` is given. If :math:`n` is 0, the next lines give the numbering, the values for the external potential, the field and field gradients for each atom. If :math:`n` is greater than 0, the :math:`n` next lines specify the sources of the external potential, each line gives three cartesian coordinates, one point charge, and (optionally) three dipole components. If Å is used as the length unit, the :kword:`ANGSTROM` keyword must be given right after :math:`n`.
  * The :kword:`NONE` word means that no external potential is given. Accordingly, the :program:`ESPF` module will compute the atomic point charges (and optionally dipoles) deriving from the electrostatic potential due to all electrons and nuclei.
  * The word is :kword:`TINKER`, which means that the current job is a QM/MM job using the |molcas|/:program:`TINKER` interface. Accordingly the external potential will be computed directly by :program:`TINKER`. Note that :program:`TINKER` requires at least two input files, ending with .xyz (coordinates) and .key (keywords). These files must share the name of the current |molcas| project. Optionally, you can add the :kword:`MULLIKEN` or :kword:`LOPROP` keyword after :kword:`TINKER`: it indicates what kind of charges are passed to :program:`TINKER`. These charges may be used during the MM microiterations. If no keyword is given, the ESPF multipoles are selected.
  * The word is :kword:`GROMACS`, which means that the current job is a QM/MM job using the |molcas|/:program:`GROMACS` interface, with the external potential computed by :program:`GROMACS`. The binary input file read by :program:`GROMACS`, the so-called tpr file, must be named as ":file:`topol.tpr`" and must be manually copied to the working directory. As above, a second keyword on the same line can be used to select the type of multipoles sent to the MM code. Default is to use the ESPF multipoles.

  * Any other word. The following characters up to the next space are taken as a file name and the rest of the line is ignored. Instead, the full input (including the first line) is read from the specified file and must follow the syntax specified above.

  .. xmldoc::  <KEYWORD MODULE="ESPF" NAME="EXTERNAL" APPEAR="External potential" INPUT="REQUIRED" KIND="CUSTOM" LEVEL="BASIC">
               <ALTERNATE KIND="CHOICE" LIST="NONE,TINKER,GROMACS" />
               <ALTERNATE KIND="STRING" />
               %%Keyword: External potential <basic>
               %%Tested ##055 ##803
               <HELP>
               Specify how the external potential is given. Can be given inline or in another file.
               </HELP>
               </KEYWORD>

Optional keywords

.. class:: keywordlist

:kword:`TITLE`
  Title of the job.

  .. xmldoc::  <KEYWORD MODULE="ESPF" NAME="TITLE" KIND="STRING" LEVEL="BASIC">
               %%Keyword: Title <basic>
               <HELP>
               One line following this one is regarded as title.
               </HELP>
               </KEYWORD>

:kword:`MULTipoleorder`
  Multipolar order of the ESPF operators. For :program:`TINKER`, allowed values are 0 (charge-like) or 1 (charge- and dipole-like). For :program:`GROMACS`, only 0 is allowed. Default value is 0.

  .. xmldoc::  <KEYWORD MODULE="ESPF" NAME="MULT" APPEAR="Multipole order" KIND="INT" DEFAULT_VALUE="0" LEVEL="ADVANCED">
               %%Keyword: MultipoleOrder <basic>
               %%Tested ##803
               <HELP>
               Give the order of the ESPF operators. Only 0 (charge) or 1 (charge and
               dipole).
               </HELP>
               </KEYWORD>

:kword:`GRID`
  Modify the grid specifications. The grid is made of points belonging to molecular surfaces defined according to the van der Waals radii of each quantum atom. Two schemes are available. The first one is the GEPOL procedure, as implemented into the PCM SCRF method. The other one is called PNT and is the default. On the next line, first select the method with the GEPOL or PNT option. On the same line, one integer number and one real number are given if PNT is selected. The first one gives the maximum number of shells around the van der Waals surface of the quantum atoms. The second one gives the distance between the shells. Note that all points within the van der Waals envelope are discarded to avoid the penetration effects. Default values are 4 shells separated by 1 Å.
  Alternatively, if GEPOL is selected, the same line must contain 1 integer indicating the number of surfaces to be computed (must be < 6).

  .. xmldoc::  <KEYWORD MODULE="ESPF" NAME="GRID" APPEAR="Grid spec." KIND="STRING" LEVEL="ADVANCED">
               %%Keyword: Grid <advanced>
               %%Tested ##803
               <HELP>
               Modify the grid specifications.
               </HELP>
               </KEYWORD>

:kword:`SHOW`
  Requires the printing of the ESPF.DATA file.

  .. xmldoc::  <KEYWORD MODULE="ESPF" NAME="SHOW" KIND="SINGLE" LEVEL="ADVANCED">
               %%Keyword: Show <advanced>
               %%Tested NONE
               <HELP>
               Printing of the ESPF.DATA file.
               </HELP>
               </KEYWORD>

:kword:`LAMOrokuma`
  Activate the Morokuma scheme for scaling the link atom positions in a QM/MM calculation. Note that in the case of :program:`TINKER`, the scaling factor is currently hard-coded and is determined from the radii of the atoms involved in the QM/MM frontier bond. This differs from the :program:`GROMACS` interface in which this factor must be provided by the user through the :kword:`LINKATOMS` keyword in :program:`GATEWAY`.

  .. xmldoc::  <KEYWORD MODULE="ESPF" NAME="LAMO" KIND="SINGLE" LEVEL="ADVANCED">
               %%Keyword: LAMorokuma <advanced>
               %%Tested ##055
               <HELP>
               Set on the Morokuma's scheme for scaling the link atom positions.
               </HELP>
               </KEYWORD>

:kword:`MMITerations`
  Maximum number of microiterations used to optimize the outer MM atoms in a |molcas|/:program:`GROMACS` run. The default is 0, which disables microiterations and leaves the outer MM atoms frozen. For the :program:`TINKER` interface, microiterations are requested in the :program:`TINKER` keyword file.

  .. xmldoc::  <KEYWORD MODULE="ESPF" NAME="MMIT" KIND="INT" LEVEL="BASIC">
               %%Keyword: MMIterations <basic>
               <HELP>
               Maximum number of microiterations to optimize the MM subsystem (with Gromacs interface).
               </HELP>
               </KEYWORD>

:kword:`MMCOnvergence`
  Convergence threshold for the MM microiterations (:program:`GROMACS` only). The optimization of the (outer) MM atoms will stop when the maximum force component is smaller than this number, in atomic units. The default is 0.001 atomic units (50 kJ/mol/nm).

  .. xmldoc::  <KEYWORD MODULE="ESPF" NAME="MMCO" KIND="REAL" LEVEL="BASIC">
               %%Keyword: MMConvergence <basic>
               <HELP>
               Convergence for the MM microiterations (with Gromacs interface).
               </HELP>
               </KEYWORD>

.. xmldoc:: </MODULE>

Examples
--------

ESPF example
............

This is a typical input for the calculation of the energy and the gradient of a glycine molecule feeling the external potential of 209 TIP3P water molecules.

.. extractfile:: ug/ESPF.input

  &Gateway
  Basis set
  C.sto-3g.....
    C1   1.11820     0.72542    -2.75821 Angstrom
    C2   1.20948     0.66728    -1.25125 Angstrom
  End of basis
  Basis set
  O.sto-3g.....
    O1   2.19794     1.10343    -0.67629 Angstrom
  End of basis
  Basis set
  H.sto-3g.....
    H1   2.02325     1.18861    -3.14886 Angstrom
    H2   0.25129     1.31794    -3.04374 Angstrom
    H3   1.02458    -0.28460    -3.15222 Angstrom
  End of basis
  Basis set
  N.sto-3g.....
    N1   0.17609     0.12714    -0.61129 Angstrom
  End of basis
  Basis set
  C.sto-3g.....
    C3   0.09389    -0.01123     0.84259 Angstrom
    C4  -1.21244    -0.67109     1.28727 Angstrom
  End of basis
  Basis set
  O.sto-3g.....
    O2  -2.06502    -1.02710     0.48964 Angstrom
  End of basis
  Basis set
  H.sto-3g.....
    H4  -0.61006    -0.21446    -1.14521 Angstrom
    H5   0.92981    -0.61562     1.19497 Angstrom
    H6   0.16338     0.97444     1.30285 Angstrom
  End of basis
  Basis set
  N.sto-3g.....
    N2  -1.41884    -0.85884     2.57374 Angstrom
  End of basis
  Basis set
  H.sto-3g.....
    H7  -0.73630    -0.57661     3.25250 Angstrom
    H8  -2.28943    -1.29548     2.82140 Angstrom
  End of basis

  &seward

  &espf
  MultipoleOrder = 0
  External = 0
  1  -0.048 -0.002 -0.006 -0.001  0.007 -0.009  0.002 -0.001  0.001 -0.001
  2  -0.047 -0.002  0.001 -0.002  0.003  0.000 -0.004  0.000 -0.001  0.000
  3  -0.053  0.004  0.000 -0.011  0.002  0.002 -0.004  0.002  0.003 -0.007
  4  -0.046  0.011 -0.009 -0.001  0.006 -0.005 -0.001  0.003  0.003 -0.004
  5  -0.042 -0.016 -0.011 -0.006  0.005 -0.007  0.003 -0.004 -0.001 -0.005
  6  -0.050  0.000  0.008  0.001  0.006 -0.006  0.000 -0.002  0.000 -0.001
  7  -0.039 -0.008  0.001  0.000  0.001 -0.002  0.001 -0.001 -0.001 -0.001
  8  -0.032 -0.007 -0.002  0.004  0.002 -0.003  0.001 -0.002  0.002 -0.001
  9  -0.011 -0.009  0.004  0.001  0.002  0.000 -0.002 -0.001  0.001  0.001
  10  0.000 -0.011  0.003  0.004  0.001  0.002 -0.003  0.001 -0.001  0.001
  11 -0.028 -0.008  0.004 -0.001 -0.001 -0.002  0.002 -0.001  0.001 -0.002
  12 -0.026  0.003 -0.008  0.014  0.002 -0.001 -0.001 -0.008  0.006 -0.009
  13 -0.037 -0.008 -0.003  0.004 -0.007  0.007  0.000  0.001  0.007 -0.001
  14 -0.016 -0.007  0.007 -0.008  0.003  0.003 -0.006  0.000  0.002  0.002
  15 -0.025  0.003  0.012 -0.007  0.003 -0.001 -0.002 -0.006  0.005  0.009
  16 -0.010 -0.011  0.000 -0.014  0.001  0.007 -0.008  0.001  0.000 -0.001

  &scf
  Charge = 0

  &alaska

|molcas|/:program:`Tinker` example
..................................

A typical start for a QM/MM calculation with the |molcas|/:program:`TINKER` interface is given in the following input. It is quite general since all the information related to the QM and MM subsystem definitions are already included in the :program:`TINKER` key file.

.. extractfile:: ug/QMMM.input

  > EXPORT TINKER=$MOLCAS/tinker/bin_qmmm
  > COPY $PATH_TO/$Project.xyz $WorkDir/$Project.xyz
  > COPY $PATH_TO/$Project.key $WorkDir/$Project.key

   &Gateway
  Tinker
  Basis = STO-3G
  Group = Nosym

   &Seward

   &Espf
  External = Tinker
  LAMorok

This can be used, e.g. with the following :program:`TINKER` files. In this example, the asparate anion is cut into two pieces, the QM subsystem contains the end of the side-chain until the :math:`\beta` carbon atom. There is a link atom between the QM :math:`\beta` and MM :math:`\alpha` carbon atoms.

QMMM.xyz

.. extractfile:: ug/QMMM.xyz

  16  ASP
   1 N3    -0.040452    0.189961    0.173219   448     2     6    14    15
   2 CT    -0.011045   -0.060807    1.622395   449     1     3     7    11
   3 C      1.446535   -0.110535    2.028518   450     2     4     5
   4 O      1.902105    0.960982    2.409042   452     3
   5 O      2.137861   -0.898168    1.387158   452     3
   6 H      0.559257   -0.496270   -0.262338   451     1
   7 CT    -0.789906   -1.336520    1.982558   216     2     8    12    13
   8 C     -2.256402   -1.184505    1.571038   218     7     9    10
   9 O2    -2.460769   -0.949098    0.356151   219     8
  10 O2    -3.120135   -1.188969    2.465678   219     8
  11 H1    -0.478878    0.773493    2.145163   453     2
  12 HC    -0.356094   -2.194944    1.466324   217     7
  13 HC    -0.720511   -1.505463    3.058628   217     7
  14 H     -0.996208    0.061130   -0.151911   451     1
  15 H      0.304306    1.116522   -0.018698   451     1
  16 HLA   -0.283317   -0.506767    1.748300  2999     2     7

QMMM.key

.. extractfile:: ug/QMMM.key

  * Change $PATH_TO_TINKER
  parameters $PATH_TO_TINKER/params/amber99.prm
  QMMM 8
  QM -8 10 7 12 13
  MM 2
  LA 16
  * Add the atom type for the LA
  atom   2999    99    HLA     "Hydrogen Link Atom"        1      1.008     0
  charge -2  0.0
  charge -11 0.0
  QMMM-MICROITERATION ON

|molcas|/:program:`Gromacs` example
...................................

To be provided soon.
