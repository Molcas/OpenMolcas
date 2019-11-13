.. index::
   single: Program; Seward
   single: Seward

.. _UG\:sec\:seward:

:program:`seward`
=================

.. only:: html

  .. contents::
     :local:
     :backlinks: none

.. xmldoc:: <MODULE NAME="SEWARD">
            %%Description:
            <HELP>
            The Seward module generates one- and two-electron integrals needed
            by other programs. The input contains an additional and optional
            embedded input section for numerical quadrature options for the
            computation of integrals associated with DFT calculations. The
            embedded section starts and ends with the keywords "Grid Input"
            and "End of Grid", respectively. Keywords associated with this
            embedded input section is labelled "(NQ)" below.
            </HELP>

The :program:`SEWARD` module generates one- and two-electron integrals needed
by other programs. The two-electron integrals may optionally be
Cholesky decomposed. In addition, it will serve as the input parser
for parameters related to the specification of the quadrature grid
used in numerical integration in association with DFT and reaction
field calculations.

.. .. figure:: seward.*
      :width: 50%
      :align: center

      H. W. Seward, secretary of State in the Lincoln administration, who suggested and supervised
      the 1867 purchase of Alaska from tzar Russia. Price: 2 cents an acre.

In the following three subsections we will in detail describe the input parameters
for analytic integration, numerical integration, and relativistic operators.

.. ..include:: ../integrals.inc

Analytic integration
--------------------

Any conventional ab initio quantum chemical calculation starts by
computing overlap, kinetic energy, nuclear attraction and electron
repulsion integrals. These are used repeatedly to determine the
optimal wave function and the total energy of the system under
investigation. Finally, to compute various properties of the system
additional integrals may be needed, examples include multipole moments
and field gradients.

.. _UG\:sec\:seward_description:

Description
...........

:program:`SEWARD` is able to compute the following integrals:

* kinetic energy,

* nuclear attraction,

* two electron repulsion (optionally Cholesky decomposed),

* :math:`n`\th (default :math:`n`\=2) order moments (overlap, dipole moment, etc.),

* electric field (generated at a given point by all charges in the system),

* electric field gradients (spherical gradient operators),

* linear momentum (velocity),

* orbital angular momentum,

* relativistic mass-velocity correction (1st order),

* one-electron Darwin contact term,

* one-electron relativistic no-pair Douglas--Kroll,

* diamagnetic shielding,

* spherical well potential (Pauli repulsion),

* ECP and PP integrals,

* modified kinetic energy and multipole moment integrals
  (integration on a finite sphere centered at the origin) for use in
  the variational :math:`R`\-matrix approach,

* external field (represented by a large number of charges and dipoles),

* angular momentum products, and

* atomic mean-field integrals (AMFI) for spin--orbit coupling.

*Note that while* :program:`SEWARD` *compute these integrals the input to
select them and their settings are put in the input of* :program:`GATEWAY`.

In general these integrals will be written to a file, possibly in
the form of Cholesky vectors (two-electron integrals only). However,
:program:`SEWARD` can also compute the orbital contributions and total components of
these properties if provided with orbital coefficients and
orbital occupation numbers.

To generate the one- and two-electron integrals
:program:`SEWARD` uses two different integration schemes. Repulsion type integrals (two-electron
integrals, electric field integrals, etc.) are evaluated by
the reduced multiplication scheme of the Rys quadrature :cite:`seward`.
All other integrals are computed by the Gauss--Hermite quadrature.
:program:`SEWARD` use spherical Gaussians as basis functions,
the only exception to this are the diffuse/polarization functions
of the 6-31G family of basis sets.
The double coset :cite:`dcf` formalism is used to treat symmetry.
:program:`SEWARD` is especially designed to handle ANO-type basis sets, however, segmented basis
sets are also processed.

At present the following limitations are built into :program:`SEWARD`:

.. include:: ../limitations.inc

.. _UG\:sec\:seward_dependencies:

Dependencies
............

:program:`SEWARD` usually runs after program :program:`GATEWAY`. In the same time, any input used
in :program:`GATEWAY` can be placed into :program:`SEWARD` input. However, it is recommended to
specify all information about the molecule and the basis set in :program:`GATEWAY` input.

:program:`SEWARD` does normally not depend on any other code, except of :program:`GATEWAY`.
There are two exceptions however.
The first one is when :program:`SEWARD` is used as a property module. In these cases the file
:file:`INPORB` has to have been generated by a wave function code. The second case, which is
totally transparent to the user, is when :program:`SEWARD` picks up the new Cartesian coordinates
generated by :program:`SLAPAF` during a geometry optimization.

.. _UG\:sec\:seward_files:

Files
.....

Input Files
:::::::::::

Apart form the standard input file
:program:`SEWARD` will use the following input files: :file:`RYSRW`, :file:`ABDATA`,
:file:`RUNFILE`, :file:`INPORB` (for calculation of properties) (:numref:`UG:sec:files_list`).
In addition, :program:`SEWARD` uses the following files:

.. class:: filelist

:file:`BASLIB`
  The default directory for one-particle basis set information.
  This directory contains files which are part
  of the program system and could
  be manipulated by the user in accordance with the instructions in
  :numref:`UG:sec:the_basis_set_libraries` and following subsections.
  New basis set files can be added to this directory by the local
  |molcas| administrator.

:file:`QRPLIB`
  Library for numerical mass-velocity plus Darwin potentials (used for ECPs).

Output files
::::::::::::

In addition to the standard output file
:program:`SEWARD` may generate the following files:
:file:`ONEINT`, :file:`ORDINT`, :file:`CHVEC`, :file:`CHRED`, :file:`CHORST`,
:file:`CHOMAP`, :file:`CHOR2F` (:numref:`UG:sec:files_list`).

.. _UG\:sec\:seward_input:

Input
.....

Below follows a description of the input to :program:`SEWARD`.
Observe that if
nothing else is requested
:program:`SEWARD` will by default compute the overlap, the
dipole, the quadrupole, the nuclear attraction, the kinetic energy,
the one-electron Hamiltonian, and the two-electron repulsion
integrals.

The input for each module is preceded by its name like: ::

  &SEWARD

Argument(s) to a keyword, either individual or composed by several entries,
can be placed in a separated line or in the same line separated by a semicolon.
If in the same line, the first argument requires an equal sign after the
name of the keyword.

General keywords
::::::::::::::::

.. class:: keywordlist

:kword:`TITLe`
  One line of title card follows.

  .. xmldoc:: <KEYWORD MODULE="SEWARD" NAME="TITL" APPEAR="Title" KIND="CUSTOM" LEVEL="BASIC">
              <HELP>
              Enter one optional title card.
              </HELP>
              %%Keyword: Title <basic>
              One line of title card follows.
              </KEYWORD>

:kword:`TEST`
  :program:`SEWARD` will only process the input and generate a non-zero return code.

  .. xmldoc:: <KEYWORD MODULE="SEWARD" NAME="TEST" APPEAR="Test" KIND="SINGLE" LEVEL="BASIC">
              %%Keyword: Test <basic>
              <HELP>
              SEWARD will only process the input and generate a non-zero
              return code.
              </HELP>
              </KEYWORD>

:kword:`ONEOnly`
  :program:`SEWARD` will not compute the two-electron integrals.

  .. xmldoc:: <KEYWORD MODULE="SEWARD" NAME="ONEO" APPEAR="Only 1-el. Integrals" KIND="SINGLE" LEVEL="BASIC">
              %%Keyword: Oneonly <basic>
              <HELP>
              SEWARD will not compute the two-electron integrals.
              </HELP>
              </KEYWORD>

:kword:`NODKroll`
  :program:`SEWARD` will not compute Douglas--Kroll integrals.

  .. xmldoc:: <KEYWORD MODULE="SEWARD" NAME="NODK" APPEAR="Skip Douglas-Kroll Integrals" KIND="SINGLE" LEVEL="BASIC">
              %%Keyword: NoDK <basic>
              <HELP>
              SEWARD will not compute Douglas-Kroll integrals.
              </HELP>
              </KEYWORD>

:kword:`DIREct`
  Prepares for later integral-direct calculations. As with keyword
  :kword:`OneOnly`, :program:`SEWARD` will evaluate no two-electron integrals.

  .. xmldoc:: <KEYWORD MODULE="SEWARD" NAME="DIRE" APPEAR="Direct option" KIND="SINGLE" LEVEL="BASIC">
              %%Keyword: Direct <basic>
              <HELP>
              Prepares for later integral-direct calculations. As with keyword
              oneonly, SEWARD will evaluate no two-electron integrals.
              </HELP>
              </KEYWORD>

:kword:`EXPErt`
  Sets "expert mode", in which various default settings are
  altered. Integral-direct calculations will be carried out
  if the two-electron integral file is unavailable.

  .. xmldoc:: <KEYWORD MODULE="SEWARD" NAME="EXPE" APPEAR="Expert mode" KIND="SINGLE" LEVEL="ADVANCED">
              <HELP>
              Activates expert mode which will deactivate some checks. This will, for
              example, allow the used to combine relativistic and non-relativistic basis sets.
              </HELP>
              %%Keyword: Expert <advanced>
              Sets "expert mode", in which various default settings are
              altered. Integral-direct calculations will be carried out
              if the two-electron integral file is unavailable.
              </KEYWORD>

:kword:`CHOLesky`
  :program:`SEWARD` will Cholesky decompose the two-electron integrals using
  default configuration (in particular, the decomposition threshold is
  1.0d-4) of the decomposition driver.
  The decomposition threshold can be changed using keyword :kword:`THRC`.
  Default is to not decompose.

  .. xmldoc:: <GROUP MODULE="SEWARD" KIND="BOX" NAME="AUX" APPEAR="CD options" LEVEL="BASIC">

  .. xmldoc:: <KEYWORD MODULE="SEWARD" NAME="CHOL" APPEAR="Cholesky" KIND="SINGLE" EXCLUSIVE="RIJ,RIJK,RIC,RICD,LOW,MEDI,HIGH" LEVEL="BASIC">
              %%Keyword: Cholesky <basic>
              <HELP>
              Cholesky decompose the two-electron integrals using default settings.
              </HELP>
              </KEYWORD>

:kword:`1CCD`
  :program:`SEWARD` will Cholesky decompose the two-electron integrals using
  the one-center approximation.
  The decomposition threshold can be changed using keyword :kword:`THRC`.
  Default is to not decompose.

  .. xmldoc:: <KEYWORD MODULE="SEWARD" NAME="1CCD" APPEAR="1CCD" KIND="SINGLE" LEVEL="BASIC">
              %%Keyword: 1CCD <basic>
              <HELP>
              One-center Cholesky decomposition of the two-electron integrals.
              </HELP>
              </KEYWORD>

:kword:`THRCholesky`
  Specify decomposition threshold for Cholesky decomposition of two-electron integrals
  on the next line. Default value: 1.0d-4.

  .. xmldoc:: <KEYWORD MODULE="SEWARD" NAME="THRC" APPEAR="ThrCholesky" KIND="REAL" LEVEL="BASIC">
              %%Keyword: ThrCholesky <basic>
              <HELP>
              Specify decomposition threshold for Cholesky decomposition of two-electron integrals.
              </HELP>
              </KEYWORD>

:kword:`SPAN`
  Specify span factor (0 :math:`<` span :math:`\leq` 1) for Cholesky decomposition of two-electron integrals on the next line.
  Span=1 implies full pivoting, Span=0 means no pivoting. If the span factor is too low, numerical errors may cause
  negative diagonal elements and force the program to quit; if the span factor is too large, the execution time may
  increase. Default value: 1.0d-2.

  .. xmldoc:: <KEYWORD MODULE="SEWARD" NAME="SPAN" APPEAR="Span factor" KIND="REAL" LEVEL="BASIC">
              %%Keyword: Span <basic>
              <HELP>
              Specify span factor (in between, but not equal to, 0 and 1) for Cholesky decomposition of two-electron integrals on the next line.
              </HELP>
              </KEYWORD>

:kword:`LOW Cholesky`
  :program:`SEWARD` will Cholesky decompose the two-electron integrals using
  low accuracy (threshold 1.0d-4)
  configuration of the decomposition driver.
  Default is to not decompose.

  .. xmldoc:: <KEYWORD MODULE="SEWARD" NAME="LOW" APPEAR="Low" KIND="SINGLE" EXCLUSIVE="RIJ,RIJK,RIC,RICD,CHOL,MEDI,HIGH" LEVEL="BASIC">
              %%Keyword: Low <basic>
              <HELP>
              Cholesky decompose the two-electron integrals using low accuracy settings. Recommended.
              </HELP>
              </KEYWORD>

:kword:`MEDIum Cholesky`
  :program:`SEWARD` will Cholesky decompose the two-electron integrals using
  medium accuracy (threshold 1.0d-6)
  configuration of the decomposition driver.
  Default is to not decompose.

  .. xmldoc:: <KEYWORD MODULE="SEWARD" NAME="MEDI" APPEAR="Medium" KIND="SINGLE" EXCLUSIVE="RIJ,RIJK,RIC,RICD,CHOL,LOW,HIGH" LEVEL="BASIC">
              %%Keyword: Medium <basic>
              <HELP>
              Cholesky decompose the two-electron integrals using medium accuracy settings.
              </HELP>
              </KEYWORD>

:kword:`HIGH Cholesky`
  :program:`SEWARD` will Cholesky decompose the two-electron integrals using
  high-accuracy (threshold 1.0d-8)
  configuration of the decomposition driver.
  Default is to not decompose.

  .. xmldoc:: <KEYWORD MODULE="SEWARD" NAME="HIGH" APPEAR="High" KIND="SINGLE" EXCLUSIVE="RIJ,RIJK,RIC,RICD,CHOL,LOW,MEDI" LEVEL="BASIC">
              %%Keyword: High <basic>
              <HELP>
              Cholesky decompose the two-electron integrals using high accuracy settings.
              </HELP>
              </KEYWORD>

:kword:`LDF`
  Local Density Fitting: :program:`SEWARD` will fit each AO product using only auxiliary functions centered on the two parent atoms.
  Fitting coefficients are computed and stored for use in later modules (only non-hybrid KS-DFT at the moment).
  Requires an auxiliary basis set (generated by RICD in Gateway or externally defined).
  Can not be combined with Cholesky and (global) DF/RI. LDF is turned off by default.

  .. xmldoc:: <KEYWORD MODULE="SEWARD" NAME="LDF" APPEAR="Local Density Fitting" KIND="SINGLE" EXCLUSIVE="CHOL,LOW,MEDI,HIGH" LEVEL="BASIC">
              %%Keyword: LDF <basic>
              <HELP>
              Local Density Fitting using auxiliary functions centered on the two parent atoms of each AO product.
              </HELP>
              </KEYWORD>

:kword:`LDF1`
  Local Density Fitting using auxiliary functions centered on the two parent atoms of each AO product. Equivalent to keyword :kword:`LDF`.

  .. xmldoc:: <KEYWORD MODULE="SEWARD" NAME="LDF1" APPEAR="Local Density Fitting" KIND="SINGLE" EXCLUSIVE="CHOL,LOW,MEDI,HIGH" LEVEL="BASIC">
              %%Keyword: LDF1 <basic>
              <HELP>
              Local Density Fitting using auxiliary functions centered on the two parent atoms of each AO product. Equivalent to keyword LDF.
              </HELP>
              </KEYWORD>

:kword:`LDF2`
  LDF using auxiliary functions centered on the two parent atoms of each AO product as well as the two-center functions required to achieve the given target accuracy.

  .. xmldoc:: <KEYWORD MODULE="SEWARD" NAME="LDF2" APPEAR="Local Density Fitting with two-center auxiliary functions" KIND="SINGLE" EXCLUSIVE="CHOL,LOW,MEDI,HIGH" LEVEL="BASIC">
              %%Keyword: LDF2 <basic>
              <HELP>
              Local Density Fitting using auxiliary functions centered on the two parent atoms of each AO product as well as the two-center functions required to achieve the given target accuracy.
              </HELP>
              </KEYWORD>

:kword:`TARGet accuracy`
  Specify the target accuracy for LDF2 on the next line. Default value: threshold used to generate the aCD or acCD auxiliary basis set
  (or 1.0d-4 in case of externally defined auxiliary basis sets).

  .. xmldoc:: <KEYWORD MODULE="SEWARD" NAME="TARG" APPEAR="Target accuracy for LDF2" KIND="REAL" EXCLUSIVE="CHOL,LOW,MEDI,HIGH" LEVEL="BASIC">
              %%Keyword: TARG <basic>
              <HELP>
              Specify the target accuracy for LDF2 on the next line. Default value: threshold used to generate the aCD or acCD auxiliary basis set (or 1.0d-4 in case of externally defined auxiliary basis sets).
              </HELP>
              </KEYWORD>

:kword:`APTH`
  Specify the LDF/LDF2 atom pair prescreening threshold on the next line, thus defining which atom pairs are considered significant.
  Default value: the target accuracy or the cutoff threshold for primitive integral evaluation, whichever is smaller.

  .. xmldoc:: <KEYWORD MODULE="SEWARD" NAME="APTH" APPEAR="Atom pair threshold" KIND="REAL" EXCLUSIVE="CHOL,LOW,MEDI,HIGH" LEVEL="BASIC">
              %%Keyword: APTH <basic>
              <HELP>
              Specify the LDF/LDF2 atom pair prescreening threshold on the next line, thus defining which atom pairs are considered significant. Default value: the target accuracy or the cutoff threshold for primitive integral evaluation, whichever is smaller.
              </HELP>
              </KEYWORD>

:kword:`CLDF`
  Constrained LDF/LDF2: specify constraint order on the next line (-1 for unconstrained, 0 for charge constraint). Default: unconstrained LDF.

  .. xmldoc:: <KEYWORD MODULE="SEWARD" NAME="CLDF" APPEAR="Constrained LDF" KIND="INT" EXCLUSIVE="CHOL,LOW,MEDI,HIGH" LEVEL="BASIC">
              %%Keyword: CLDF <basic>
              <HELP>
              Constrained LDF/LDF2: specify constraint order on the next line (-1 for unconstrained, 0 for charge constraint). Default: unconstrained LDF.
              </HELP>
              </KEYWORD>

  .. xmldoc:: </GROUP>

:kword:`FAKE CD/RI`
  If CD/RI vectors are already available, :program:`SEWARD` will not redo work!

  .. xmldoc:: <KEYWORD MODULE="SEWARD" NAME="FAKE" APPEAR="Skip CD/RI vectors generation" KIND="SINGLE" LEVEL="BASIC">
              %%Keyword: FAKE <basic>
              <HELP>
              If CD/RI vectors are already available, SEWARD will not redo work!
              </HELP>
              </KEYWORD>

  .. :kword:`MOLPro`
       The integrals are normalized as in the program package
       :program:`MOLPRO`. The default is the normalization
       according to the double coset representatives (DCR) formalism
       for conventional calculations; for Cholesky decomposition, the
       default is the :kword:`MOLECULE` normalization (except when AMFI or
       Douglas--Kroll is specified, where the DCR formalism is again
       default).

  ..   .. xmldoc:: %%Keyword: MOLPRO <advanced>
                   The integrals are normalized as in the program package
                   MOLPRO. The default is the normalization
                   according to the double coset representatives (DCR) formalism
                   for conventional calculations; for Cholesky decomposition, the
                   default is the MOLECULE normalization (except when AMFI or
                   Douglas-Kroll is specified, where the DCR formalism is again
                   default).

  .. :kword:`MOLEcule`
       The integrals are normalized as in the integral program
       :program:`MOLECULE`. The default is the normalization
       according to the double coset representatives (DCR) formalism
       for conventional calculations; for Cholesky decomposition, the
       default is the :kword:`MOLECULE` normalization (except when AMFI or
       Douglas--Kroll is specified, where the DCR formalism is again
       default).

  ..   .. xmldoc:: %%Keyword: MOLECULE <advanced>
                   The integrals are normalized as in the integral program
                   MOLECULE. The default is the normalization
                   according to the double coset representatives (DCR) formalism
                   for conventional calculations; for Cholesky decomposition, the
                   default is the MOLECULE normalization (except when AMFI or
                   Douglas-Kroll is specified, where the DCR formalism is again
                   default).

  .. xmldoc:: <KEYWORD MODULE="SEWARD" NAME="MOLECULE" KIND="SINGLE" LEVEL="UNDOCUMENTED" />


  .. :kword:`MOLCas`
       The integrals are normalized according to the double coset
       representatives (DCR) formalism. This is the default for
       conventional calculations. The default for Cholesky decompositions
       is the :kword:`MOLECULE` normalization (except when AMFI or
       Douglas--Kroll is specified, where the DCR formalism is again
       default).

  ..   .. xmldoc:: %%Keyword: MOLCAS <advanced>
                   The integrals are normalized according to the double coset
                   representatives (DCR) formalism. This is the default
                   for conventional calculations; for Cholesky decomposition, the
                   default is the MOLECULE normalization (except when AMFI or
                   Douglas-Kroll is specified, where the DCR formalism is again
                   default).

:kword:`JMAX`
  The integer entry on the next line is the highest rotational quantum
  number for which
  :program:`SEWARD` will compute the rotational energy within the
  rigid rotor model. The default value is 5.

  .. xmldoc:: %%Keyword: JMAX <advanced>
              The integer entry on the next line is the highest rotational quantum
              number for which
              SEWARD will compute the rotational energy within the
              rigid rotor model. The default value is 5.

:kword:`SYMMetry`
  See the the description in the manual for the program :program:`GATEWAY`.

  .. :kword:`SKIP`
       The subsequent entry contains a sequence of integers, one for
       each irreducible representation. By setting an integer to a nonzero
       value all basis functions from the corresponding
       irreducible representation will be removed.
       The default is to skip no basis functions.

   ..  .. xmldoc:: %%Keyword: Skip <advanced>
                   The subsequent entry contains a sequence of integers, one for
                   each irreducible representation. By setting an integer to a nonzero
                   value all basis functions from the corresponding
                   irreducible representation will be removed.

  .. :kword:`NOTAbles`
       Indicates that no tables will be used for the the roots and weights
       of the Rys polynomials. The default is to use tables.
       This option is used to verify the data base for
       the calculation of roots and weights. Observe that the activation of
       this option will make the program perform rather slow.

  ..   .. xmldoc:: %%Keyword: Notables <advanced>
                   Indicates that no tables will be used for the the roots and weights
                   of the Rys polynomials. The default is to use tables.

:kword:`BASIs Set`
  See the the description in the manual for the program :program:`GATEWAY`.

:kword:`ZMAT`
  See the the description in the manual for the program :program:`GATEWAY`.

:kword:`XBAS`
  See the the description in the manual for the program :program:`GATEWAY`.

:kword:`XYZ`
  See the the description in the manual for the program :program:`GATEWAY`.

:kword:`NOGUessorb`
  Disable automatic generation of starting orbitals with the GuessOrb procedure.

  .. xmldoc:: <KEYWORD MODULE="SEWARD" NAME="NOGUESSORB" APPEAR="No guess orbitals" KIND="SINGLE" LEVEL="BASIC">
              %%Keyword: NOGUESSORB <basic>
              <HELP>
              Disable automatic generation of starting orbitals with the GuessOrb procedure.
              </HELP>
              </KEYWORD>

:kword:`NODElete`
  Do not delete any orbitals automatically.

  .. xmldoc:: <GROUP MODULE="SEWARD" KIND="BOX" NAME="DELOPT" APPEAR="Delete options" LEVEL="BASIC">

  .. xmldoc:: <KEYWORD MODULE="SEWARD" NAME="NODE" APPEAR="No orbitals deleted" KIND="SINGLE" LEVEL="BASIC">
              %%Keyword: NODElete <basic>
              <HELP>
              Do not delete any orbitals automatically.
              </HELP>
              </KEYWORD>

:kword:`SDELete`
  Set the threshold for deleting orbitals based on the eigenvalues of the overlap matrix.
  All eigenvalues with eigenvectors below this threshold will be deleted.
  If you want no orbitals deleted use keyword :kword:`NODElete`.

  .. xmldoc:: <KEYWORD MODULE="SEWARD" NAME="SDEL" APPEAR="S delete threshold" KIND="REAL" LEVEL="BASIC">
              %%Keyword: SDELete <basic>
              <HELP>
              Set the threshold for deleting orbitals based on the eigenvalues of the overlap matrix.
              </HELP>
              </KEYWORD>

:kword:`TDELete`
  Set the threshold for deleting orbitals based on the eigenvalues of the kinetic energy matrix.
  All eigenvalues with eigenvectors above this threshold will be deleted.
  If you want no orbitals deleted use keyword :kword:`NODElete`.

  .. xmldoc:: <KEYWORD MODULE="SEWARD" NAME="TDEL" APPEAR="T delete threshold" KIND="REAL" LEVEL="BASIC">
              %%Keyword: TDELete <basic>
              <HELP>
              Set the threshold for deleting orbitals based on the eigenvalues of the kinetic energy matrix.
              </HELP>
              </KEYWORD>

  .. xmldoc:: </GROUP>

:kword:`VERBose`
  Force :program:`SEWARD` to print a bit more verbose.

  .. xmldoc:: <KEYWORD MODULE="SEWARD" NAME="VERB" APPEAR="Verbose output" KIND="SINGLE" LEVEL="BASIC">
              %%Keyword: Verbose <basic>
              <HELP>
              Force SEWARD to print a bit more verbose.
              </HELP>
              </KEYWORD>

:kword:`EMBEdding`
  Reads in an embedding potential from a file. It can also write out the
  density and the electrostatic potential on a grid. It is a block keyword
  which *must* be ended with :kword:`ENDEmbedding`.

  .. xmldoc:: <GROUP MODULE="SEWARD" NAME="EMBE" APPEAR="External embedding potential" KIND="BLOCK" LEVEL="ADVANCED">
              %%Keyword: Embedding <advanced>
              <HELP>
              Reads in an embedding potential from a file and can output the density and ESP on a grid.
              </HELP>

  The subkeywords are:

  :kword:`EMBInput` ---
    Specifies the path to the file which contains the embedding potential
    (e.g. ``EMBI=myEmbPot.dat``). The file contains a potential given on a grid. It
    has the number of grid points in the first line. Then in five columns
    data for each grid point is given (x, y, z, weight of this grid point,
    value of the potential). Default is ``EMBPOT``.

  :kword:`OUTGrid` ---
    Specifies the path to a file containing a grid on which the output is
    produced. It is only needed if you want to have the data on a grid different
    from the one given in :kword:`EMBInput`. The columns for the potential and weights
    are not needed in this file (and also not read in).

  :kword:`WRDEnsity` ---
    Switches on the calculation of the final electron density on a grid. The
    output file path must be specified along with this keyword (e.g. ``WRDE=myDens.dat``).

  :kword:`WREPotential` ---
    Switches on the calculation of the electrostatic potential on a grid. The
    output file path must be specified along with this keyword (e.g. ``WREP=myESP.dat``).

  :kword:`ENDEmbedding` ---
    Ends the :kword:`EMBEdding` section. This keyword *must* be present.

  The :kword:`EMBEdding` feature is currently only supported by the SCF part of |molcas|.

  .. xmldoc:: </GROUP>

Keywords associated to one-electron integrals
:::::::::::::::::::::::::::::::::::::::::::::

.. class:: keywordlist

:kword:`MULTipoles`
  Entry which specifies the highest order of the
  multipole for which integrals will be generated. The default center
  for the dipole moment operator is the origin.
  The default center for the higher order operators is the
  center of the nuclear mass.
  The default is to do up to quadrupole moment integrals (2).

  .. xmldoc:: <GROUP MODULE="SEWARD" KIND="BOX" NAME="ONE" APPEAR="1-electron Integral options" LEVEL="BASIC">

  .. xmldoc:: <KEYWORD MODULE="SEWARD" NAME="MULTIPOLE" APPEAR="Set Max order of Multipole Moment Integrals" KIND="INT" LEVEL="ADVANCED" DEFAULT_VALUE="2">
              %%Keyword: Multipoles <basic>
              <HELP>
              Entry which specifies the highest order of the
              multipole for which integrals will be generated.
              For details consult the manual.
              </HELP>
              </KEYWORD>

:kword:`CENTer`
  This option is used to override the default selection of the origin
  of the multipole moment operators.
  On the first entry add an integer entry specifying the number of
  multipole moment operators for which the origin of expansion
  will be defined.
  Following this, one entry for each operator, the order of
  the multipole operator and the coordinates of the center (in au) of
  expansion are specified.
  The default is the origin for 0th and 1st multipoles, and the center of mass for higher-order multipoles.

  .. xmldoc:: <KEYWORD MODULE="SEWARD" NAME="CENT" KIND="REALS_COMPUTED" APPEAR="Multipole moments origins" SIZE="4" LEVEL="ADVANCED">
              %%Keyword: Center <advanced>
              <HELP>
              This option is used to override the default selection of the origin
              of the multipole moment operators.
              On the first entry specify the number of such modifications.
              On the following entries enter first the order of the multipole moment operators
              followed by the coordinates of the origin of the operator.
              </HELP>
              </KEYWORD>

:kword:`RELInt`
  Requests the computation of mass-velocity and one-electron Darwin
  contact term integrals for the calculation of a first order correction
  of the energy with respect to relativistic effects.

  .. xmldoc:: <KEYWORD MODULE="SEWARD" NAME="RELI" APPEAR="1st order rel. corr." KIND="SINGLE" LEVEL="ADVANCED">
              %%Keyword: Relint <basic>
              <HELP>
              Requests the computation of mass-velocity and one-electron Darwin
              contact term integrals for the calculation of a first order correction
              of the energy with respect to relativistic effects.
              </HELP>
              </KEYWORD>

  .. :kword:`BSSMethod`
       Request that the one-electron Hamiltonian include the scalar relativistic
       effects according to the so-called Barysz--Sadlej--Snijders transformation.

  ..   .. xmldoc:: %%Keyword: BSSmethod <basic>
                   Request that the one-electron Hamiltonian include the scalar relativistic
                   effects according to the so-called Barysz-Sadlej-Snijders transformation.

  .. xmldoc:: <KEYWORD MODULE="SEWARD" NAME="BSSMETHOD" KIND="SINGLE" LEVEL="UNDOCUMENTED" />

:kword:`RELAtivistic`
  Request arbitrary scalar relativistic Douglas--Kroll--Hess (DKH) correction to the one-electron Hamiltonian
  and the so-called picture-change correction to the property integrals (multipole moments
  and electronic potential related properties). An argument of the form ``RXXPyy`` follows.
  Here XX represents the order of the DKH correction to the one-electron Hamiltonian and
  yy the order of the picture-change correction. The character P denotes the parameterization
  used in the DKH procedure.
  The possible parametrizations P of the unitary transformation used
  in the DKH transformation supported by |molcas| are:

  .. container:: list

    ``P=O``: Optimum parametrization (OPT)

    ``P=E``: Exponential parametrization (EXP)

    ``P=S``: Square-root parametrization (SQR)

    ``P=M``: McWeeny parametrization (MCW)

    ``P=C``: Cayley parametrization (CAY)

  Hence, the proper keyword for the 4th order relativistically corrected one-electron
  Hamiltonian and 3rd order relativistically corrected
  property integrals in the EXP parameterization would read as R04E03. If yy is larger than XX it is set to
  XX. If yy is omitted it will default to same value as XX. Recommended orders and parametrization is
  R02O02.
  Since the EXP parameterization employs a fast algorithm, it is
  recommended for high order DKH transformation.

  .. xmldoc:: <KEYWORD MODULE="SEWARD" NAME="RELATIVISTIC" APPEAR="Relativistic correction order and parametrization (RXXPyy)" KIND="CUSTOM" LEVEL="BASIC">
              %%Keyword: RELAtivistic <basic>
              <HELP>
              Request arbitrary scalar relativistic Douglas-Kroll-Hess (DKH) correction to the one-electron Hamiltonian
              and the so-called picture-change correction to the property integrals (multipole moments
              and electronic potential related properties). An argument of the form RXXPyy follows
              Here XX represents the order of the DKH correction to the one-electron Hamiltonian and
              yy the order of the picture-change correction. The character P denotes the parameterization
              used in the DKH procedure.
              The possible parametrizations P of the unitary transformation used
              in the DKH transformation supported by MOLCAS are:
              ||(P=O) Optimum parametrization (OPT);
              ||(P=E) Exponential parametrization (EXP);
              ||(P=S) Square-root parametrization (SQR);
              ||(P=M) McWeeny parametrization (MCW);
              ||(P=C) Cayley parametrization (CAY).
              Hence, the proper keyword for the 4th order relativistically corrected one-electron
              Hamiltonian and 3rd order relativistically corrected
              property integrals in the EXP parameterization would read as R04E03. If yy is larger than XX it is set to
              XX. If yy is omitted it will default to same value as XX. Recommended orders and parametrization is
              R02O02.
              Since the EXP parameterization employs a fast algorithm, it is
              recommended for high order DKH transformation.
              </HELP>
              </KEYWORD>

:kword:`RLOCal`
  Request local approximation to the relativistic exact decoupling approaches such
  as X2C, BSS and DKH. This option cannot be used together with point group symmetry.

  .. xmldoc:: <KEYWORD MODULE="SEWARD" NAME="RLOCAL" APPEAR="Relativistic local approximation" KIND="SINGLE" LEVEL="BASIC">
              %%Keyword: Rlocal <basic>
              <HELP>
              Request local approximation to the relativistic exact decoupling approaches such
              as X2C, BSS and DKH. This option cannot be used together with point group symmetry.
              </HELP>
              </KEYWORD>

  .. xmldoc:: </GROUP>

:kword:`Grid Input`
  Specification of numerical quadrature parameters, consult the numerical quadrature section of this
  manual.

Additional keywords for property calculations
:::::::::::::::::::::::::::::::::::::::::::::

.. class:: keywordlist

:kword:`VECTors`
  Requests a property calculation. For this purpose a file,
  :file:`INPORB`,
  must be available, which contains the MO's and occupation numbers
  of a wave function. A custom filename can be given with :kword:`FileOrb`.

  .. xmldoc:: <GROUP MODULE="SEWARD" NAME="PROPERTY" APPEAR="Property calculations options" KIND="BOX" LEVEL="BASIC">

  .. xmldoc:: <KEYWORD MODULE="SEWARD" NAME="VECT" APPEAR="Activate property calculation" KIND="SINGLE" LEVEL="BASIC">
              %%Keyword: Vectors <basic>
              <HELP>
              Requests a property calculation. For this purpose an INPORB file must be available.
              </HELP>
              </KEYWORD>

:kword:`FILEorb`
  The next line specifies the filename containing the orbitals and occupation numbers
  from which the properties will be computed. By default a file named :file:`INPORB`
  will be used.

  .. xmldoc:: <KEYWORD MODULE="SEWARD" NAME="FILE" APPEAR="Orbitals file" KIND="STRING" REQUIRE="VECT" LEVEL="BASIC">
              %%Keyword: FileOrb <basic>
              <HELP>
              The next line specifies the filename containing the orbitals and occupation numbers
              from which the properties will be computed. By default a file named INPORB
              will be used.
              </HELP>
              </KEYWORD>

:kword:`ORBCon`
  The keyword will force :program:`SEWARD`
  to produce a list of the orbital
  contributions to the properties being computed. The default is to
  generate a compact list.

  .. xmldoc:: <KEYWORD MODULE="SEWARD" NAME="ORBC" APPEAR="List Orbital contributions" KIND="SINGLE" REQUIRE="VECT" LEVEL="ADVANCED">
              %%Keyword: Orbcon <advanced>
              <HELP>
              The keyword will force SEWARD to produce a list of the orbital
              contributions to the properties being computed. The default is to
              generate a compact list.
              </HELP>
              </KEYWORD>

:kword:`THRS`
  The real entry on the following line specifies the threshold for
  the occupation number of an orbital in order for the
  :kword:`OrbCon`
  option to list the contribution of that orbital to a property.
  The default is 1.0d-6.

  .. xmldoc:: <KEYWORD MODULE="SEWARD" NAME="THRS" APPEAR="Property threshold" KIND="REAL" REQUIRE="ORBC" DEFAULT_VALUE="1.D-6" LEVEL="ADVANCED">
              <HELP>
              Specify the threshold for the occupation number of a orbital to allow that
              this orbital is included when listing orbital contributions to properties.
              </HELP>
              %%Keyword: Thrs <advanced>
              The real entry on the following line specifies the threshold for
              the occupation number of an orbital in order for the
              keyword "OrbCon" option to list the contribution of that orbital to a
              property. The default is 1.0d-6.
              </KEYWORD>

:kword:`ORBAll`
  When :kword:`OrbCon` is present, the keyword will force
  :program:`SEWARD` to produce a list of the computed properties
  of all orbitals (including unoccupied orbitals), and the
  properties are not weighted by orbital occupation numbers.
  The total electronic and nuclear
  contributions printed are the same as those printed by using
  :kword:`OrbCon` without :kword:`OrbAll`.

  .. xmldoc:: <KEYWORD MODULE="SEWARD" NAME="ORBA" APPEAR="List properties of all orbitals" KIND="SINGLE" REQUIRE="ORBC" LEVEL="ADVANCED">
              %%Keyword: OrbAll <advanced>
              <HELP>
              Print computed properties for all orbitals (including unoccupied orbitals),
              and not weighted by orbital occupation numbers.
              </HELP>
              Requires "OrbCon".
              </KEYWORD>

  .. xmldoc:: </GROUP>

Keywords for two-electron integrals
:::::::::::::::::::::::::::::::::::

.. class:: keywordlist

:kword:`NOPAck`
  The two-electron integrals will not be packed. The default is to
  pack the two-electron integrals.

  .. xmldoc:: <GROUP MODULE="SEWARD" NAME="2EL" APPEAR="2-electron integral options" KIND="BOX" LEVEL="BASIC">

  .. xmldoc:: <KEYWORD MODULE="SEWARD" NAME="NOPA" APPEAR="Unpacked 2-el. integrals" KIND="SINGLE" LEVEL="BASIC">
              %%Keyword: Nopack <advanced>
              <HELP>
              The two-electron integrals will not be packed. The default is to
              pack the two-electron integrals.
              </HELP>
              </KEYWORD>

:kword:`PKTHre`
  An entry specifies the desired accuracy for the packing
  algorithm, the default is 1.0d-14.

  .. xmldoc:: <KEYWORD MODULE="SEWARD" NAME="PKTH" APPEAR="Packing threshold" KIND="REAL" DEFAULT_VALUE="1.0D-14" LEVEL="BASIC">
              <HELP>
              Specifies the desired accuracy for the packing algorithm.
              </HELP>
              %%Keyword: Pkthre <advanced>
              An entry specifies the desired accuracy for the packing
              algorithm, the default is 1.0d-14.
              </KEYWORD>

:kword:`STDOut`
  Generate a two-electron integral file according to the standard of
  version 1 of |molcas|. The default is to generate the
  two-electron integrals according to the standard used since version 2 of
  |molcas|.

  .. xmldoc:: %%Keyword: Stdout <advanced>
              Generate a two-electron integral file according to the standard of
              version 1 of MOLCAS. The default is to generate the
              two-electron integrals according to the standard used since version 2 of
              MOLCAS.

:kword:`THREshold`
  Threshold for writing integrals to disk follows. The default is 1.0d-14.

  .. xmldoc:: <KEYWORD MODULE="SEWARD" NAME="THRE" APPEAR="To disk threshold" KIND="REAL" DEFAULT_VALUE="1.0D-14" LEVEL="BASIC">
              <HELP>
              Specify the threshold for writing integrals to disk.
              </HELP>
              %%Keyword: Threshold <advanced>
              Threshold for writing integrals to disk follows. The default is 1.0d-14.
              </KEYWORD>

:kword:`CUTOff`
  Threshold for ignoring the calculation of integrals based on the
  pair prefactor follows. The default is 1.0d-16.

  .. xmldoc:: <KEYWORD MODULE="SEWARD" NAME="CUTO" APPEAR="Integral threshold" KIND="REAL" DEFAULT_VALUE="1.0D-16" LEVEL="BASIC">
              <HELP>
              Specify the threshold for ignoring the calculation of integrals based on the pair prefactors.
              </HELP>
              %%Keyword: Cutoff <advanced>
              Threshold for ignoring the calculation of integrals based on the
              pair prefactor follows on the next line. The default is 1.0d-16.
              </KEYWORD>

  .. xmldoc:: </GROUP>

Keywords associated to electron-molecule scattering calculations within the framework of the :math:`R`\-matrix method
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

This section contains keyword which control the radial numerical integration
of the diffuse basis functions describing the scattered electrons in the
variational :math:`R`\-matrix approach. The activation of this option is controlled
by that the center of the diffuse basis is assigned the unique atom label DBAS.

.. class:: keywordlist

:kword:`RMAT`
  Radius of the :math:`R`\-matrix sphere (in bohr). This sphere is centered at the
  coordinate origin. The default is 10 bohr.

  .. xmldoc:: %%Keyword: RMAT <basic>
              Radius of the R-matrix sphere (in bohr). This sphere is centered
              at the coordinate origin. Default value is set to 10 bohr.

:kword:`RMEA`
  Absolute precision in radial integration.
  The default is 1d-9.

  .. xmldoc:: %%Keyword: RMEA <advanced>
              Absolute precision in radial integration for R-matrix calculation.
              The default is 1d-9.

:kword:`RMER`
  Relative precision in radial integration.
  The default is 1d-14.

  .. xmldoc:: %%Keyword: RMER <advanced>
              Relative precision in radial integration for R-matrix calculation.
              The default is 1d-14.

:kword:`RMQC`
  Effective charge of the target molecule.
  This is the effective charge seen by the incident electron outside
  of the :math:`R`\-matrix sphere. The default is 0d0.

  .. xmldoc:: %%Keyword: RMQC <advanced>
              Effective charge of the target molecule for R-matrix calculation.
              This is the effective charge seen by the incident electron outside
              of the R-matrix sphere. The default is 0d0.

:kword:`RMDI`
  Effective dipole of the target molecule.
  This is the effective dipole seen by the incident electron outside
  of the :math:`R`\-matrix sphere. The default is (0d0,0d0,0d0).

  .. xmldoc:: %%Keyword: RMDI <advanced>
              Effective dipole of the target molecule for R-matrix calculation.
              This is the effective dipole seen by the incident electron outside
              of the R-matrix sphere. The default is (0d0,0d0,0d0).

:kword:`RMEQ`
  Minimal value of the effective charge of the target molecule
  to be considered. This is also the
  minimal value of the components of the effective dipole
  to be considered. Default is 1d-8

  .. xmldoc:: %%Keyword: RMEQ <advanced>
              Minimal value of the effective charge of the target molecule
              to be considered in R-matrix calculation. This is also the
              minimal value of the components of the effective dipole
              to be considered in R-matrix calculation. Default is 1d-8.

:kword:`RMBP`
  Parameter used for test purposes in the definition of the
  Bloch term. Default is 0d0.

  .. xmldoc:: %%Keyword: RMBP <advanced>
              Parameter used for test purposes in the definition of the
              Bloch term in R-matrix calculation. Default is 0d0.

:kword:`CELL`
  Defines the three vectors of the unit cell (:math:`\vec{e}_1`, :math:`\vec{e}_2`, :math:`\vec{e}_3`).
  The optional keyword
  *Angstrom* before the definition of vectors would read data in Å.
  Must consist of three entries (four in the case of Å) which correspond to coordinates of the vectors. All the atoms which
  are defined after that key are considered as the atoms of the cell.

  .. xmldoc:: %%Keyword: CELL <advanced>
              Defines the three vectors of the unit cell (e1,e2,e3).
              The optional keyword
              Angstrom before the definition of vectors would read data in A.
              Must consist of three entries (four in the case of A) which correspond to coordinates of the vectors.
              All the atoms which
              are defined after that key are considered as the atoms of the cell.

:kword:`SPREad`
  Three integer numbers :math:`n_1`, :math:`n_2`, :math:`n_3` which define the spread of the unit cell along the unit cell vectors.
  For example, ``0 0 2`` would add all cell's atoms translated on :math:`-2\vec{e}_3`, :math:`-\vec{e}_3`, :math:`\vec{e}_3`, :math:`2\vec{e}_3`.
  This key must be placed **before** the definition of the unit cell atoms.

  .. xmldoc:: %%Keyword: SPREAD <advanced>
              Three integer numbers n_1, n_2, n_3 which define the spread of the unit cell along the unit cell vectors.
              For example, 0 0 2 would add all cell's atoms translated on -2*e3, -e3, e3, 2*e3.
              This key must be placed before the definition of the unit cell atoms.

Below follows an input for the calculation of integrals of a carbon
atom. The comments in the input gives a brief explanation of the
subsequent keywords.

.. extractfile:: ug/SEWARD.C.input

  &SEWARD
  Title= This is a test deck!
  * Remove integrals from a specific irreps
  Skip= 0 0 0 0 1 1 1 1
  * Requesting only overlap integrals.
  Multipole= 0
  * Request integrals for diamagnetic shielding
  DSHD= 0.0 0.0 0.0; 1; 0.0 0.0 0.0
  * Specify a title card
  * Request only one-electron integrals to be computed
  OneOnly
  * Specify group generators
  Symmetry= X Y Z
  * Enable an inline basis set
  Expert
  * Specify basis sets
  Basis set
  C.ANO-L...6s5p3d2f.
  Contaminant d
  C  0.0 0.0 0.0
  End of basis

The basis set label and the all electron basis set library
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

.. compound::

  The label, which defines the basis set for a given atom or set of atoms,
  is given as input after the
  keyword :kword:`Basis set`. It has the following general
  structure (notice that the last
  character is a period): ::

    atom.type.author.primitive.contracted.aux.

  where the different identifiers have the following meaning:

.. class:: whatnotlist

``atom``
  Specification of the atom by its chemical symbol.

``type``
  Gives the type of basis set (ANO, STO, ECP, etc.) according to
  specifications given in the basis set library, *vide supra*.
  Observe that the upper cased character of the type label
  defines the file name in the basis directory.

``author``
  First author in the publication where that basis set appeared.

``primitive``
  Specification of the primitive set (e.g. 14s9p4d3f).

``contracted``
  Specification of the contracted set to be selected. Some basis
  sets allow only one type of contraction, others all types up
  to a maximum. The first basis functions for each angular
  momentum is then used.
  **Note**, for the basis set types ANO and ECP, on-the-fly decontraction of the most
  diffuse functions are performed in case the number of contracted functions specified in this field
  exceeds what formally is specified in the library file.

``aux``
  Specification of the type of AIMP, for instance, to choose between
  non-relativistic and relativistic core AIMP's.

Only the identifiers ``atom``, ``type``,
and ``contracted`` have to be included in
the label. The others can be left out. However, the periods have to be
kept. Example --- the basis set label
"``C.ano-s...4s3p2d.``"
will by |molcas| be
interpreted as
"``C.ano-s.Pierloot.10s6p3d.4s3p2d.``",
which is the first
basis set in the ANO-S file in the
basis directory that fulfills the specifications given.

More information about basis set format can be found in
the section Advanced examples.

.. index::
   single: Numerical Integration

.. _UG\:sec\:dft:

Numerical integration
---------------------

Various Density Functional Theory (DFT) models can be used in |molcas|.
Energies and analytical gradients are available for all DFT models.
In DFT the exact exchange present in HF theory is replaced by a more general
expression, the exchange-correlation functional, which accounts for both the exchange energy, :math:`E_{\text{X}} [P]`
and the electron correlation energy, :math:`E_{\text{C}} [P]`.

.. _UG\:sec\:df_description:

Description
...........

We shall now describe briefly how the exchange and correlation energy terms look like.
The functionals used in DFT are integrals of some function of the electron density and optionally the gradient
of the electron density

.. math:: E_{\text{X}}[P] = \int f(\rho_{\alpha}(r), \rho_{\beta}(r), \nabla \rho_{\alpha}(r), \nabla \rho_{\beta}(r))\;dr

The various DFT methods differ in which function, :math:`f`, is used for :math:`E_{\text{X}}[P]` and for :math:`E_{\text{C}}[P]`.
In |molcas| pure DFT methods are supported, together with hybrid methods, in which the exchange functional
is a linear combination of the HF exchange and a functional integral of the above form.
The latter are evaluated by numerical quadrature.
In the :program:`SEWARD` input the parameters for the numerical integration can be set up.
In the :program:`SCF` and :program:`RASSCF` inputs the keywords for
using different functionals can be specified.
Names for the various pure DFT models are given by combining the names for the exchange and correlation functionals.

The DFT gradients has been implemented for both the fixed and the moving grid approach :cite:`Johnson:93,Handy:93,Baker:94`.
The latter is known to be translationally invariant by definition and is recommended in geometry optimizations.

.. .. _UG\:sec\:dft_files:

   Files
   .....

   .. class:: filelist

   :file:`RUNFILE`
     The run file will contain the parameters defining and controlling the numerical integration.

Input
.....

Below follows a description of the input to the numerical integration utility in the
:program:`SEWARD` input.

Compulsory keywords

.. class:: keywordlist

:kword:`GRID Input`
  This marks the beginning of the input to the numerical integration utility.

  .. xmldoc:: <GROUP MODULE="SEWARD" NAME="GRIDBLOCK" KIND="BLOCK" APPEAR="Numerical Quadrature Options" LEVEL="ADVANCED">
              %%Keyword: Grid input <basic>
              <HELP>
              Specification of numerical quadrature parameters.
              </HELP>

:kword:`END Of Grid-Input`
  This marks the end of the input to the numerical integration utility.

Optional keywords

.. class:: keywordlist

:kword:`GRID`
  It specifies the quadrature quality.
  The possible indexes that can follow are
  COARSE, SG1GRID, FINE, ULTRAFINE
  following the Gaussian98 convention.
  Default is FINE.

  .. xmldoc:: <KEYWORD MODULE="SEWARD" NAME="GRID" APPEAR="Generic Grids" KIND="CHOICE" LEVEL="ADVANCED" LIST="FINE,COARSE,SG1GRID,ULTRAFINE">
              %%Keyword: Grid (NQ) <advanced>
              <HELP>
              It specifies the quadrature quality.
              The possible indexes that can follow are
              COARSE, SG1GRID, FINE, ULTRAFINE
              following the Gaussian98 convention.
              Default is FINE.
              </HELP>
              </KEYWORD>

  .. :kword:`GLOBal`
       It activates the use of the global partitioning technique.

       .. .. xmldoc:: %%Keyword: Global <advanced>
                      It activates the use of the global partitioning technique.

:kword:`RQUAd`
  It specifies the radial quadrature scheme.
  Options are LOG3 (Mura and Knowles) :cite:`Mura:96`, BECKE
  (Becke) :cite:`BeckeG:88`, MHL (Murray et al.) :cite:`Murray:93`, TA (Treutler and
  Ahlrichs, defined for :math:`\ce{H}`\--\ :math:`\ce{Kr}`) :cite:`Treutler:95`, and LMG (Lindh et
  al.) :cite:`LMG:01`, respectively. The default is MHL.

  .. xmldoc:: <KEYWORD MODULE="SEWARD" NAME="RQUA" APPEAR="Radial grid type" KIND="CHOICE" LEVEL="ADVANCED" LIST="MHL,LOG3,BECKE,TA,LMG">
              %%Keyword: RQuad (NQ) <advanced>
              <HELP>
              It specifies the radial quadrature scheme.
              Options are LOG3 (Mura and Knowles), BECKE (Becke) , MHL (Murray et a.), TA (Treutler and
              Ahlrichs, defined for H-Kr), and LMG (Lindh et al.), respectively. The default is MHL.
              </HELP>
              </KEYWORD>

  .. xmldoc:: <SELECT MODULE="SEWARD" NAME="AGRID" APPEAR="Angular Grid" CONTAINS="LEBEDEV,LOBATTO,GGL" LEVEL="ADVANCED">
              <HELP>
              Specifies the type of angular grid. Options are Gauss-Gauss-Legendre (GGL), Lobatto, and
              Lebedev. Lebedev is the default.
              </HELP>

  .. xmldoc:: <KEYWORD MODULE="SEWARD" NAME="LEBEDEV" APPEAR="Lebedev angular grid" KIND="SINGLE" LEVEL="ADVANCED" EXCLUSIVE="LOBATTO,GGL" />

  .. xmldoc:: <KEYWORD MODULE="SEWARD" NAME="LOBATTO" APPEAR="Lobatto angular grid" KIND="SINGLE" LEVEL="ADVANCED" EXCLUSIVE="LEBEDEV,GGL" />

  .. xmldoc:: <KEYWORD MODULE="SEWARD" NAME="GGL" APPEAR="GGL angular grid" KIND="SINGLE" LEVEL="ADVANCED" EXCLUSIVE="LEBEDEV,LOBATTO" />

  .. xmldoc:: </SELECT>

:kword:`GGL`
  It activates the use of Gauss and Gauss--Legendre angular quadrature.
  Default is to use the Lebedev angular grid.

  .. xmldoc:: %%Keyword: GGL (NQ) <advanced>
              It activates the use of Gauss and Gauss-Legendre angular quadrature.
              Default is to use the Lebedev angular grid.

:kword:`LEBEdev`
  It turns on the Lebedev angular grid.

  .. xmldoc:: %%Keyword: Lebedev (NQ) <advanced>
              It turns on the Lebedev angular grid.

:kword:`LOBAtto`
  It activates the use of Lobatto angular quadrature.
  Default is to use the Lebedev angular grid.

  .. xmldoc:: %%Keyword: Lobatto (NQ) <advanced>
              It activates the use of Lobatto angular quadrature.
              Default is to use the Lebedev angular grid.

:kword:`LMAX`
  It specifies the angular grid size. Default is 29.

  .. xmldoc:: <KEYWORD MODULE="SEWARD" NAME="LMAX" APPEAR="Angular Grid Size" KIND="INT" DEFAULT_VALUE="29" LEVEL="ADVANCED">
              %%Keyword: LMax (NQ) <advanced>
              <HELP>
              Specifies the highest order of a real spherical harmonic which the angular grid will integrate exact. Default is 29.
              </HELP>
              </KEYWORD>

:kword:`NGRId`
  It specifies the maximum number of grid points to process at one instance.
  Default is 128 grid points.

  .. xmldoc:: <KEYWORD MODULE="SEWARD" NAME="NGRI" APPEAR="Batch size" KIND="INT" DEFAULT_VALUE="128" LEVEL="ADVANCED">
              %%Keyword: nGrid (NQ) <advanced>
              <HELP>
              It specifies the maximum number of grid points to process at one instance.
              Default is 128 grid points.
              </HELP>
              </KEYWORD>

:kword:`NOPRunning`
  It turns off the the angular prunning. Default is to prune.

  .. xmldoc:: <KEYWORD MODULE="SEWARD" NAME="NOPR" APPEAR="No Prunning of Grid" KIND="SINGLE" LEVEL="ADVANCED">
              %%Keyword: Noprunning (NQ) <advanced>
              <HELP>
              It turns off the the angular prunning. Default is to prune.
              </HELP>
              </KEYWORD>

:kword:`NR`
  It is followed by the number of radial grid points. Default is 75 radial grid points.

  .. xmldoc:: <KEYWORD MODULE="SEWARD" NAME="NR" APPEAR="Radial Grid Size" KIND="INT" DEFAULT_VALUE="75" LEVEL="ADVANCED">
              %%Keyword: nR (NQ) <advanced>
              <HELP>
              It is followed by the number of radial grid points. Default is 75 radial grid points.
              </HELP>
              </KEYWORD>

  .. :kword:`WHOLe`
       It activates the use of routines which scan the whole atomic grid for
       each sub block. Default is to only scan the relevant part.

  ..   .. xmldoc:: %%Keyword: Whole (NQ) <advanced>
                   It activates the use of routines which scan the whole atomic grid for
                   each sub block. Default is to only scan the relevant part.

:kword:`FIXEd grid`
  Use a fixed grid in the evaluation of the gradient. This corresponds to
  using the grid to numerically evaluate the analytic gradient expression.
  Default is to use a moving grid.

  .. xmldoc:: <SELECT MODULE="SEWARD" NAME="GRID_TYPE" APPEAR="Grid type" CONTAINS="MOVING,FIXED" LEVEL="ADVANCED">
              <HELP>
              Specify if the grid should be fixed or moving. A moving grid is default.
              </HELP>

  .. xmldoc:: <KEYWORD MODULE="SEWARD" NAME="FIXE" APPEAR="Fixed grid" KIND="SINGLE" LEVEL="ADVANCED">
              %%Keyword: Fixed Grid (NQ) <advanced>
              <HELP>
              Use a fixed grid in the evaluation of the gradient. This corresponds to
              using the grid to numerically evaluate the analytic gradient expression.
              Default is to use a moving grid.
              </HELP>
              </KEYWORD>

:kword:`MOVIng grid`
  Use a moving grid in the evaluation of the gradient. This correspond to
  evaluating the gradient of the numerical expression of the DFT energy. This is the default.

  .. xmldoc:: <KEYWORD MODULE="SEWARD" NAME="MOVI" APPEAR="Moving grid" KIND="SINGLE" LEVEL="ADVANCED">
              %%Keyword: Moving grid (NQ) <advanced>
              <HELP>
              Use a moving grid in the evaluation of the gradient. This correspond to
              evaluating the gradient of the numerical expression of the DFT energy. This is the default.
              </HELP>
              </KEYWORD>

  .. xmldoc:: </SELECT>

:kword:`RTHReshold`
  It is followed by the value for the the radial threshold.
  Default value is 1.0D-13.

  .. xmldoc:: <KEYWORD MODULE="SEWARD" NAME="RTHR" APPEAR="Radial Grid Threshold" KIND="REAL" DEFAULT_VALUE="1.0D-13" LEVEL="ADVANCED">
              %%Keyword: RThreshold (NQ) <advanced>
              <HELP>
              Follows the value for the the radial threshold.
              Default value is 1.0D-13.
              </HELP>
              </KEYWORD>

:kword:`T_X`
  Threshold for screening in the assembling of the density on the grid.
  Default value is 1.0D-18.

  .. xmldoc:: <KEYWORD MODULE="SEWARD" NAME="T_X" KIND="REAL" DEFAULT_VALUE="1.0D-18" LEVEL="ADVANCED">
              %%Keyword: T_X (NQ) <advanced>
              <HELP>
              Threshold for screening in the assembling of the density on the grid.
              Default value is 1.0D-18.
              </HELP>
              </KEYWORD>

:kword:`T_Y`
  Threshold for screening in the assembling of the integrals.
  Default value is 1.0D-11.

  .. xmldoc:: <KEYWORD MODULE="SEWARD" NAME="T_Y" KIND="REAL" DEFAULT_VALUE="1.0D-11" LEVEL="ADVANCED">
              %%Keyword: T_Y (NQ) <advanced>
              <HELP>
              Threshold for screening in the assembling of the integrals.
              Default value is 1.0D-11.
              </HELP>
              </KEYWORD>

:kword:`NOSCreening`
  Turn off any screening in the numerical integration.

  .. xmldoc:: <KEYWORD MODULE="SEWARD" NAME="NOSC" APPEAR="NoScreening" KIND="SINGLE" LEVEL="ADVANCED">
              %%Keyword: NOSCreening (NQ) <basic>
              <HELP>
              Turn off any screening in the numerical integration.
              </HELP>
              </KEYWORD>

:kword:`CROWding`
  The crowding factor, according to MHL, used in the pruning of the angular
  grid close to the nuclei. Default value 3.0.

  .. xmldoc:: <KEYWORD MODULE="SEWARD" NAME="CROW" APPEAR="CROWDING" KIND="REAL" DEFAULT_VALUE="3.0" LEVEL="ADVANCED">
              %%Keyword: CROWding (NQ) <basic>
              <HELP>
              The crowding factor, according to MHL, used in the pruning of the angular
              grid close to the nuclei. Default value 3.0.
              </HELP>
              </KEYWORD>

  .. xmldoc:: <KEYWORD MODULE="SEWARD" NAME="NORO" KIND="SINGLE" LEVEL="UNDOCUMENTED" />

  .. xmldoc:: </GROUP>

The :program:`SCF` and :program:`RASSCF` programs have
their own keywords to decide which functionals to use in a DFT calculation.

Below follows an example of a DFT calculation with two different functionals.

.. extractfile:: ug/DFT.H.input

  &GATEWAY
  Basis set
  H.3-21G.....
  H1 0.0  0.0 0.0
  End of basis

  &SEWARD
  Grid input
   RQuad= Log3; nGrid= 50000; GGL; lMax= 26; Global
  End of Grid Input

  &SCF; Occupations=1; KSDFT=LDA5; Iterations= 1 1

  &SCF; Occupations= 1; KSDFT=B3LYP; Iterations= 1 1

Relativistic operators
----------------------

The current different implementation of all relativistic operators in |molcas| as
described in the following sections has been programmed and tested in
Ref. :cite:`R2C_review`

Using the Douglas--Kroll--Hess Hamiltonian
..........................................

For all-electron calculations, the preferred way is to use the scalar-relativistic
Douglas--Kroll--Hess (DKH) Hamiltonian, which, in principle, is available up to arbitrary
order in Molcas; for actual calculations, however,
the standard 2nd order is usually fine, but one may use a higher order than 8th order
by default to be on the safe side.

.. compound::

  The arbitrary-order Hamiltonian is activated by setting ::

    RXXPyy

  somewhere in the :program:`SEWARD` input, where
  the ``XX`` denotes the order of the DKH Hamiltonian in the external potential.
  I.e., for the standard 2nd-order Hamiltonian you may use ``R02O``.
  Note
  in particular that the parametrization ``P`` does not affect the Hamiltonian up to
  fourth order. Therefore, as long as you run calculations with DKH Hamiltonians below
  5th order you may use any symbol for the parametrization as they would all yield the
  same results.

The possible parametrizations ``P`` of the unitary transformation used
in the DKH transformation supported by |molcas| are:

.. container:: list

  ``P=O``: Optimum parametrization (OPT)

  ``P=E``: Exponential parametrization (EXP)

  ``P=S``: Square-root parametrization (SQR)

  ``P=M``: McWeeny parametrization (MCW)

  ``P=C``: Cayley parametrization (CAY)

Up to fourth order (``XX``\=04) the DKH Hamiltonian is independent of the chosen
parametrization.
Higher-order DKH Hamiltonians depend slightly on the chosen parametrization of the unitary
transformations applied in order to decouple the Dirac Hamiltonian.
Since the EXP parameterization employs a fast algorithm :cite:`DKH_polynomial`, it is recommended
for high-order DKH transformation.

For details on the arbitrary-order DKH Hamiltonians see :cite:`DKH_Theory` with respect to theory,
:cite:`DKH_Implementation` with respect to aspects of implementation, and
:cite:`DKH_Principles` with respect to general principles of DKH.
The current version of |molcas| employs different algorithms,
but the polynomial cost scheme of the DKH implementation as described in :cite:`DKH_polynomial` is
used as the default algorithm. The implementation in MOLCAS has been presented in :cite:`R2C_review`.

For details on the different parametrizations of the unitary transformations see :cite:`DKH_Parameterization`.

Douglas--Kroll--Hess transformed properties
...........................................

As mentioned above, four-component molecular property operators need to be DKH
transformed as well when going from a four-component to a two- or one-component description;
the results do not coincide with the well-known corresponding nonrelativistic expressions
for a given property but are properly picture change corrected.

.. compound::

  The transformation of electric-field-like molecular property operators can be carried
  out for any order smaller or equal to the order chosen for the scalar-relativistic DKH Hamiltonian.
  In order to change the default transformation of order 2, you may concatenate the
  input for the DKH Hamiltonian by 2 more numbers specifying the order in the property, ::

    RxxPyy

  where ``yy`` denotes the order of the Hamiltonian starting with first order 01.
  The DKH transformation is then done automatically for all one-electron electric-field-like
  one-electron property matrices.

Also note that the current implementation of both the Hamiltonian and the property
operators is carried out in the full, completely decontracted basis set of the
molecule under consideration. The local nature of the relativistic contributions is not
yet exploited and hence large molecules may require considerable computing time for
all higher-order DKH transformations.

For details on the arbitrary-order DKH properties see :cite:`DKH_Theory2` with respect to theory
and :cite:`DKH_Implementation2,R2C_review` with respect to implementation aspects.

Using the X2C/Barysz--Sadlej--Snijders Hamiltonian
..................................................

Exact decoupling of the relativistic Dirac Hamiltonian can be achieved with infinite-order
approaches, such as the so-called X2C (exact-two-component) and BSS (Barysz--Sadlej--Snijders)
method. In Molcas, both methods are available for all-electron calculations. The
evaluation of transformation matrices employs a non-iterative scheme.

.. which is superior to the iterative scheme in convergence behavior and computational cost.

The exact decoupling Hamiltonian is activated by setting either ``RX2C`` or ``RBSS``
somewhere in the :program:`SEWARD` input, where ``RX2C`` and ``RBSS`` denote using the
scalar (one-component) X2C and BSS Hamiltonian respectively. The one-electron Hamiltonian
as well as the property integrals will be transformed according to the given exact decoupling method.
In other words, all property integrals are by default picture change corrected.

The computation time of the X2C/BSS method is almost the same as of the DKH method at 8th order,
while X2C is a little bit faster than BSS since the additional free-particle Foldy--Wouthuysen
transformation is skipped in the X2C approach :cite:`R2C_review`. For molecules including only light atoms,
the DKH method with low orders (< 8) is enough to account for the relativistic effects.

.. As for heavy
   elements containing molecules, DKH with low orders is not accurate enough while DKH with high
   orders is computational expensive, the X2C method is then recommended in this case.

The differences between different exact decoupling relativistic methods are very small compared
to errors introduced by other approximations, such as the basis set incompleteness, approximate
density functionals, etc. Therefore, any exact decoupling model is acceptable for the
treatment of relativistic effects in molecular calculations.

For details on the exact decoupling approaches see :cite:`R2C_review` with respect to theories and
comparison of numerical results, :cite:`X2C_kut05jcp,X2C_liu06jcp,X2C_pen07jcp` for the X2C method, and
:cite:`BSS_bar97ijqc,BSS_ked07cpl` for the BSS method.

Local approximation to relativistic decoupling
..............................................

The computational cost for relativistic transformations increases rapidly if the molecule becomes
larger. The local DLU scheme :cite:`DLU_Theory` was proposed to reduce the computational cost
based on the atomic decomposition of the (unitary) decoupling transformation. It is important to note that
the DLU scheme can be applied to any kind of relativistic approaches mentioned above (i.e., DKH, BSS, and X2C). It was found :cite:`DLU_Theory`
that the DLU approximation introduces very small errors for total energies, which are less than
0.01 milli-hartree for molecules including heavy atoms.
The local approximation is activated by setting :kword:`RLOCal` somewhere in the :program:`SEWARD` input.

The direct local approximation to the Hamiltonian, called DLH in Ref. :cite:`DLU_Theory` may be activated by setting ``RLOCal=DLH``.
However, as DLH is not superior to the DLU scheme, but introduces slightly larger errors, it is not recommended.

The picture-change effect for molecular properties is automatically taken care of when a local approximation is employed for the transformed
operator. The default order of DKH transformed properties is set to the same as the order of the DKH Hamiltonian.

It is important to note the local approximation cannot be used together with point group symmetry
in the current implementation. Because the relativistic transformation is applied in the molecular
orbital (MO) representation instead of the atomic orbital (AO) representation. Thus, the program will
report an error and exit if symmetry is used.

.. xmldoc:: <KEYWORD MODULE="SEWARD" NAME="FOOC" KIND="SINGLE" LEVEL="UNDOCUMENTED" />

.. xmldoc:: <KEYWORD MODULE="SEWARD" NAME="UNCONTRACTED" KIND="SINGLE" LEVEL="UNDOCUMENTED" />

.. xmldoc:: <KEYWORD MODULE="SEWARD" NAME="NEMO" KIND="SINGLE" LEVEL="UNDOCUMENTED" />

.. xmldoc:: <KEYWORD MODULE="SEWARD" NAME="CLIGHT" KIND="REAL" LEVEL="UNDOCUMENTED" />

.. xmldoc:: <KEYWORD MODULE="SEWARD" NAME="PAMFI" KIND="INT" LEVEL="UNDOCUMENTED" />

.. xmldoc:: <KEYWORD MODULE="SEWARD" NAME="DOANALYTICAL" KIND="SINGLE" LEVEL="UNDOCUMENTED" />

.. xmldoc:: <KEYWORD MODULE="SEWARD" NAME="PRINT" KIND="INTS_COMPUTED" SIZE="2" LEVEL="UNDOCUMENTED" />

.. xmldoc:: <GROUP MODULE="SEWARD" NAME="CHOINPUT" KIND="BLOCK" LEVEL="UNDOCUMENTED">
            <KEYWORD MODULE="SEWARD" NAME="THRC" KIND="REAL" LEVEL="UNDOCUMENTED" />
            <KEYWORD MODULE="SEWARD" NAME="PRIN" KIND="INT" LEVEL="UNDOCUMENTED" />
            <KEYWORD MODULE="SEWARD" NAME="BUFF" KIND="INT" LEVEL="UNDOCUMENTED" />
            <KEYWORD MODULE="SEWARD" NAME="THRD" KIND="REAL" LEVEL="UNDOCUMENTED" />
            <KEYWORD MODULE="SEWARD" NAME="DMP1" KIND="REAL" LEVEL="UNDOCUMENTED" />
            <KEYWORD MODULE="SEWARD" NAME="DMP2" KIND="REAL" LEVEL="UNDOCUMENTED" />
            <KEYWORD MODULE="SEWARD" NAME="SPAN" KIND="REAL" LEVEL="UNDOCUMENTED" />
            <KEYWORD MODULE="SEWARD" NAME="MINQ" KIND="INT" LEVEL="UNDOCUMENTED" />
            <KEYWORD MODULE="SEWARD" NAME="MAXQ" KIND="INT" LEVEL="UNDOCUMENTED" />
            <KEYWORD MODULE="SEWARD" NAME="SCRE" KIND="SINGLE" LEVEL="UNDOCUMENTED" />
            <KEYWORD MODULE="SEWARD" NAME="QUAL" KIND="INT" LEVEL="UNDOCUMENTED" />
            <KEYWORD MODULE="SEWARD" NAME="THRN" KIND="REAL" LEVEL="UNDOCUMENTED" />
            <KEYWORD MODULE="SEWARD" NAME="WARN" KIND="REAL" LEVEL="UNDOCUMENTED" />
            <KEYWORD MODULE="SEWARD" NAME="TOON" KIND="REAL" LEVEL="UNDOCUMENTED" />
            <KEYWORD MODULE="SEWARD" NAME="NOAB" KIND="SINGLE" LEVEL="UNDOCUMENTED" />
            <KEYWORD MODULE="SEWARD" NAME="IOVE" KIND="INT" LEVEL="UNDOCUMENTED" />
            <KEYWORD MODULE="SEWARD" NAME="FRAC" KIND="INTS" SIZE="2" LEVEL="UNDOCUMENTED" />
            <KEYWORD MODULE="SEWARD" NAME="PARA" KIND="SINGLE" LEVEL="UNDOCUMENTED" />
            <KEYWORD MODULE="SEWARD" NAME="1-CE" KIND="SINGLE" LEVEL="UNDOCUMENTED" />
            <KEYWORD MODULE="SEWARD" NAME="ONES" KIND="SINGLE" LEVEL="UNDOCUMENTED" />
            <KEYWORD MODULE="SEWARD" NAME="ADDR" KIND="INT" LEVEL="UNDOCUMENTED" />
            <KEYWORD MODULE="SEWARD" NAME="RSTD" KIND="SINGLE" LEVEL="UNDOCUMENTED" />
            <KEYWORD MODULE="SEWARD" NAME="RSTC" KIND="SINGLE" LEVEL="UNDOCUMENTED" />
            <KEYWORD MODULE="SEWARD" NAME="REOR" KIND="SINGLE" LEVEL="UNDOCUMENTED" />
            <KEYWORD MODULE="SEWARD" NAME="CHOM" KIND="INT" LEVEL="UNDOCUMENTED" />
            </GROUP>

.. xmldoc:: <INCLUDE MODULE="GATEWAY" />

.. xmldoc:: </MODULE>
