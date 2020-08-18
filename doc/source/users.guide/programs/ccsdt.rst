.. index::
   single: Program; CCSDT
   single: CCSDT

.. _sec\:ccsdt:

:program:`ccsdt`
================

.. only:: html

  .. contents::
     :local:
     :backlinks: none

.. xmldoc:: <MODULE NAME="CCSDT">
            %%Description:
            <HELP>
            The CCSDT set of programs performs the iterative ROHF CCSD
            procedure, optionally followed by the (T) calculation contribution.
            It requires the JOBIPH file produced by RASSCF, and TRAONE and TRAINT
            files produced by MOTRA.
            </HELP>

:program:`CCSDT` performs the iterative single determinant CCSD procedure for
open shell systems and the noniterative triple contribution calculation to
the CCSD energy.
For further details the reader is referred to
:numref:`Sections %s <TUT:sec:ccsdt>` and
:numref:`%s <TUT:sec:rp_wf>` of the tutorials and examples manual.

.. index::
   pair: Dependencies; CCSDT

.. _sec\:ccsdt_dependencies:

Dependencies
------------

:program:`CCSDT` requires a previous run of the :program:`RASSCF` program
to produce orbital energies, Fock matrix elements, wave function
specification, and some other parameters stored in file :file:`JOBIPH`.
The :program:`RASSCF` program should be run with the options that produce
canonical output orbitals, which is not default.
:program:`CCSDT` also requires transformed integrals produced by :program:`MOTRA`
and stored in the files :file:`TRAONE` and :file:`TRAINT`.

It is well known that the CCSD procedure brings the spin
contamination into the final
wave function :math:`\ket{\Psi}` even in the case where the reference function
:math:`\ket{\Phi}` is the proper
spin eigenfunction. The way how to reduce the spin
contamination and mainly the number of independent amplitudes is to introduce
the spin adaptation.
Besides the standard nonadapted (spinorbital) CCSD procedure this program
allows to use different levels of spin
adaptation of CCSD amplitudes (the recommended citations are Refs.
:cite:`ccsd_neo2,ccsd_neo1`):

* DDVV T2 adaptation.

  This is the most simple and most universal scheme, in which only the dominant
  part of T2 amplitudes, namely those where both electrons are excited from
  *doubly occupied (inactive)* to *virtual (secondary)* orbitals, are adapted.
  The remaining types of amplitudes are left unadapted, i.e. in the spinorbital form.
  This alternative is an excellent approximation to the full adaptation and
  can be used for any multiplet.

* Full T1 and T2 adaptation (only for doublet states yet).

  In this case full spin adaptation of all types of amplitudes is performed.
  In the present implementation this version is limited to systems with
  the single unpaired electrons, i.e. to the doublet states only.

Besides these two possibilities there are also available some
additional partial ones (see keyword
:kword:`ADAPTATION` in :numref:`sec:ccsdt_input`). These adaptations are
suitable only for some specific purposes. More details on spin adaptation in
the CCSD step can be found in Refs. :cite:`ccsd_neo1,ccsd_neo2,ccsd_kno`.
The current implementation of the spin adaptation saves no computer time. A more
efficient version is under development.

The noniterative triples calculation can follow these approaches:

* CCSD + T(CCSD) --- according to Urban et al. :cite:`t3_urban`
* CCSD(T) --- according to Raghavachari et al. :cite:`t3_ragh`
* CCSD(T) --- according e.g. to Watts et al. :cite:`t3_watts`

Actual implementation and careful analysis and discussion of these
methods is described in Ref. :cite:`t3_neo`, which is a recommended reference
for this program.

The first alternative represents the simplest noniterative T3 treatment and contains
only pure :math:`\braket{T3}{W T2}` term. Second possibility represents the well known
extension to the first one by the :math:`\braket{T3}{W T1}` term
(:math:`W` is the two electron perturbation). For closed shell
systems this is the most popular and most frequently used noniterative triples
method. For single determinant open shell systems, described by the
ROHF reference
function standard (Raghavachari et. al.) method needs to be extended by the
additional fourth order energy term, namely
:math:`\braket{T3}{U T2}` (:math:`U` is the off-diagonal part of the Fock operator).

In contrast to the iterative CCSD procedure, noniterative approaches are not
invariant with respect to the partitioning of the Hamiltonian.
Hence, we obtain
different results using orbital energies, Fock matrix elements
or some other quantities in the
denominator. According to our experiences :cite:`t3_neo`,
diagonal Fock matrix elements in the
denominator represent the best choice. Using of other alternatives
requires some experience.
Since the triple excitation contribution procedure works strictly within the restricted formalism, resulting
noniterative triples contributions depend also on the choice of the reference
function. However, differences between this approach (with the reference
function produced by a single determinant RASSCF procedure and the diagonal
Fock matrix elements considered in the denominator) and the corresponding
invariant treatment (with the semicanonical orbitals)
are found to be chemically negligible.

For noniterative T3 contribution both non-adapted (spin-orbital) and spin-adapted
CCSD amplitudes can be used. For more details, see Ref. :cite:`t3_neo`.

.. index::
   pair: Files; CCSDT

.. _sec\:ccsdt_files:

Files
-----

Input files
...........

:program:`CCSDT` will use the following input
files: :file:`TRAONE`, :file:`TRAINT`, :file:`RUNFILE`, :file:`JOBIPH`,
(for more information see :numref:`UG:sec:files_list`).

Output files
............

.. class:: filelist

:file:`RSTART`
  file with CC amplitudes and CC energy.
  The name of the file can be changed using keyword :kword:`RESTART`.
  It contains restart information, like
  T1aa, T1bb, T2aaaa, T2bbbb, T2abab, CC energy and the number of iterations.

:file:`T3hfxyy`
  These files contain integrals of :math:`\braket{ia}{bc}` type where *x*
  represents
  the symmetry and *yy* the value of the given index :math:`i`.
  The number of
  these files is equal to the number of :math:`\alpha` occupied orbitals
  (*inactive + active*).

.. index::
   pair: Input; CCSDT

.. _sec\:ccsdt_input:

Input
-----

The input for each module is preceded by its name like: ::

  &CCSDT

.. class:: keywordlist

:kword:`TITLe`
  This keyword should be followed by precisely one title line.
  It should not begin with a blank (else it will not be printed!)
  This keyword is *optional*.

  .. xmldoc:: <KEYWORD MODULE="CCSDT" NAME="TITLE" APPEAR="Title" KIND="STRING" LEVEL="BASIC">
              %%Keyword: TITLe <basic>
              <HELP>
              Followed by precisely one title line, not beginning with a blank.
              </HELP>
              </KEYWORD>

:kword:`CCSD`
  This keyword specifies that only CCSD calculation will follow
  and the integrals will be prepared for the CCSD procedure only.
  This keyword is *optional*. (Default=OFF)

  .. xmldoc:: <SELECT MODULE="CCSDT" NAME="ANYTRIP" APPEAR="Any triples?" CONTAINS="CCSD,CCT" LEVEL="BASIC">

  .. xmldoc:: <KEYWORD MODULE="CCSDT" NAME="CCSD" APPEAR="CCSD only" KIND="SINGLE" LEVEL="BASIC" EXCLUSIVE="CCT">
              %%Keyword: CCSD <basic>
              <HELP>
              Specifies that only CCSD calculation will follow.
              </HELP>
              </KEYWORD>

:kword:`CCT`
  This keyword specifies that after CCSD calculation also noniterative
  T3 step will follow. For such calculations this key must
  be switched on. The integrals for the triple contribution calculation
  will then be prepared.
  This keyword is *optional*. (Default=ON)

  .. xmldoc:: <KEYWORD MODULE="CCSDT" NAME="CCT" APPEAR="Triples (default)" KIND="SINGLE" LEVEL="BASIC" EXCLUSIVE="CCSD">
              %%Keyword: CCT <basic>
              <HELP>
              Specifies that after CCSD, also the noniterative T3 calculation will follow.
              </HELP>
              </KEYWORD>

  .. xmldoc:: </SELECT>

:kword:`ADAPtation`
  The parameter on the following line defines the type of spin adaptations
  of CCSD amplitudes.

  .. container:: list

    0 --- no spin adaptation --- full spinorbital formalism

    1 --- T2 DDVV spin adaptation

    2 --- T2 DDVV + T1 DV spin adaptation (only recommended for specific purposes,
    since the adaptation of T1 included incompletely)

    3 --- full T2 and T1 spin adaptation (in current implementations
    limited to doublets only)

    4 --- full T2 adaptation without SDVS coupling (for doublets only)

  This keyword is *optional*. (Default=0)

  .. xmldoc:: <KEYWORD MODULE="CCSDT" NAME="ADAP" APPEAR="Adaptation" KIND="CHOICE" LIST="0: No adaptation,1: T2 DDVV,2: T2 DDVV and T1 DV,3: Full T2 and T1,4: Full T2" LEVEL="BASIC" DEFAULT_VALUE="0">
              <HELP>
              Choose how CCSD amplitudes should be spin adapted (if at all).
              </HELP>
              </KEYWORD>
              %%Keyword: ADAPtation <basic>
              Sets the type of CCSD amplitudes spin adaptation.
              ||0 - None
              ||1 - T2 DDVV
              ||2 - T2 DDVV + T1 DV
              ||3 - Full T2 and T1 spin adaptation (doublets only)
              ||4 - Full T2 adaptation without SDVS coupling (doublets only)

:kword:`DENOminators`
  The parameter on the following line specifies the type of denominators that
  will be used in the CCSD procedure.

  .. container:: list

    0 --- diagonal Fock matrix elements (different for :math:`\alpha` and :math:`\beta`
    spins)

    1 --- spin averaged diagonal Fock matrix elements ---
    :math:`\frac{f_{\alpha\alpha}+f_{\beta\beta}}{2}`

    2 --- orbital energies

  In some cases alternatives 1 and 2 are identical.
  For nonadapted CCSD calculations the resulting CCSD energy
  is invariant with respect to the selection of denominators.
  However, convergence may be affected.

  In the present implementation a symmetric denominators
  (i.e. the input 1 or 2) should be used for spin adapted CCSD calculations.
  This keyword is *optional*. (Default=0)

  .. xmldoc:: <KEYWORD MODULE="CCSDT" NAME="DENO" APPEAR="Denominators" KIND="CHOICE" LIST="0: Diagonal Fock matrix elements,1: Spin averaged diagonal,2: Orbital energies" LEVEL="BASIC" DEFAULT_VALUE="0">
              <HELP>
              Choose the type of denominators in the CCSD procedure.
              </HELP>
              </KEYWORD>
              %%Keyword: DENOminators <basic>
              Sets the type of denominators in the CCSD procedure.
              ||0 - Diagonal Fock matrix elements
              ||1 - Spin averaged diagonal Fock matrix elements
              ||2 - Orbital energies

:kword:`SHIFts`
  Following line contains *socc* and *svirt* levelshift values for occupied and
  virtual orbitals respectively. Typical values are in the range 0.0--0.5 (in *a.u.*) ::

    dp(occ)=dp(occ)-socc
    dp(virt)=dp(virt)+svirt

  For spin adaptations 3 and 4 only inactive (D) and active (V) orbitals
  will be shifted, due to the character of the adaptation scheme. For other cases all
  orbitals are shifted.

  This keyword is *optional*. (Defaults: *socc* = 0.0, *svirt* = 0.0)

  .. xmldoc:: <KEYWORD MODULE="CCSDT" NAME="SHIFT" APPEAR="Shifts" KIND="REALS" SIZE="2" LEVEL="ADVANCED" DEFAULT_VALUE="0.0">
              <HELP>
              Enter level shift values for occupied and virtual orbitals.
              </HELP>
              </KEYWORD>
              %%Keyword: SHIFts <advanced>
              On the following line the level shift values for occupied and virtual
              orbitals needs to be specified, typically around 0.0 - 0.5.

:kword:`TRIPles`
  The parameter on the following line specifies the
  type of noniterative triples
  procedure. There are three different types of perturbative triples available
  (see :numref:`sec:ccsdt`).

  .. container:: list

    0 --- CCSD approach (no triples step)

    1 --- CCSD+T(CCSD) according to Urban et. al :cite:`t3_urban`

    2 --- CCSD(T) according to Raghavachari et. al. :cite:`t3_ragh`

    3 --- CCSD(T) according e.g. to Watts et. al. :cite:`t3_watts`

  This keyword is *optional*. (Default=3)

  .. xmldoc:: <KEYWORD MODULE="CCSDT" NAME="TRIPLES" APPEAR="What triples" KIND="CHOICE" LIST="0: CCSD,1: CCSD + T(CCSD),2: CCSD(T) Raghavachari,3: CCSD(T) Watts" LEVEL="BASIC" DEFAULT_VALUE="3">
              <HELP>
              Choose the type of triples contribution calculation.
              </HELP>
              </KEYWORD>
              %%Keyword: TRIPles <basic>
              Sets the type of triples contribution calculation.
              ||0 - CCSD
              ||1 - CCSD + T(CCSD)   (Urban et al.)
              ||2 - CCSD(T)          (Raghavachari et al.)
              ||3 - CCSD(T)          (Watts et al.)

:kword:`T3DEnominators`
  The parameter on the following line specifies the type of denominators that
  will be used in noniterative triples procedure.

  .. container:: list

    0 --- diagonal Fock matrix elements (different for :math:`\alpha` and :math:`\beta`
    spins)

    1 --- spin averaged diagonal Fock matrix elements ---
    :math:`\frac{f_{\alpha\alpha}+f_{\beta\beta}}{2}`

    2 --- orbital energies

  In some cases alternatives 1 and 2 are identical.
  This keyword is *optional*. (Default=0)

  .. xmldoc:: <KEYWORD MODULE="CCSDT" NAME="T3DEN" APPEAR="T3 denominators" KIND="CHOICE" LIST="0: Diagonal,1: Spin averaged,2: Orbital energies" LEVEL="ADVANCED" DEFAULT_VALUE="0">
              <HELP>
              Choose the type of denominators used in the (T) calculation procedure.
              </HELP>
              </KEYWORD>
              %%Keyword: T3DEnominators <advanced>
              Sets the type of denominators used in the (T) calculation procedure.
              ||0 - Diagonal Fock matrix elements
              ||1 - Spin averaged Fock matrix elements
              ||2 - Orbital energies

:kword:`T3SHifts`
  The following line contains *socc* and *svirt* levelshift values for
  occupied and virtual orbitals respectively.
  Typical values are in the range 0.0--0.5 (in *a.u.*) ::

    dp(occ)=dp(occ)-socc
    dp(virt)=dp(virt)+svirt

  In contrast to the iterative CCSD procedure, in noniterative T3 step results are
  not invariant with respect to the denominator shifting. It is extremely dangerous
  to use any other than 0.0 0.0 shifts here, since resulting T3 energy may have
  no physical meaning. This keyword may be useful only in estimating some
  trends in resulting energy, however, using of default values is strongly
  recommended.

  This keyword is *optional*. (Defaults: *socc* = 0.0, *svirt* = 0.0)

  .. xmldoc:: <KEYWORD MODULE="CCSDT" NAME="T3SH" APPEAR="T3 Shifts" KIND="REALS" SIZE="2" LEVEL="ADVANCED">
              <HELP>
              Enter two numbers with level shifts for occupied and virtual orbitals
              in (T) calculations. Use with care, if at all, and consult the manual.
              </HELP>
              </KEYWORD>
              %%Keyword: T3SHifts <advanced>
              This keyword is followed by two numbers that set the levelshift values
              for occupied and virtual orbitals in (T) calculations. The default values
              (0,0) should not normally be changed.

:kword:`ITERations`
  This keyword is followed on the next line by the maximum number
  of iterations in the CCSD procedure. In the case of the RESTART run this is the
  number of last allowed iteration, since counting of iterations in
  RESTART run starts from the value taken from the :file:`RSTART` file.
  This keyword is *optional*. (Default=30)

  .. xmldoc:: <KEYWORD MODULE="CCSDT" NAME="ITER" APPEAR="MAX iter" KIND="INT" LEVEL="BASIC" MIN_VALUE="0" DEFAULT_VALUE="30">
              %%Keyword: ITERations <basic>
              <HELP>
              Sets the maximum number of CCSD iterations (Default:30).
              </HELP>
              </KEYWORD>

:kword:`ACCUracy`
  The real value on the following line defines the convergence criterion on
  CCSD energy. This keyword is *optional*. (Default=1.0d-7)

  .. xmldoc:: <KEYWORD MODULE="CCSDT" NAME="ACCU" APPEAR="Accuracy" KIND="REAL" LEVEL="BASIC" MIN_VALUE="0.0" DEFAULT_VALUE="1.0d-7">
              <HELP>
              Change the default convergence criterion (1.0D-7) on CCSD energy.
              </HELP>
              </KEYWORD>
              %%Keyword: ACCUracy <basic>
              This keyword sets the convergence criterion on CCSD energy.

:kword:`END of input`
  This keyword indicates that there is no more input
  to be read.
  This keyword is *optional*.

:kword:`EXTRapolation`
  This keyword switches on the DIIS extrapolation. This keyword is followed
  by two additional parameters on the next line *n1* and *n2*.

  .. container:: list

    *n1* --- specifies the first iteration, in which DIIS extrapolation procedure
    will start for the first time. This value must not be less then *n2*,
    recommended
    value is 5--7.

    *n2* --- specifies the size of the DIIS procedure, i.e. the number of previous
    CCSD steps which will be used for new prediction. In the present implementation
    *n2* is limited to 2--4.

  This keyword is *optional*. (Default=OFF)

  .. xmldoc:: <KEYWORD MODULE="CCSDT" NAME="EXTR" APPEAR="Extrapolation" KIND="INTS" SIZE="2" LEVEL="BASIC">
              <HELP>
              Switch on DIIS extrapolation, and set two parameters:
              The first iteration to employ DIIS, and the number of previous iterations
              to use for new prediction.
              </HELP>
              </KEYWORD>
              %%Keyword: EXTRapolation <basic>
              Switches the DIIS extrapolation on. Two additional parameters are required
              on the next line: the first iteration to employ DIIS and the number of
              previous iterations to use for new prediction.

:kword:`PRINt`
  The parameter on the next line specifies the level of output printing

  .. container:: list

    0 --- minimal level of printing

    1 --- medium level of printing

    2 --- full output printing (useful for debugging purposes)

  This keyword is *optional*. (Default=0)

  .. xmldoc:: <KEYWORD MODULE="CCSDT" NAME="PRINT" APPEAR="Print level" KIND="CHOICE" LIST="0: Minimal,1: Medium,2: Full" LEVEL="ADVANCED" DEFAULT_VALUE="0">
              <HELP>
              Sets the amount of the program verbosity.
              </HELP>
              </KEYWORD>
              %%Keyword: PRINtlevel <advanced>
              Sets the amount of the program verbosity as 0..2. Default: 0.

:kword:`LOAD`
  This keyword is followed by the line which specifies the
  name of the CCSD amplitudes and energy file. The default name is :file:`RSTART`,
  but it can be changed in CCSD step using :kword:`RESTART` keyword.
  This keyword is *optional*. (Default=:file:`RSTART`)

  .. xmldoc:: <KEYWORD MODULE="CCSDT" NAME="LOAD" APPEAR="Load" KIND="STRING" LEVEL="ADVANCED" DEFAULT_VALUE="RSTART">
              <HELP>
              Alter the file name used to save restart information (Default: RSTART)
              </HELP>
              </KEYWORD>
              %%Keyword: LOAD <advanced>
              This keyword is followed by the line that specifies the name, where the
              restart information was saved.

:kword:`RESTart`
  This keyword defines the restart conditions and modifies the name of the file,
  in which restart information (CC amplitudes, CC energy and the number
  of iterations) is saved. On the following two lines there
  are control key *nn* and the name of restart information storing file
  *name*.

  *nn* --- restart status key

  .. container:: list

    0 --- restart informations will be not saved

    1 --- restart informations will be saved after each iteration in
    *name*.

    2 --- restart run. CC amplitudes and energy will be taken from
    *name* file and the CCSD procedure will continue with
    these values as an estimate.

  *name* --- specifies the restart information storing key. The name is limited
  to 6 characters.

  This keyword is *optional*. (Defaults: *nn* = 1,
  *name* = RSTART)

  .. xmldoc:: <KEYWORD MODULE="CCSDT" NAME="REST" APPEAR="Restart" KIND="STRINGS" SIZE="2" LEVEL="ADVANCED">
              <HELP>
              LINE 1: Choose restart conditions. 0=nothing saved, 1=just save restart info,
              2=also start using restart info. LINE2: The restart file name (at most 6 char).
              </HELP>
              </KEYWORD>
              %%Keyword: RESTart <advanced>
              Followed by two lines.
              LINE 1: Choose restart conditions. 0=nothing saved, 1=just save restart info,
              2=also start using restart info. LINE2: The restart file name (at most 6 char).

:kword:`IOKEy`
  This keyword specifies the input-output file handling.

  .. container:: list

    1 --- Internal Fortran file handling

    2 --- |molcas| DA file handling

  The default (1) is recommended in majority of cases, since when calculating relatively
  large systems with low symmetry, the size of some intermediate files produced may become large,
  what could cause some troubles on 32-bit machines (2 GB file size limit).

  .. xmldoc:: <KEYWORD MODULE="CCSDT" NAME="IOKEY" APPEAR="I/O key" KIND="CHOICE" LIST="0: Fortran I/O,1: MOLCAS I/O" LEVEL="ADVANCED" DEFAULT_VALUE="0">
              <HELP>
              Specify the file type handling. Fortran I/O is default.
              </HELP>
              </KEYWORD>
              %%Keyword: IOKEy <advanced>
              Specifies the file type handling, with Fortran I/O being the default.
              ||1 - Fortran I/O
              ||2 - MOLCAS DA I/O

:kword:`MACHinetyp`
  This keyword specifies which type of matrix multiplication is preferred on a given
  machine. The following line contains two parameters *nn*, *limit*.

  .. container:: list

    *nn* = 1 --- standard multiplication :math:`A B` is preferred

    *nn* = 2 --- transposed multiplication :math:`A^{\text{T}} B` is preferred

  Parameter *limit* specifies the limit for using :math:`A^{\text{T}} B`
  multiplication, when *nn* = 2. (It has no meaning for *nn* = 1.)

  If *size(A)/size(B)* :math:`\geq` *limit* --- standard multiplication is performed,
  *size(A)/size(B)* :math:`<` *limit* --- transposed multiplication is
  performed.

  (*size(A,B)* --- number of elements in matrix A,B).

  Recommended value for *limit* is 2--3.

  Using of transposed matrix (*nn* = 2)
  multiplication may bring some computer time reduction only in special
  cases, however, it requires some additional work space. Default is optimal
  for absolute majority of cases.

  This keyword is *optional*. (Default=1).

  .. xmldoc:: <KEYWORD MODULE="CCSDT" NAME="MACH" APPEAR="Machine" KIND="INTS" SIZE="2" LEVEL="ADVANCED">
              <HELP>
              Use two integers to specify preferred matrix multiply type.
              Usually default is good, and input requires care: Consult manual!
              </HELP>
              </KEYWORD>
              %%Keyword: MACHinetyp <advanced>
              This keyword sets the preferred type of matrix multiplication.
              On the following line n, limit must be specified:
              ||n=1 - standard matrix multiplication is performed
              ||n=2 - A(T)*B matrix multiplication is performed, if
              ||      size(A)/size(B) is less than limit. See manual!

Note, that :kword:`CCSD` and :kword:`CCT` keywords are mutually exclusive.

.. index::
   single: CCSDT; Closed-shell

.. _sec\:ccsdt_cs:

How to run closed shell calculations using ROHF CC codes
--------------------------------------------------------

First of all it should be noted here, that it is not advantageous
to run closed shell calculations using ROHF CC codes, since
in the present implementation it will require the same number of
arithmetical operations and the core and disk space like corresponding
open shell calculations.

Since ROHF CC codes are connected to the output of :program:`RASSCF` code (through the
:file:`JOBIPH` file), it is necessary to run closed shell Hartree--Fock using
the :program:`RASSCF` program. This can be done by setting the number of active orbitals
and electrons to zero (also by including only doubly occupied orbitals into the
active space; this has no advantage but increases the computational effort).
to guarantee the single reference character of the wave function.

The CC program will recognize the closed shell case automatically and will reorganize
all integrals in a required form.
For more information the reader is referred to the tutorials and examples manual.

Below is an input file for :math:`\ce{HF+}` CCSD(T) calculation. ::

  &CCSDT
  Title
   HF(+) CCSD(T) input example
  CCT
  Triples
  3

.. xmldoc:: </MODULE>
