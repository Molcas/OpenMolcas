.. index::
   single: Program; CASPT2
   single: CASPT2

.. _UG\:sec\:caspt2:

:program:`caspt2`
=================

.. only:: html

  .. contents::
     :local:
     :backlinks: none

.. xmldoc:: <MODULE NAME="CASPT2">
            %%Description:
            <HELP>
            The CASPT2 program is used to compute a dynamic correlation correction to
            a RASSCF energy, and optionally to the density matrix and properties.
            It can be regarded as a Møller-Plesset perturbation method through
            second order (MP2), except that is is equally useful for general RASSCF
            root functions, also excited states, regardless of symmetry and multiplicity.
            </HELP>

Second order multiconfigurational perturbation theory is used in the program
:program:`CASPT2` :cite:`Andersson:90,Andersson:92a` to compute the (dynamic)
:index:`correlation energy <single: Perturbation theory; CASPT2>`. The reference state is
usually of the CAS type, but the program has been extended to also accept
RAS reference states :cite:`Malmqvist:2008,Sauri:2011`.
The first step is therefore a RASSCF calculation and the CASPT2
calculation gives a second order estimate of the difference between the RASSCF
and the full CI energy. For calculations using a true RAS reference, benchmark
calculations were reported by Sauri et al. :cite:`Sauri:2011`. For CASSCF
references, the CASPT2 method has been tested in a large number
of applications :cite:`Roos:95a,Roos:96b`. Here follows a brief summary of
results.

.. index::
   single: CASPT2; Precision
   single: Bond length; CASPT2
   single: Bond energy; CASPT2
   single: Heat of reaction; CASPT2

Bond distances are normally obtained with an accuracy of better that 0.01
Å for bonds between first and second row atoms. With the standard Fock
matrix formulation, bond energies are normally underestimated with between 2
and 5 kcal/mol for each bond formed. This is due to a systematic error in the
method :cite:`Andersson:93a`. In every process where the number of paired
electrons is changed, an error of this size will occur for each electron pair.
For example, the singlet-triplet energy difference in the methylene radical
(:math:`\ce{CH2}`) is overestimated with about 3 kcal/mol :cite:`Andersson:92a`. Heats of
reactions for isogyric reactions are predicted with an accuracy of |+-|\ 2
kcal/mol. These results have been obtained with saturated basis sets and all
valence electrons active. The use of smaller basis sets and other types of
active spaces may, of course, affect the error.

These systematic errors have recently been considerably reduced by the
introduction of a modified zeroth order Hamiltonian :cite:`Ghigo:04a`. The method
introduces a shift (the IPEA shift) that modifies the energies of active
orbitals such that they become closer to ionization energies when excited from
and closer to electron affinities when excited out of. The approach has been
tested for 49 diatomic molecules, reducing the mean error in :math:`D_0` from 0.2 to
0.1 eV. For the triply bonded molecules :math:`\ce{N2}`, :math:`\ce{P2}`, and :math:`\ce{As2}` it was reduced
from 0.45 eV to less than 0.15 eV. Similar improvements were obtained for
excitation and ionization energies. The IPEA modified :math:`H_0` (with a shift
parameter of 0.25) is default in |molcas| from version 6.4.

An alternative to IPEA is to use the options, called ":math:`g_1`", ":math:`g_2`", and
":math:`g_3`" (See Ref. :cite:`Andersson:95a`), that stabilizes the energies of the
active orbitals. The remaining error is no longer systematic, and is generally
reduced. For example, the error in the singlet-triplet separation of :math:`\ce{CH2}` is
reduced to 1 kcal/mol :cite:`Andersson:95a`. This option is, however, not
recommended any longer because it has been replaced by the IPEA Hamiltonian.

The CASPT2 method can be used in any case where a valid reference function can
be obtained with the CASSCF method. There is thus no restriction in the number
of open shells or the spin coupling of the electrons. Excited states can be
treated at the same level as ground states. Actually one of the major
successes with the method has been in the calculation of excitation energies.
A large number of applications have been performed for conjugated organic
molecules. Both :index:`Rydberg <single: Rydberg states; CASPT2>` and valence excited states can be treated and the
:index:`error in computed excitation energies <single: Excitation energy; CASPT2>` is normally in the range 0.0--0.2 eV.
Similar results have been obtained for ligand field and charge-transfer
excitations in transition metal compounds. From |molcasvi| it is possible to use
the CASPT2 method in conjunction with the Douglas--Kroll--Hess relativistic
Hamiltonian, which has made possible calculations on heavy element compounds
such a third row transition metal compounds and actinides with accurate results.

The CASPT2 method can also be used in combination with the :program:`FFPT`
program to compute dynamic correlation contributions to properties with good
results in most cases. Numerical gradients are available with the
:program:`slapaf` module.

The CASPT2 method is based on second order perturbation theory. To be
successful, the perturbation should be small. A correct selection of
the active space in the preceding CASSCF calculation is therefore of
utmost importance. All near-degeneracy effects leading to configurations
with large weights must be included at this stage of the calculation.
If this is not done, the first order wave function will contain large
coefficients. When this occurs, the :program:`CASPT2` program issues a warning.
If the energy contribution from such a configuration is large, the results is
not to be trusted and a new selection of the active space should be made.

Especially in calculations on excited states, :index:`intruder <single: Intruders; CASPT2>` states may occur in the
first order wave function. Warnings are then issued by
the program that an energy denominator is small or negative. Such intruder
states often arise from Rydberg orbitals, which have not been included in the
active space. Even if this sometimes leads to large first order CI coefficients,
the contribution to the second order energy is usually very small, since the
interaction with the intruding Rydberg state is small. It might then be
safe to neglect the warning. A safer procedure is to include the Rydberg
orbital into the active space. It can sometimes be deleted from the MO space.

Calculations on compounds with heavy atoms (transition metals, actinides, etc.)
may yield many virtual orbitals with low energies. The interaction energies for
excitations to states where these orbitals are occupied are often very small and
the low denominators can then be removed by a suitable level shift (see below).
But it is always safer to include such orbitals in the active space.

Two keywords have been introduced to deal with this fairly common
situation, for excited states, that weakly coupled intruders cause
spurious singularities, "spikes" in e.g. a potential curve. The two
keywords :kword:`SHIFT` and :kword:`IMAGINARY SHIFT` (mutually exclusive) will introduce a :index:`shift <single: Level shift; CASPT2>`
in the energy denominators,
thus avoiding singularities, and will also correct the energy for the use of
this shift. The net effect is that the energy is almost unaffected except in
the vicinity of the weak singularity, which is removed. The :kword:`SHIFT` keyword adds
a real shift, and the use of this procedure is well tested
:cite:`Roos:95b,Roos:96a`. The :kword:`IMAGINARY SHIFT` adds an imaginary quantity, and
then uses the real value of the resulting second-order energy
:cite:`Forsberg:96`. This offers some advantage, in particular for weak intruder
states.

In some cases, where one can expect strong interaction between different CASSCF
wave functions, it is advisable to use the Multi-State (MS) CASPT2 method
:cite:`Finley:98b`, the extended Multi-State (XMS) method :cite:`Granovsky2011,Shiozaki2011`
or the new extended dynamically weighted CASPT2 :cite:`Battaglia2020`.
A second order effective Hamiltonian is constructed for a
number of CASSCF wave functions obtained in a state-average calculation. This
introduces interaction matrix elements at second order between the different
CASSCF states. The effective Hamiltonian is diagonalized to obtain the final
second order energies. The program also produces a file, :file:`JOBMIX`, with the new
effective zeroth order wave functions, which are linear combinations of the
original CASSCF states. This method has been used successfully to separate
artificially mixed valence and Rydberg states and for transition metal compounds
with low lying excited states of the same symmetry as the ground state.
In the original multi-state method,
perturbed wave functions are computed for each of several root functions,
separately; these are used to compute the effective Hamiltonian.
In the XMS-CASPT2 method, the perturbations are computed with one
common zeroth-order Hamiltonian, and the eigenstates of
the effective Hamiltonian are written onto the :file:`JOBIPH` file rather than used
to generate a new :file:`JOBMIX` file.
The new XDW-CASPT2 method interpolates between the MS and XMS variants,
retaining the advantages of both approaches.

It is clear from the discussion above that it is not a "black box" procedure
to perform CASPT2 calculations on excited states. It is often necessary to
iterate the procedure with modifications of the active space and the selection
of roots in the CASSCF calculation until a stable result is obtained. Normally,
the CASSCF calculations are performed as average calculations over the number
of electronic states of interest, or a larger number of states.
It is imperative that the result is checked
before the CASPT2 calculations are performed. The solutions should contain
the interesting states. If all of them are not there, the number of roots in
the CASSCF calculation has to be increased. Suppose for example, that four
states of a given symmetry are required. Two of them are valence excited states
and two are Rydberg states. A CASSCF calculation is performed as an average
over four roots. Inspection of the solution shows only one valence excited
state, the other three are Rydberg states. After several trials it turns out
that the second valence excited state occurs as root number seven in the
CASSCF calculation. The reason for such a behavior is, of course, the
very different dynamic correlation energies of the valence excited states as
compared to the Rydberg states. It is important that the AO basis set is
chosen to contain a good representation of the Rydberg orbitals, in order to
separate them from the valence excited states. For more details on how to
perform calculations on excited states we refer to the
literature :cite:`Roos:95b,Roos:96a` and :numref:`TUT:sec:excited` of
the examples manual.

The first order wave function is obtained in the :program:`CASPT2` program as an
iterative solution to a large set of linear equations. The size of the
equation system is approximately :math:`n^2 m^2/2` where :math:`n` is the sum of inactive
and active orbitals and :math:`m` is the sum of active and secondary orbitals.
Symmetry will reduce the size with approximately a factor :math:`g_{\text{sym}}`, the
number of irreps of the point group.

:program:`CASPT2` produces a set of :index:`molecular orbitals <pair: Orbitals; CASPT2>` that can be used
as start orbitals for other programs or further calculations.
A minimal CASSCF and CASPT2 gives orbitals and occupation numbers
which can be used to design a proper larger calculation.
By default, the orbitals are natural orbitals obtained from the
density matrix of the (normalized) wave function through first order.
However, the active/active block of that density matrix is not computed
exactly. An approximation has been designed in such a way that the trace
is correct, and the natural occupation numbers of active orbitals are
between zero and two. Due to the approximation, any properties computed
using these orbitals are inexact and can be used only qualitatively. An exact
first order density matrix can be computed but this is more time-consuming. It
is controlled by the keyword :kword:`DENSity`. Use this keyword to compute
properties like dipole moments, etc. The most secure accurate way to do that is.
however, to use finite field perturbation theory (FFPT).

.. index::
   pair: Dependencies; CASPT2

.. _UG\:sec\:caspt2_dependencies:

Dependencies
------------

The :program:`CASPT2` program needs the :file:`JOBIPH` file from a :program:`RASSCF`
calculation, and in addition one- and two-electron integrals and some auxiliary
files from :program:`SEWARD`.

.. index::
   pair: Files; CASPT2

.. _UG\:sec\:caspt2_files:

Files
-----

Input files
...........

:program:`CASPT2` will use the following input
files: :file:`ONEINT`, :file:`ORDINT`, :file:`RUNFILE`, :file:`JOBIPH`
(for more information see :numref:`UG:sec:files_list`).

Output files
............

.. class:: filelist

:file:`PT2ORB`
  Molecular orbitals.

.. index::
   pair: Input; CASPT2

.. _UG\:sec\:caspt2_input:

Input
-----

This section describes the input to the :program:`CASPT2` program, starting with its name: ::

  &CASPT2

.. index::
   pair: Keywords; CASPT2

Keywords
........

.. class:: keywordlist

:kword:`TITLe`
  This keyword is followed by one title line.

  .. xmldoc:: <KEYWORD MODULE="CASPT2" NAME="TITLE" APPEAR="Title" KIND="STRING" LEVEL="BASIC">
              %%Keyword: Title <basic>
              <HELP>
              Enter one title line for this job.
              </HELP>
              </KEYWORD>

:kword:`MULTistate`
  Enter number of root states, and a list of which CI vector from
  the CASSCF calculation to use for each state, for example "``2 1 2``"
  would specify the first and second root.
  Also used for single-state calculations, when the root state is not
  the ground state, for example "``1 2``" would specify the second root.
  The special value "``all``" can be used if all the states included
  in the CASSCF orbital optimization (keyword :kword:`CIRoot` in :program:`RASSCF`)
  are desired.
  Please note thet this is different from an extended multi-state calculation,
  see also :kword:`XMULtistate`.

  .. xmldoc:: <KEYWORD MODULE="CASPT2" NAME="MULTISTATE" APPEAR="Multi-State" KIND="INTS_COMPUTED" SIZE="1" LEVEL="BASIC">
              <ALTERNATE KIND="CUSTOM" />
              %%Keyword: Multistate <basic> GUI:list
              <HELP>
              Enter the number of states for CASPT2 to compute, and a list of numbers
              showing which CASSCF state to use as root state for each.
              Alternatively, enter "all" for all the states included in the CASSCF
              orbital optimization.
              </HELP>
              </KEYWORD>

:kword:`XMULtistate`
  Perform an extended MS-CASPT2 calculation according to :cite:`Granovsky2011,Shiozaki2011`.
  Enter number of root states, and a list of which CI vector from
  the CASSCF calculation to use for each state in the same exact way
  as done for :kword:`MULTistate`. For example "``2 1 2``"
  would specify the first and second root.
  The special value "``all``" can be used if all the states included
  in the CASSCF orbital optimization (keyword :kword:`CIRoot` in :program:`RASSCF`)
  are desired.
  This keyword is mutually exclusive with :kword:`MULTistate`.

  It can be used for an XMS-PDFT calculation (which needs :program:`RASSCF`,
  :program:`CASPT2` and :program:`MCPDFT` modules).
  To carry out an XMS-PDFT calculation, one needs to rotate the SA-CASSCF 
  or SA-RASSCF states to intermediate states (using the same rotation as in 
  XMS-CASPT2), and the rotation matrix can be obtained in the :program:`CASPT2` module 
  with this keyword.
  With this, two files are generated in the scratch directory, :file:`Do_Rotate.txt`,
  which stores the XMS rotation vector, and :file:`H0_Rotate.txt`, which stores the 
  Hamiltonian for the XMS rotated states. If the user wants to skip the 
  expensive perturbation-theory calculation, this keyword can be combined 
  with :kword:`XROH` to skip the perturbation part.

  .. xmldoc:: <KEYWORD MODULE="CASPT2" NAME="XMULTISTATE" APPEAR="Extended Multi-State" KIND="INTS_COMPUTED" SIZE="1" LEVEL="BASIC">
              <ALTERNATE KIND="CUSTOM" />
              %%Keyword: XMultistate <basic> GUI:list
              <HELP>
              Enter the number of states for CASPT2 to compute, and a list of numbers
              showing which CASSCF state to use as root state for each.
              Alternatively, enter "all" for all the states included in the CASSCF
              orbital optimization.
              </HELP>
              </KEYWORD>

:kword:`XROH`
  This keyword can be used in an XMS-PDFT calculation (which needs :program:`RASSCF`, :program:`CASPT2` and :program:`MCPDFT` modules), together with :kword:`XMUL`. 
  When this keyword is used, the :program:`CASPT2` module will not perform perturbation theory calculations; instead, it will only print the rotation matrix and the Hamiltonian matrix of the intermediate states.
  More information can be found on the Minnesota OpenMolcas page\ [#fn1]_.

  .. [#fn1] https://comp.chem.umn.edu/openmolcas/

  .. xmldoc:: <KEYWORD MODULE="CASPT2" NAME="XROH" KIND="SINGLE" LEVEL="BASIC">
              %%Keyword: XROH <basic> 
              <HELP>
              Skips PT2 calculation. Only effective when XMUL is used.
              </HELP>
              </KEYWORD>

:kword:`DWMS`
  It constructs the Fock matrices used in the zeroth-order Hamiltonian
  using dynamically weighted densities. Used in conjunction with
  :kword:`XMULtistate` it performs a XDW-CASPT2 calculation
  according to :cite:`Battaglia2020`.
  It is possible to use this option also with :kword:`MULTistate`, in
  such case the original CASSCF states are used as inputs for the dynamically
  weighted densities, rather than the rotated references as in XDW-CASPT2.
  An integer number for the exponential factor :math:`\zeta` can be specified,
  if not, the default value of 50 is used. By specifying any negative integer
  number, the limit :math:`\zeta\to\infty` is taken, resulting in the
  same weights as in MS-CASPT2. The other limiting case is :math:`\zeta=0`,
  for which equal weights are assigned to all states and thus XDW-CASPT2
  is exactly equivalent to XMS-CASPT2.

  .. xmldoc:: <KEYWORD MODULE="CASPT2" NAME="DWMS" APPEAR="Dynamically Weighted Multi-State" KIND="INT" LEVEL="BASIC">
              %%Keyword: DWMS <basic> GUI:number
              <HELP>
              Enter an integer value specifying the exponent zeta used to
              compute the weights. A negative value corresponds to taking
              the limit to infinity, completely avoiding any mixing of
              the densities.
              </HELP>
              </KEYWORD>

:kword:`IPEAshift`
  This shift corrects the energies of the active orbitals and is
  specified in atomic units. It will be weighted by a function of the
  diagonal density matrix element :math:`D_{pp}`.
  This option is used to modify the standard definition of the
  zeroth order Hamiltonian (:math:`H_0`), which includes an IPEA shift of 0.25
  :cite:`Ghigo:04a`. The modification of :math:`H_0` has been introduced (Nov 2005) to
  reduce the systematic error which leads to a relative overestimation of the
  correlation energy for open shell system. It also reduces the intruder problems.
  Default is to use an IPEA shift of 0.25.

  .. xmldoc:: <KEYWORD MODULE="CASPT2" NAME="IPEA" APPEAR="IPEA shift" KIND="REAL" LEVEL="ADVANCED" DEFAULT_VALUE="0.25">
              %%Keyword: IPEAshift <basic> GUI:number
              <HELP>
              Parameter (Default 0.25), adds a shift dependent on density matrix for active
              orbitals, reducing overestimated correlation energy for open shells.
              </HELP>
              </KEYWORD>

:kword:`IMAGinary`
  Add an imaginary shift to the external part of the zero order
  Hamiltonian. The correlation energy computed is the real part
  of the resulting complex perturbation energy.
  Also, a corrected
  value, obtained by Hylleraas' variational formula, is computed.
  See Ref. :cite:`Forsberg:96`.
  As with the real shift, this option is used to eliminate intruder
  problems.

  .. xmldoc:: <KEYWORD MODULE="CASPT2" NAME="IMAGINARY" APPEAR="Imaginary shift" KIND="REAL" LEVEL="BASIC" DEFAULT_VALUE="0.0">
              %%Keyword: Imaginary <advanced> GUI:number
              <HELP>
              Add an imaginary shift (Default 0.0) to eliminate weak intruders.
              </HELP>
              </KEYWORD>

:kword:`SHIFt`
  Add a shift to the external part of the zero order Hamiltonian.
  See Refs. :cite:`Forsberg:96,Roos:95b,Roos:96b`.
  In addition to the conventionally computed second order energy
  value, another energy obtained by Hylleraas' variational formula
  is computed. This energy is then very close to the unshifted
  energy, except close to singularities due to intruders.
  This option should only be used to eliminate intruder state problems.

  .. xmldoc:: <KEYWORD MODULE="CASPT2" NAME="SHIFT" APPEAR="Real shift" KIND="REAL" LEVEL="ADVANCED" DEFAULT_VALUE="0.0">
              %%Keyword: Shift  <advanced> GUI:number
              <HELP>
              Add a shift to the external part of the zero order Hamiltonian,
              which may shift away weak intruders. Imaginary shift is better.
              </HELP>
              </KEYWORD>

:kword:`AFREeze`
  This keyword is used to select atoms for defining the correlation orbital
  space for the CASPT2 calculation. Assume that you have a large molecule where
  the activity takes place in a limited region (the active site). It could be a
  metal atom with its surrounding ligands. You can then use this option to reduce
  the size of the CASPT2 calculation by freezing and deleting orbitals that have
  only a small population in the active site. An example: The cobalt imido complex
  :math:`\ce{Co^{III}(nacnac)(NPh)}` has 43 atoms. The active site was cobalt and the
  surrounding ligand atoms. Using the AFRE option reduces the time for the CASPT2
  calculation from 3 h to 3 min with a loss of accuracy in relative energies for
  24 electronic states of less than 0.1 eV. The first line after the keyword
  contains the number of selected atoms then the selection thresholds (the
  recommended value is 0.1 or less). An additional line gives the names of the
  atoms as defined in the Seward input. Here is a sample input for the cobalt
  complex mentioned above. ::

    AFREeze
     6 0.10 0.00
     Co N1 N2 C5 C6 C7

  This input means that inactive orbitals with less than 0.1 of the density on
  the active sites will be frozen, while no virtual orbitals will be deleted.

  .. xmldoc:: <KEYWORD MODULE="CASPT2" NAME="AFREEZE" LEVEL="ADVANCED" KIND="CUSTOM">
              %%Keyword: AFREeze <advanced>
              <HELP>
              This keyword is used to select atoms for defining the correlation orbital
              space for the CASPT2 calculation. Inactive orbitals with Mulliken populations
              smaller than a given threshold on the selected atoms will be frozen and
              virtual orbitals will be deleted. The next line give the number of atoms and
              selection thresholds. An additional line gives the names of the atoms as
              defined in the Seward input. Use with care! Not much tested yet, but is very
              effective in reducing the computational time for CASPT2 in large molecules.
              </HELP>
              </KEYWORD>

:kword:`LOVCaspt2`
  "Freeze-and-Delete" type of CASPT2, available only in connection with Cholesky or RI.
  Needs (pseudo)canonical orbitals from RASSCF. An example of input for the keyword :kword:`LOVC` is the following: ::

    LovCASPT2
     0.3
    DoMP2  (or DoEnv)

  In this case, both occupied and virtual orbitals (localized by the program) are divided in two groups: those mainly located on
  the region determined (automatically) by the spatial extent of the active orbitals ("active site"),
  and the remaining ones, which are obviously "outside" this region.
  The value of the threshold (between 0 and 1) is used to perform this selection
  (in the example, 30% of the gross Mulliken population of a given orbital on the active site).
  By default, the CASPT2 calculation is performed only for the correlating orbitals associated with the active site.
  The keyword :kword:`DoMP2` is optional and forces the program to perform also an MP2 calculation on
  the "frozen region".
  Alternatively, one can specify the keyword :kword:`VirAll` in order to use all virtual orbitals as correlating space for the
  occupied orbitals of the active site.
  A third possibility is to use the keyword :kword:`DoEnv` to compute the energy of the environment as total MP2 energy
  minus the MP2 energy of the active site.

  .. xmldoc:: <KEYWORD MODULE="CASPT2" NAME="LOVCASPT2" APPEAR="Localized occupied-virtual CASPT2" LEVEL="ADVANCED" KIND="REAL">
              %%Keyword: LOVC <advanced>
              <HELP>
              "Freeze-and-Delete" type of CASPT2, available only in connection with Cholesky or RI.
              Needs (pseudo)canonical orbitals from RASSCF. An example of input for the keyword LOVC is the following:
              ||
              ||LovCASPT2
              || 0.3
              ||DoMP2  (or DoEnv)
              ||
              In this case, both occupied and virtual orbitals (localized by the program) are divided in two groups: those mainly located on
              the region determined (automatically) by the spatial extent of the active orbitals ("active site"),
              and the remaining ones, which are obviously "outside" this region.
              The value of the threshold (between 0 and 1) is used to perform this selection
              (in the example, 30% of the gross Mulliken population of a given orbital on the active site).
              By default, the CASPT2 calculation is performed only for the correlating orbitals associated with the active site.
              The keyword DoMP2 is optional and forces the program to perform also an MP2 calculation on
              the "frozen region".
              Alternatively, one can specify the keyword VirAll in order to use all virtual orbitals as correlating space for the
              occupied orbitals of the active site.
              A third possibility is to use the keyword DoEnv to compute the energy of the environment as total MP2 energy
              minus the MP2 energy of the active site.
              </HELP>
              </KEYWORD>

  .. xmldoc:: <KEYWORD MODULE="CASPT2" NAME="DOMP2" LEVEL="UNDOCUMENTED" KIND="SINGLE" />

  .. xmldoc:: <KEYWORD MODULE="CASPT2" NAME="VIRALL" LEVEL="UNDOCUMENTED" KIND="SINGLE" />

  .. xmldoc:: <KEYWORD MODULE="CASPT2" NAME="DOENV" LEVEL="UNDOCUMENTED" KIND="SINGLE" />

:kword:`FNOCaspt2`
  Performs a Frozen Natural Orbital (FNO) CASPT2 calculation, available only in combination with Cholesky or RI integral representation.
  Needs (pseudo)canonical orbitals from RASSCF. An example of input for the keyword :kword:`FNOC` is the following: ::

    FNOCaspt2
     0.4
    DoMP2

  The keyword :kword:`FNOC` has one compulsory argument (real number in ]0,1]) specifying the fraction of virtual orbitals
  (in each irrep) to be retained in the FNO-CASPT2 calculation.
  The keyword :kword:`DoMP2` is optional and used to compute the (estimated) correction for the truncation error.

  .. xmldoc:: <KEYWORD MODULE="CASPT2" NAME="FNOCASPT2" APPEAR="Frozen natural orbital CASPT2" LEVEL="ADVANCED" KIND="REAL">
              %%Keyword: FNOC <advanced>
              <HELP>
              Performs a Frozen Natural Orbital (FNO) CASPT2 calculation, available only in combination with Cholesky or RI integral representation.
              Needs (pseudo)canonical orbitals from RASSCF. An example of input for the keyword FNOC is the following:
              ||
              ||FNOCaspt2
              || 0.4
              ||DoMP2
              ||
              The keyword FNOC has one compulsory argument (real number in ]0,1]) specifying the fraction of virtual orbitals
              (in each irrep) to be retained in the FNO-CASPT2 calculation.
              The keyword DoMP2 is optional and used to compute the (estimated) correction for the truncation error.
              </HELP>
              </KEYWORD>

:kword:`FOCKtype`
  Use an alternative Fock matrix. The default Fock matrix is described in
  :cite:`Andersson:90,Andersson:92a` and the other original CASPT2 references.
  The three different modifications named G1, G2 and G3 are described in
  :cite:`Andersson:95a`.
  Note: from 6.4 it is not recommended to use this keyword but
  stay with the IPEA modified :math:`H_0`, which is default.

  .. xmldoc:: <KEYWORD MODULE="CASPT2" NAME="FOCKTYPE" APPEAR="Fock type" KIND="CHOICE" LIST="G1,G2,G3" LEVEL="ADVANCED">
              %%Keyword: Focktype <basic> GUI:select(G1,G2,G3)
              <HELP>
              Present choices: G1, G2, G3. Refers to modified Fock matrices
              (See manual). It is better to use the default IPEA shift.
              </HELP>
              </KEYWORD>

:kword:`FROZen`
  This keyword is used to specify the number of frozen orbitals,
  i.e. the orbitals that are not correlated in the calculation.
  The next line contain the number of frozen orbitals per
  symmetry. The default is to freeze the maximum of those that were frozen in the
  :program:`RASSCF` calculation and the deep core orbitals.
  The frozen orbitals are always the first ones in each symmetry.

  .. xmldoc:: <KEYWORD MODULE="CASPT2" NAME="FROZEN" APPEAR="Frozen" KIND="INTS_LOOKUP" SIZE="NSYM" LEVEL="ADVANCED" MIN_VALUE="0">
              %%Keyword: Frozen <advanced> GUI:list
              <HELP>
              Replace default number of frozen orbitals of each symmetry type with user input.
              Default: Those that were frozen in the RASSCF, or standard table dependent on
              basis set, whichever is larger.
              </HELP>
              </KEYWORD>

:kword:`DELEted`
  This keyword is used to specify the number of deleted orbitals,
  i.e. the orbitals that are not used as correlating orbitals in
  the calculation. The next line contain the number deleted orbitals per symmetry.
  The default is to delete those that were deleted in the :program:`RASSCF`
  calculation.
  The deleted orbitals are always the last ones in each symmetry.

  .. xmldoc:: <KEYWORD MODULE="CASPT2" NAME="DELETED" APPEAR="Deleted" KIND="INTS_LOOKUP" SIZE="NSYM" LEVEL="ADVANCED" MIN_VALUE="0">
              %%Keyword: Deleted <advanced>
              <HELP>
              Replace default number of deleted orbitals of each symmetry type with user input.
              Default: Those that were deleted in the RASSCF.
              </HELP>
              </KEYWORD>

:kword:`DENSity`
  Computes the full density matrix from the first order wave function,
  rather than approximated as is the (faster) default option. Used to
  compute :program:`CASPT2` properties, such as dipole moments, etc.

  .. xmldoc:: <KEYWORD MODULE="CASPT2" NAME="DENSITY" APPEAR="Exact density" KIND="SINGLE" LEVEL="ADVANCED">
              %%Keyword: Density <advanced>
              <HELP>
              Force calculation of accurate density matrix from the
              CASPT2 wave function. Used for dipole moments, etc.
              </HELP>
              </KEYWORD>

:kword:`RFPErt`
  This keyword makes the program add reaction field effects to the energy
  calculation. This is done by adding the reaction field effects to the
  one-electron Hamiltonian as a constant perturbation, i.e. the reaction field
  effect is not treated self consistently. The perturbation is extracted from RUNOLD,
  if that file is not present if defaults to RUNFILE.

  .. xmldoc:: <KEYWORD MODULE="CASPT2" NAME="RFPERT" APPEAR="Reaction field perturbation" KIND="SINGLE" LEVEL="ADVANCED">
              %%Keyword: RFPert  <advanced>
              <HELP>
              Add reaction field from environment as static perturbation
              to the one-electron Hamiltonian.
              The perturbation is extracted from RUNOLD if it exists, otherwise from RUNFILE.
              </HELP>
              </KEYWORD>

:kword:`RLXRoot`
  Specifies which root to be relaxed in a geometry optimization of a
  multi-state CASPT2 wave function. Defaults to the highest root or
  root defined by the same keyword in the :program:`RASSCF` module.

  .. xmldoc:: <KEYWORD MODULE="CASPT2" NAME="RLXROOT" APPEAR="Relaxed root" KIND="INT" LEVEL="ADVANCED" MIN_VALUE="1">
              %%Keyword: RLXRoot <advanced>
              <HELP>
              Which root to use in a geometry optimization of a
              multi-state CASPT2 wave function. Default: root
              defined by RLXROOT in the RASSCF module, if any,
              else the highest root.
              </HELP>
              </KEYWORD>

  .. :kword:`HZERo`
       (No official variants. Perhaps in later versions.)

  .. xmldoc:: <KEYWORD MODULE="CASPT2" NAME="HZERO" KIND="STRING" LEVEL="UNDOCUMENTED" />

:kword:`THREsholds`
  On next line, enter two
  thresholds: for removal of zero-norm components in the
  first-order perturbed wave function, and for removal of near linear
  dependencies in the first-order perturbed wave function. Default
  values are 1.0d-10 and 1.0d-08 respectively.

  .. xmldoc:: <KEYWORD MODULE="CASPT2" NAME="THRESHOLD" APPEAR="Thresholds" KIND="REALS" SIZE="2" LEVEL="ADVANCED" DEFAULT_VALUES="1.0D-10,1.0D-8" MIN_VALUE="0.0">
              %%Keyword: Thresholds  <advanced>
              <HELP>
              The first threshold is for removing redundant excitations, the second is
              for removing linear dependences of standardized linear equation system.
              Default: 1.0d-10 and 1.0d-08.
              </HELP>
              </KEYWORD>

:kword:`MAXIter`
  On next line, enter the maximum allowed number of iterations
  in a procedure for solving a system of
  linear equations using a conjugate gradient method. Default is 20.
  A gradient norm is reported. This gradient is a residual error from the
  CASPT2 equation solution and should be small, else the number of iterations
  must be increased.

  .. xmldoc:: <KEYWORD MODULE="CASPT2" NAME="MAXITER" APPEAR="Maximum iterations" KIND="INT" LEVEL="ADVANCED" DEFAULT_VALUE="20" MIN_VALUE="0">
              %%Keyword: MaxIter <advanced>
              <HELP>
              The maximum allowed number of iterations.
              (Zero iterations gives diagonal approximation. Default:20)
              </HELP>
              </KEYWORD>

:kword:`CONVergence`
  On next line, enter the convergence threshold for the procedure described above.
  The iterative procedure is repeated until the norm of the residual
  (RNORM) is less than this convergence threshold. Default is 1.0d-06.

  .. xmldoc:: <KEYWORD MODULE="CASPT2" NAME="CONV" APPEAR="Convergence" KIND="REAL" LEVEL="ADVANCED" DEFAULT_VALUE="1.0D-6" MIN_VALUE="0.0">
              %%Keyword: Convergence <advanced>
              <HELP>
              Convergence threshold for norm of residual vector. Default 1.0d-06
              </HELP>
              </KEYWORD>

:kword:`NOMIx`
  Normally, an (X)MS-CASPT2 calculation produces a new jobiph file named :file:`JOBMIX`.
  It has the same CASSCF wave functions as the original ones, except that those CI vectors
  that were used in the (Extended) Multi-State CASPT2 calculation have been mixed,
  using the eigenvectors of the effective Hamiltonian matrix as transformation coefficients.
  Keyword :kword:`NOMIX` prevents creation of this :file:`JOBMIX` file.

  .. xmldoc:: <KEYWORD MODULE="CASPT2" NAME="NOMIX" APPEAR="No JobMix" KIND="SINGLE" LEVEL="ADVANCED">
              %%Keyword: NoMix <advanced>
              <HELP>
              Do not produce a JobMix file, even if this is a multi-state calculation.
              </HELP>
              </KEYWORD>

:kword:`NOMUlt`
  This keyword removes the multi-state part of the calculation and only runs a
  series of independent CASPT2 calculations for the roots specified by the
  :kword:`MULTistate` or :kword:`XMULtistate` keyword. Useful when many roots are required,
  but multi-state is not needed, or desired. Note that a :file:`JOBMIX` file is produced
  anyway, but the vectors will not be mixed, and the energies will be single-state CASPT2
  energies. If used with the :kword:`XMULtistate` keyword, the zeroth-order Hamiltonian
  will be constructed with the state-average density and therefore will be the same for
  all the states.

  .. xmldoc:: <KEYWORD MODULE="CASPT2" NAME="NOMULTI" APPEAR="No Multi-State" KIND="SINGLE" LEVEL="BASIC">
              %%Keyword: NoMultistate <basic>
              <HELP>
              Do just a series of independent CASPT2 runs, without any multi-state coupling.
              Useful when many roots are required, but multi-state is not needed, or desired.
              </HELP>
              </KEYWORD>

:kword:`ONLY`
  This keyword requires the :kword:`MULTistate` or :kword:`XMULtistate` keyword,
  and is followed by an integer specifying one of the roots.
  In a (Extended) Multistate calculation, it requests to compute the energy of
  only the specified root. However, the effective Hamiltonian coupling terms
  between this root and all the others included in the (Extended) Multistate
  treatment will be computed and printed out.
  This output will be used in a subsequent calculation, in conjunction
  with the :kword:`EFFE` keyword.

  .. xmldoc:: <KEYWORD MODULE="CASPT2" NAME="ONLY" APPEAR="Only root" KIND="INT" LEVEL="ADVANCED">
              %%Keyword: ONLY <advanced>
              <HELP>
              This keyword requires the MULTistate or XMULtistate keyword, and is
              followed by an integer specifying one of the roots.
              In a Multistate calculation, it requests to compute the energy of only
              the specified root. However, the effective Hamiltonian coupling terms
              between this root and all the others included in the Multistate
              treatment will be computed and printed out.
              This output will be used in a subsequent calculation, in conjunction
              with the EFFE keyword.
              </HELP>
              </KEYWORD>

:kword:`EFFE`
  This keyword requires the :kword:`MULTistate` or :kword:`XMULtistate` keyword.
  It is followed by the number of states and a matrix of real numbers,
  specifying the effective Hamiltonian couplings, as provided in a previous
  calculation using the :kword:`ONLY` keyword.
  In a (Extended) Multistate calculation over, e.g., 3 states, 3 separate
  calculations with the :kword:`ONLY` keyword will be performed, possibly
  on separate computing nodes, so as to speed up the overall process.
  The three couplings vectors will be given to the :kword:`EFFE`
  keyword in matrix form, i.e. the first column is made by the
  couplings of the first computed root, etc.
  The program will then quickly compute the (Extended) Multistate energies.

  .. xmldoc:: <KEYWORD MODULE="CASPT2" NAME="EFFE" APPEAR="Effective Hamiltonian couplings" KIND="CUSTOM" LEVEL="ADVANCED">
              %%Keyword: EFFE <advanced>
              <HELP>
              This keyword requires the MULTistate or XMULtistate keyword. It is
              followed by the number of states and a matrix of real numbers,
              specifying the effective Hamiltonian couplings, as provided in
              a previous calculation using the ONLY keyword.
              In a (Extended) Multistate calculation over, e.g., 3 states, 3 separate
              calculations with the ONLY keyword will be performed, possibly
              on separate computing nodes, so as to speed up the overall process.
              The three couplings vectors will be given to the EFFE
              keyword in matrix form, i.e. the first column is made by the
              couplings of the first computed root, etc.
              The program will then quickly compute the (Extended) Multistate energies.
              </HELP>
              </KEYWORD>

:kword:`NOORbitals`
  In calculations with very many orbitals, use this keyword to skip the
  printing of the MO orbitals.

  .. xmldoc:: <KEYWORD MODULE="CASPT2" NAME="NOORBITALS" APPEAR="No orbitals" KIND="SINGLE" LEVEL="BASIC">
              %%Keyword: NoOrbitals <basic>
              <HELP>
              Skip printing of the MO orbitals.
              </HELP>
              </KEYWORD>

:kword:`PROPerties`
  Normally, a CASPT2 calculation does not produce any density matrix,
  natural orbitals or properties in order to save time and memory
  (especially for large calculations).
  Keyword :kword:`PROP` activates these calculations, at the expense of (some)
  extra time and memory (especially if used together with the :kword:`DENS` keyword).

  .. xmldoc:: <KEYWORD MODULE="CASPT2" NAME="PROPERTIES" APPEAR="Properties" KIND="SINGLE" LEVEL="BASIC">
              %%Keyword: Properties <basic>
              <HELP>
              Compute (approximate) density matrix, natural orbitals and properties.
              </HELP>
              </KEYWORD>

:kword:`NOTRansform`
  This keyword specifies that the wave function should not be transformed
  to use quasi-canonical orbitals, even if :program:`CASPT2` does not know if this
  was done or not and by default would do such a transformation.
  Effectively, the Fock matrix is replaced by a diagonal
  approximation in the input orbital system.

  .. xmldoc:: <KEYWORD MODULE="CASPT2" NAME="NOTRANSFORM" APPEAR="No transform" KIND="SINGLE" LEVEL="ADVANCED">
              %%Keyword: NoTransform <advanced>
              <HELP>
              Prevent transformation to pseudo-canonical orbitals, even if CASPT2
              would assumed this is needed (Default is: transform when assumed necessary.)
              </HELP>
              </KEYWORD>

:kword:`TRANsform`
  This keyword specifies that the wave function should be transformed
  to use pseudo-canonical orbitals, even if this was specified
  as option to the CASSCF calculation and should be unnecessary.
  (Default is: to transform when necessary, and not else.)

  .. xmldoc:: <KEYWORD MODULE="CASPT2" NAME="TRANSFORM" APPEAR="Transform" KIND="SINGLE" LEVEL="ADVANCED">
              %%Keyword: Transform <advanced>
              <HELP>
              Demand transformation to pseudo-canonical orbitals
              even if this was specified as option of CASSCF so it ought to
              be unnecessary. (Default is: transform only when assumed necessary.)
              </HELP>
              </KEYWORD>

:kword:`OFEMbedding`
  Adds an Orbital-Free Embedding potential to the Hamiltonian. Available only in combination with Cholesky or RI integral representation.
  No arguments required. The runfile of the environment subsystem (:file:`AUXRFIL`) must be available.

  .. xmldoc:: <KEYWORD MODULE="CASPT2" NAME="OFEMBEDDING" APPEAR="Orbital-free embedding" KIND="SINGLE" LEVEL="ADVANCED">
              %%Keyword: OFEM <advanced>
              <HELP>
              Adds an Orbital-Free Embedding potential to the Hamiltonian. Available only in combination with Cholesky or RI integral representation.
              No arguments required. The runfile of the environment subsystem (AUXRFIL) must be available.
              </HELP>
              </KEYWORD>

:kword:`GHOStdelete`
  Excludes from PT2 treatment orbitals localized on ghost atoms. A threshold for this selection must be specified.

  .. xmldoc:: <KEYWORD MODULE="CASPT2" NAME="GHOSTDELETE" APPEAR="Ghost delete" KIND="REAL" LEVEL="ADVANCED">
              %%Keyword: GHOS <advanced>
              <HELP>
              Excludes from PT2 treatment orbitals localized on ghost atoms. A threshold for this selection must be specified.
              </HELP>
              </KEYWORD>

:kword:`OUTPut`
  Use this keyword, followed by any of the words :kword:`BRIEF`, :kword:`DEFAULT`, or :kword:`LONG`, to
  control the extent of orbital listing.
  :kword:`BRIEF` gives a very short orbital listing,
  :kword:`DEFAULT` a normal output, and :kword:`LONG` a detailed listing.

  .. xmldoc:: <KEYWORD MODULE="CASPT2" NAME="OUTPUT" APPEAR="Orbital printing" KIND="CHOICE" LIST="BRIEF,DEFAULT,LONG" LEVEL="BASIC">
              %%Keyword: OUTPut <basic>
              <HELP>
              BRIEF gives a very short orbital listing,
              DEFAULT a normal output, and LONG a detailed listing.
              </HELP>
              </KEYWORD>

:kword:`PRWF`
  This keyword is used to specify the threshold for printing the
  CI coefficients, default is 0.05.

  .. xmldoc:: <KEYWORD MODULE="CASPT2" NAME="PRWF" APPEAR="Print threshold" KIND="REAL" LEVEL="ADVANCED" DEFAULT_VALUE="0.05" MIN_VALUE="0.0" MAX_VALUE="1.0">
              %%Keyword: PRWF <advanced>
              <HELP>
              Threshold for printing CI coefficients. Default 0.05.
              </HELP>
              </KEYWORD>

:kword:`PRSD`
  This keyword is used to request that not only CSFs are printed with
  the CI coefficients, but also the determinant expansion.

  .. xmldoc:: <KEYWORD MODULE="CASPT2" NAME="PRSD" APPEAR="Print determinant expansion" KIND="SINGLE" LEVEL="ADVANCED">
              %%Keyword: PRSD <advanced>
              <HELP>
              Activate printing of CSFs in terms of determinants.
              </HELP>
              </KEYWORD>

:kword:`CHEMps2`
  Activate DMRG-CASPT2 calculation with |molcas|--CheMPS2 interface.
  The keyword :kword:`3RDM` must be used in :program:`RASSCF`.
  The program will skip the calculations of the :math:`n`-particle reduced density matrix.
  Note that multi-state calculations are not supported, the calculation will run but produce wrong CASPT2 total energy.
  Always specify :kword:`MULTi` = 1 *iroot*, where *iroot* is the root index.

  .. xmldoc:: <KEYWORD MODULE="CASPT2" NAME="CHEMPS2" APPEAR="DMRG-CASPT2 (CheMPS2)" KIND="SINGLE" LEVEL="BASIC">
              %%Keyword: CHEMps2 <basic>
              <HELP>
              Activate DMRG-CASPT2 calculation with Molcas-CheMPS2 interface.
              </HELP>
              </KEYWORD>

:kword:`CUMUlant`
  Activate DMRG-cu(4)-CASPT2 calculation with |molcas|--Block interface.
  The keyword :kword:`3RDM` must be used in :program:`RASSCF`.
  The program will skip the calculations of the 3-particle reduced density matrix and approximate
  the 4-particle reduced density matrix.

  .. xmldoc:: <KEYWORD MODULE="CASPT2" NAME="CUMULANT" APPEAR="DMRG-cu(4)-CASPT2 (Block)" KIND="SINGLE" LEVEL="BASIC">
              %%Keyword: CUMUlant <basic>
              <HELP>
              Activate DMRG-cu(4)-CASPT2 calculation with Molcas-Block interface.
              </HELP>
              </KEYWORD>

The given default values for the keywords
:kword:`Convergence` and
:kword:`Thresholds` normally give a second order energy which is correct
in eight decimal places.

Input example
.............

::

  &CASPT2
  Title
   The water molecule
  Density matrix

The CASPT2 energy and density matrix is computed for the water molecule with the
O(1s) orbital frozen. The standard IPEA-:math:`H_0` is used.

Input example for SS-DMRG-CASPT2 with |molcas|--CheMPS2 interface

::

  &RASSCF
  Title    = Water molecule. Ground state
  Spin     = 1
  Symmetry = 1
  Inactive = 2 0 1 0
  Ras2     = 2 2 0 0
  DMRG     = 500
  LUMOrb
  3RDM

  &CASPT2
  CHEMps2

Input example for SA-DMRG-CASPT2 with |molcas|--CheMPS2 interface

::

  &RASSCF
  Title    = Water molecule. Averaging two states
  Spin     = 1
  Symmetry = 1
  Inactive = 2 0 1 0
  Ras2     = 2 2 0 0
  CIROot   = 2 2 1
  DMRG     = 500

  &RASSCF
  Title    = Ground state
  Spin     = 1
  Symmetry = 1
  Inactive = 2 0 1 0
  Ras2     = 2 2 0 0
  CIROot   = 1 1 ; 1
  DMRG     = 500
  LUMOrb
  CIONly
  3RDM

  &CASPT2
  CHEMps2
  MULTistate = 1 1

  &RASSCF
  Title    = First excited state
  Spin     = 1
  Symmetry = 1
  Inactive = 2 0 1 0
  Ras2     = 2 2 0 0
  CIROot   = 1 2 ; 2
  DMRG     = 500
  CIONly
  3RDM

  &CASPT2
  CHEMps2
  MULTistate = 1 2

.. xmldoc:: <KEYWORD MODULE="CASPT2" NAME="LROOT" KIND="INT" LEVEL="UNDOCUMENTED" />

.. xmldoc:: <KEYWORD MODULE="CASPT2" NAME="FILE" KIND="STRING" LEVEL="UNDOCUMENTED" />

.. xmldoc:: <KEYWORD MODULE="CASPT2" NAME="RHSD" KIND="SINGLE" LEVEL="UNDOCUMENTED" />

.. xmldoc:: <KEYWORD MODULE="CASPT2" NAME="WTHR" KIND="REALS" SIZE="3" LEVEL="UNDOCUMENTED" />

.. xmldoc:: <KEYWORD MODULE="CASPT2" NAME="G1SE" KIND="SINGLE" LEVEL="UNDOCUMENTED" />

.. xmldoc:: </MODULE>
