Computing Excited States
========================

The calculation of electronic excited states is typically a multiconfigurational problem, and
therefore it should preferably be treated with multiconfigurational methods such as CASSCF and
CASPT2. We can start this section by computing the low-lying electronic states of the
acrolein molecule at the CASSCF level and using a minimal
basis set. The standard file with cartesian coordinates is:

.. extractfile:: problem_based_tutorials/acrolein.xyz

  8
  Angstrom
   O      -1.808864   -0.137998    0.000000
   C       1.769114    0.136549    0.000000
   C       0.588145   -0.434423    0.000000
   C      -0.695203    0.361447    0.000000
   H      -0.548852    1.455362    0.000000
   H       0.477859   -1.512556    0.000000
   H       2.688665   -0.434186    0.000000
   H       1.880903    1.213924    0.000000

We shall carry out State-Averaged (SA) CASSCF calculations, in which one single
set of molecular orbitals is used to compute all the states of a given spatial
and spin symmetry. The obtained density matrix is the average for all states
included, although each state will have its own set of optimized CI
coefficients. Different weights can be considered for each of the states,
but this should not be used except in very special cases by experts. It is
better to let the CASPT2 method to handle that. The use of a SA-CASSCF
procedure has an great advantage. For example, all states in a SA-CASSCF
calculation are orthogonal to each other, which is not necessarily true for
state specific calculations. Here, we shall include five states of singlet
character the calculation. As no symmetry is invoked all the states belong by
default to the first symmetry, including the ground state.

.. extractfile:: problem_based_tutorials/CASSCF.excited.acrolein.input

  *CASSCF SA calculation on five singlet excited states in acrolein
  *File: CASSCF.excited.acrolein
  *
  &GATEWAY
    Title= Acrolein molecule
    coord = acrolein.xyz; basis = STO-3G; group = c1
  &SEWARD; &SCF
  &RASSCF
    LumOrb
    Spin= 1; Nactel= 6 0 0; Inactive= 12; Ras2= 5
    CiRoot= 5 5 1
  &GRID_IT
    All

We have used as active all the :math:`\pi` and :math:`\pi^*` orbitals, two bonding and
two antibonding :math:`\pi` orbitals with four electrons and in addition the oxygen
lone pair (:math:`n`). Keyword :kword:`CiRoot` informs the program that we want to
compute a total of five states, the ground state and the lowest four excited
states at the CASSCF level and that all of them should have the same weight in
the average procedure. Once analyzed we find that the calculation has provided,
in this order, the ground state, two :math:`n\to\pi^*` states, and two :math:`\pi\to\pi^*` states.
It is convenient to add the :program:`GRID_IT` input in order to be able to use
the :program:`LUSCUS` interface for the analysis of the orbitals and the occupations
in the different electronic states. Such an analysis should always be made in
order to understand the nature of the different excited states.
In order to get a more detailed analysis of the nature of the obtained states it is
also possible to obtain in a graphical way the charge density differences between
to states, typically the difference between the ground and an excited state. The
following example creates five different density files:

.. extractfile:: problem_based_tutorials/CASSCF.excited_grid.acrolein.input

  *CASSCF SA calculation on five singlet excited states in acrolein
  *File: CASSCF.excited_grid.acrolein
  *
  &GATEWAY
    Title= Acrolein molecule
    coord= acrolein.xyz; basis= STO-3G; group= c1
  &SEWARD; &SCF
  &RASSCF
   LumOrb
   Spin= 1; Nactel= 6 0 0; Inactive= 12; Ras2= 5
   CiRoot= 5 5 1
   OutOrbital
   Natural= 5
  &GRID_IT
   FILEORB = $Project.RasOrb.1
   NAME = 1; All
  &GRID_IT
   FILEORB = $Project.RasOrb.2
   NAME = 2; All
  &GRID_IT
   FILEORB = $Project.RasOrb.3
   NAME = 3; All
  &GRID_IT
   FILEORB = $Project.RasOrb.4
   NAME = 4; All
  &GRID_IT
   FILEORB = $Project.RasOrb.5
   NAME = 5; All

In :program:`GRID_IT` input we have included all orbitals. It is, however,
possible and in general recommended to restrict the calculation to certain
sets of orbitals. How to do this is described in the input manual for
:program:`GRID_IT`.

Simple math operations can be performed with grids of the same size,
for example, :program:`LUSCUS` can be used to display the difference
between two densities.

CASSCF wave functions are typically good enough, but this is not the case for
electronic energies, and the dynamic correlation effects have to be included,
in particular here with the CASPT2 method. The proper input is prepared, again
including :program:`SEWARD` and :program:`RASSCF` (unnecessary if they were
computed previously), adding a :program:`CASPT2` input with the keyword
:kword:`MultiState` set to 5 1 2 3 4 5. The :program:`CASPT2` will perform four
consecutive single-state (SS) CASPT2 calculations using the SA-CASSCF roots computed
by the :program:`RASSCF` module. At the end, a multi-state CASPT2 calculation
will be added in which the five SS-CASPT2 roots will be allowed to interact.
The final MS-CASPT2 solutions, unlike the previous SS-CASPT2 states, will be
orthogonal. The :kword:`FROZen` keyword is put here as a reminder. By
default the program leaves the core orbitals frozen.

.. extractfile:: problem_based_tutorials/CASPT2.excited.acrolein.input

  *CASPT2 calculation on five singlet excited states in acrolein
  *File: CASPT2.excited.acrolein
  *
  &GATEWAY
   Title= Acrolein molecule
   coord = acrolein.xyz; basis = STO-3G; group= c1
  &SEWARD; &SCF
  &RASSCF
   Spin= 1; Nactel= 6 0 0; Inactive= 12; Ras2= 5
   CiRoot= 5 5 1
  &GRID_IT
   All
  &CASPT2
   Multistate= 5 1 2 3 4 5
   Frozen= 4

Apart from energies and state properties it is quite often necessary to compute
state interaction properties such as transition dipole moments, Einstein coefficients,
and many other. This can be achieved with the :program:`RASSI` module, a powerful
program which can be used for many purposes
(see section :ref:`UG:sec:rassi`). We can
start by simply computing the basic interaction properties

.. extractfile:: problem_based_tutorials/CASSI.excited.acrolein.input

  *RASSI calculation on five singlet excited states in acrolein
  *File: RASSI.excited.acrolein
  *
  &GATEWAY
   Title= Acrolein molecule
   coord = acrolein.xyz; basis = STO-3G; group = c1
  &SEWARD; &SCF
  &RASSCF
   LumOrb
   Spin= 1; Nactel= 6 0 0; Inactive= 12; Ras2= 5
   CiRoot= 5 5 1
  &CASPT2
   Frozen = 4
   MultiState= 5 1 2 3 4 5

  >>COPY $Project.JobMix JOB001

  &RASSI
   Nr of JobIph
   1 5
   1 2 3 4 5
   EJob

Oscillator strengths for the computed transitions and Einstein coefficients are
compiled at the end of the :program:`RASSI` output file. To obtain these values,
however, energy differences have been used which are obtained from the previous
CASSCF calculation. Those energies are not accurate because they do not include
dynamic correlation energy and it is better to substitute them by properly
computed values, such those at the CASPT2 level. This is achieved with the
keyword :kword:`Ejob`.
More information is available
in section :ref:`TUT:sec:rassi_thio`.

Now a more complex case. We want to compute vertical singlet-triplet gaps from
the singlet ground state of acrolein to different, up to five, triplet excited
states. Also, interaction properties are requested. Considering that the spin
multiplicity differs from the ground to the excited states, the spin Hamiltonian
has to be added to our calculations and the :program:`RASSI` program takes charge
of that. It is required first, to add in the :program:`SEWARD` input the keyword
:kword:`AMFI`, which introduces the proper integrals required, and to the
:program:`RASSI` input the keyword :kword:`SpinOrbit`. Additionally, as we want
to perform the calculation sequentially and :program:`RASSI` will read from
two different wave function calculations, we need to perform specific links
to save the information. The link to the first :program:`CASPT2` calculation
will saved in file :file:`$Project.JobMix.S` the data from the :program:`CASPT2`
result of the ground state, while the second link before the second :program:`CASPT2`
run will do the same for the triplet states. Later, we link these files as
:file:`JOB001` and :file:`JOB002` to become input files for :program:`RASSI`.
In the :program:`RASSI` input :kword:`NrofJobIph` will be set to two, meaning
two :file:`JobIph` or :file:`JobMix` files, the first containing one root (the ground
state) and the second five roots (the triplet states). Finally, we have added
:kword:`EJob`, which will read the CASPT2 (or MS-CASPT2) energies from the
:file:`JobMix` files to be incorporated to the :program:`RASSI` results.
The magnitude of properties computed with spin-orbit coupling (SOC) depends
strongly on the energy gap, and this has to be computed at the highest possible
level, such as CASPT2.

.. extractfile:: problem_based_tutorials/CASPT2.S-T_gap.acrolein.input

  *CASPT2/RASSI calculation on singlet-triplet gaps in acrolein
  *File: CASPT2.S-T_gap.acrolein
  *
  &GATEWAY
   Title= Acrolein molecule
   coord = acrolein.xyz; basis = STO-3G; group= c1
  &SEWARD
   AMFI
  &SCF
  &RASSCF
   Spin= 1; Nactel= 6 0 0; Inactive= 12; Ras2= 5
   CiRoot= 1 1 1
  &CASPT2
   Frozen= 4
   MultiState= 1 1
  >>COPY $Project.JobMix JOB001
  &RASSCF
   LumOrb
   Spin= 3; Nactel= 6 0 0; Inactive= 12; Ras2= 5
   CiRoot= 5 5 1
  &CASPT2
   Frozen= 4
   MultiState= 5 1 2 3 4 5
  >>COPY $Project.JobMix JOB002
  &RASSI
   Nr of JobIph= 2 1 5; 1; 1 2 3 4 5
   Spin
   EJob

As here with keyword :kword:`AMFI`,
when using command :kword:`Coord` to build a :program:`SEWARD` input
and we want to introduce other keywords, it is enough if we place them
after the line corresponding to :kword:`Coord`.
Observe that the nature of the triplet states obtained is in sequence one
:math:`n\pi^*`, two :math:`\pi\pi^*`, and two :math:`n\pi^*`. The :program:`RASSI` output is
somewhat complex to analyze, but it makes tables summarizing oscillator
strengths and Einstein coefficients, if those are the magnitudes of interest.
Notice that a table is first done with the spin-free states, while the final
table include the spin-orbit coupled eigenstates (in the CASPT2 energy order
here), in which each former triplet state has three components.

In many cases working with symmetry will help us to perform calculations
in quantum chemistry. As it is a more complex and delicate problem we direct
the reader to the examples section in this manual. However, we include here
two inputs that can help the beginners. They are based on trans-1,3-butadiene,
a molecule with a :math:`C_{2h}` ground state. If we run the next input, the
:program:`SEWARD` and :program:`SCF` outputs will help us to understand how
orbitals are classified by symmetry, whereas reading the :program:`RASSCF` output
the structure of the active space and states will be clarified.

.. extractfile:: problem_based_tutorials/CASSCF.excited.tButadiene.1Ag.input

  *CASSCF SA calculation on 1Ag excited states in tButadiene
  *File: CASSCF.excited.tButadiene.1Ag
  *
  &SEWARD
    Title= t-Butadiene molecule
    Symmetry= Z XYZ
  Basis set
  C.STO-3G...
  C1   -3.2886930 -1.1650250 0.0000000  Bohr
  C2   -0.7508076 -1.1650250 0.0000000  Bohr
  End of basis
  Basis set
  H.STO-3G...
  H1   -4.3067080  0.6343050 0.0000000  Bohr
  H2   -4.3067080 -2.9643550 0.0000000  Bohr
  H3    0.2672040 -2.9643550 0.0000000  Bohr
  End of basis

  &SCF

  &RASSCF
   LumOrb
   Title= tButadiene molecule (1Ag states). Symmetry order (ag bg bu au)
   Spin= 1; Symmetry= 1; Nactel= 4 0 0; Inactive= 7 0 6 0; Ras2= 0 2 0 2
   CiRoot= 4 4 1

  &GRID_IT
   All

Using the next input will give information about states of a different symmetry.
Just run it as a simple exercise.

.. extractfile:: problem_based_tutorials/CASSCF.excited.tButadiene.1Bu.input

  *CASSCF SA calculation on 1Bu excited states in tButadiene
  *File: CASSCF.excited.tButadiene.1Bu
  *
  &SEWARD
   Title= t-Butadiene molecule
   Symmetry= Z XYZ
  Basis set
  C.STO-3G...
  C1   -3.2886930 -1.1650250 0.0000000  Bohr
  C2   -0.7508076 -1.1650250 0.0000000  Bohr
  End of basis
  Basis set
  H.STO-3G...
  H1   -4.3067080  0.6343050 0.0000000  Bohr
  H2   -4.3067080 -2.9643550 0.0000000  Bohr
  H3    0.2672040 -2.9643550 0.0000000  Bohr
  End of basis

  &SCF

  &RASSCF
   FileOrb= $Project.ScfOrb
   Title= tButadiene molecule (1Bu states). Symmetry order (ag bg bu au)
   Spin= 1; Symmetry= 1; Nactel= 4 0 0; Inactive= 7 0 6 0
   Ras2= 0 2 0 2
   CiRoot= 4 4 1
  >COPY $Project.RasOrb $Project.1Ag.RasOrb
  >COPY $Project.JobIph JOB001

  &GRID_IT
   Name= $Project.1Ag.lus
   All

  &RASSCF
   FileOrb= $Project.ScfOrb
   Title= tButadiene molecule (1Bu states). Symmetry order (ag bg bu au)
   Spin= 1; Symmetry= 3; Nactel= 4 0 0; Inactive= 7 0 6 0; Ras2= 0 2 0 2
   CiRoot= 2 2 1
  >COPY $Project.RasOrb $Project.1Bu.RasOrb
  >COPY $Project.JobIph JOB002

  &GRID_IT
   Name= $Project.1Bu.lus
   All

  &RASSI
   NrofJobIph= 2 4 2; 1 2 3 4; 1 2

Structure optimizations can be also performed at the CASSCF, RASSCF or CASPT2
levels. Here we shall optimize the second singlet state in the first (here the
only) symmetry for acrolein at the SA-CASSCF level. It is strongly recommended
to use the State-Average option and avoid single state CASSCF calculations for
excited states. Those states are non-orthogonal with the ground state and
are typically heavily contaminated. The usual set of input commands will be
prepared, with few changes. In the :program:`RASSCF` input two states will
be simultaneously computed with equal weight (:kword:`CiRoot` 2 2 1), but,
in order to get accurate gradients for a specific root (not an averaged one),
we have to add :kword:`Rlxroot` and set it to two, which is, among the
computed roots, that we want to optimize. The proper density matrix will be
stored. The :program:`MCLR` program optimizes, using a perturbative approach,
the orbitals for the specific root (instead of using averaged orbitals), but
the program is called automatically and no input is needed.

.. extractfile:: problem_based_tutorials/CASSCF.excited_state_optimization.acrolein.input

  *CASSCF excited state optimization in acrolein
  *File: CASSCF.excited_state_optimization.acrolein
  *
   &GATEWAY
  Title= acrolein minimum optimization in excited state 2
  Basis set
  O.STO-3G...2s1p.
  O1       1.608542      -0.142162       3.240198 Angstrom
  End of basis
  Basis set
  C.STO-3G...2s1p.
  C1      -0.207776       0.181327      -0.039908 Angstrom
  C2       0.089162       0.020199       1.386933 Angstrom
  C3       1.314188       0.048017       1.889302 Angstrom
  End of basis
  Basis set
  H.STO-3G...1s.
  H1       2.208371       0.215888       1.291927 Angstrom
  H2      -0.746966      -0.173522       2.046958 Angstrom
  H3      -1.234947       0.213968      -0.371097 Angstrom
  H4       0.557285       0.525450      -0.720314 Angstrom
  End of basis
  >>> Do while

   &SEWARD

  >>> If ( Iter = 1 ) <<<

   &SCF
  Title= acrolein minimum optimization

  >>> EndIf <<<

   &RASSCF
  LumOrb
  Title= acrolein
  Spin= 1; nActEl= 4 0 0; Inactive= 13; Ras2= 4
  CiRoot= 2 2 1
  Rlxroot= 2

   &SLAPAF

  >>> EndDo

In case of performing a :program:`CASPT2` optimization for an excited
state, still the SA-CASSCF approach can be used to generate the reference
wave function, but keyword :kword:`Rlxroot` and the use of the :program:`MCLR` program
are not necessary, because :program:`CASPT2` takes care of selecting
the proper root (the last one).

A very useful tool recently included in |molcas| is the possibility to
compute minimum energy paths (MEP), representing steepest descendant minimum
energy reaction paths which are built through a series of geometry optimizations,
each requiring the minimization of the potential energy on a hyperspherical
cross section of the PES centered on a given reference geometry and characterized
by a predefined radius. One usually starts the calculation from a high energy reference
geometry, which may correspond to the Franck--Condon (FC) structure on an excited-state PES
or to a transition structure (TS). Once the first lower energy optimized structure is
converged, this is taken as the new hypersphere center, and the procedure is iterated
until the bottom of the energy surface is reached. Notice that in the TS case a pair of
steepest descent paths, connecting the TS to the reactant and product structures
(following the forward and reverse orientation of the direction defined by the transition
vector) provides the minimum energy path (MEP) for the reaction. As mass-weighted
coordinates are used by default, the MEP coordinate corresponds to the so-called Intrinsic
Reaction Coordinates (IRC). We shall compute here the MEP from the FC structure of acrolein
along the PES of the second root in energy at the CASSCF level. It is important to remember
that the CASSCF order may not be accurate and the states may reverse orders at higher
levels such as CASPT2.

.. extractfile:: problem_based_tutorials/CASSCF.mep_excited_state.acrolein.input

  *CASSCF excited state mep points in acrolein
  *File: CASSCF.mep_excited_state.acrolein
  *
   &GATEWAY
  Title = acrolein mep calculation root 2
  Basis set
  O.STO-3G...2s1p.
   O1    1.367073     0.000000     3.083333 Angstrom
  End of basis
  Basis set
  C.STO-3G...2s1p.
   C1    0.000000     0.000000     0.000000 Angstrom
   C2    0.000000     0.000000     1.350000 Angstrom
   C3    1.367073     0.000000     1.833333 Angstrom
  End of basis
  Basis set
  H.STO-3G...1s.
   H1    2.051552     0.000000     0.986333 Angstrom
   H2   -0.684479     0.000000     2.197000 Angstrom
   H3   -1.026719     0.000000    -0.363000 Angstrom
   H4    0.513360     0.889165    -0.363000 Angstrom
  End of basis

  >>> EXPORT MOLCAS_MAXITER=300
  >>> Do while

   &SEWARD
  >>> If ( Iter = 1 ) <<<
   &SCF
  >>> EndIf <<<

   &RASSCF
     Title="acrolein mep calculation root 2"; Spin=1
     nActEl=4 0 0; Inactive=13; Ras2=4; CiRoot=2 2 1; Rlxroot=2
   &SLAPAF
     MEP-search
     MEPStep=0.1

  >>> EndDo

As observed, to prepare the input for the MEP is simple, just add the keyword :kword:`MEP-search`
and specify a step size with :kword:`MEPStep`, and the remaining structure equals that of a geometry optimization.
The calculations are time consuming, because each point of the
MEP (four plus the initial one obtained here) is computed through a specific optimization.
A file named :file:`$Project.mep.molden` (read by :program:`MOLDEN` )
will be generated in $WorkDir containing only those points belonging to the MEP.

We shall now show how to perform geometry optimizations under nongeometrical
restrictions, in particular, how to compute hypersurface crossings, which are key structures
in the photophysics of molecules. We shall get those points as minimum energy crossing points in
which the energy of the highest of the two states considered is minimized under the restriction
that the energy difference with the lowest state should equal certain value (typically zero).
Such point can be named a minimum energy crossing point (MECP). If a further restriction is
imposed, like the distance to a specific geometry, and several MECP as computed at varying distances,
it is possible to obtain a crossing seam of points where the energy between the two states is
degenerated. Those degeneracy points are funnels with the highest probability for the energy
to hop between the surfaces in internal conversion or intersystem crossing photophysical processes.
There are different possibilities. A crossing between states of the same spin
multiplicity and spatial symmetry is named a conical intersection. Elements like the nonadiabatic
coupling terms are required to obtain them strictly, and they are not computed presently
by |molcas|. If the crossing occurs between states of the same
spin multiplicity and different spatial symmetry or between states of different spin multiplicity,
the crossing is an hyperplane and its only requirement is the energetic degeneracy and the
proper energy minimization.

Here we include an example with the crossing between the lowest singlet (ground) and triplet
states of acrolein. Notice that two different states are computed, first by using
:program:`RASSCF` to get the wave function and then :program:`ALASKA` to get the gradients
of the energy. Nothing new on that, just the information needed in any geometry optimizations.
The :program:`GATEWAY` input requires to add as constraint an energy
difference between both states equal to zero. A specific instruction is required after
calculating the first state. We have to copy the communication file :file:`RUNFILE`
(at that point contains the information about the first state) to :file:`RUNFILE2`
to provide later :program:`SLAPAF` with proper information about both states:

.. extractfile:: problem_based_tutorials/CASSCF.S-T_crossing.acrolein.input

  *CASSCF singlet-triplet crossing in acrolein
  *File: CASSCF.S-T_crossing.acrolein
  *
   &GATEWAY
  Title= Acrolein molecule
  Basis set
  O.sto-3g....
   O1             1.5686705444       -0.1354553340        3.1977912036  Angstrom
  End of basis
  Basis set
  C.sto-3g....
   C1            -0.1641585340        0.2420235062       -0.0459895824  Angstrom
   C2             0.1137722023       -0.1389623714        1.3481527296  Angstrom
   C3             1.3218729238        0.1965728073        1.9959513294  Angstrom
  End of basis
  Basis set
  H.sto-3g....
   H1             2.0526602523        0.7568282320        1.4351034056  Angstrom
   H2            -0.6138178851       -0.6941171027        1.9113821810  Angstrom
   H3            -0.8171509745        1.0643342316       -0.2648232855  Angstrom
   H4             0.1260134708       -0.4020589690       -0.8535699812  Angstrom
  End of basis
  Constraints
     a = Ediff
    Value
     a = 0.000
  End of Constraints

  >>> Do while

   &SEWARD

  >>> IF ( ITER = 1 ) <<<
   &SCF
  >>> ENDIF <<<

   &RASSCF
     LumOrb
     Spin= 1; Nactel= 4 0 0; Inactive= 13; Ras2= 4
     CiRoot= 1 1; 1
   &ALASKA
  >>COPY $WorkDir/$Project.RunFile $WorkDir/RUNFILE2

   &RASSCF
     LumOrb
     Spin= 3; Nactel= 4 0 0; Inactive= 13; Ras2= 4
     CiRoot= 1 1; 1
   &ALASKA
   &SLAPAF

  >>> EndDo

Solvent effects can be also applied to excited states, but first the reaction
field in the ground (initial) state has to be computed. This is because solvation in
electronic excited states is a non equilibrium situation in with the electronic
polarization effects (fast part of the reaction field) have to treated apart
(they supposedly change during the excitation process) from the orientational
(slow part) effects. The slow fraction of the reaction field is maintained from
the initial state and therefore a previous calculation is required.
From the practical point of view the input is simple as illustrated in the next
example. First, the proper reaction-field
input is included in :program:`SEWARD`, then a :program:`RASSCF` and :program:`CASPT2`
run of the ground state, with keyword :kword:`RFPErt` in :program:`CASPT2`,
and after that another SA-CASSCF calculation of five roots to get the wave function
of the excited states. Keyword :kword:`NONEequilibrium` tells the program to extract
the slow part of the reaction field from the previous calculation of the ground
state (specifically from the :file:`JOBOLD` file, which may be stored for other
calculations) while the fast part is freshly computed. Also, as it is a SA-CASSCF
calculation (if not, this is not required) keyword :kword:`RFRoot` is introduced
to specify for which of the computed roots the reaction field is generated. We have
selected here the fifth root because it has a very large dipole moment, which is
also very different from the ground state dipole moment. If you compare the excitation
energy obtained for the isolated and the solvated system, a the large red shift is
obtained in the later.

.. extractfile:: problem_based_tutorials/CASPT2.excited_solvent.acrolein.input

  *CASPT2 excited state in water for acrolein
  *File: CASPT2.excited_solvent.acrolein
  *
  &GATEWAY
    Title= Acrolein molecule
    coord = acrolein.xyz; basis = STO-3G; group= c1
    RF-input
     PCM-model; solvent= water
    End of RF-input
  &SEWARD
  &RASSCF
    Spin= 1; Nactel= 6 0 0; Inactive= 12; Ras2= 5
    CiRoot= 1 1 1
  &CASPT2
    Multistate= 1 1
    RFPert
  &RASSCF
    Spin= 1; Nactel= 6 0 0; Inactive= 12; Ras2= 5
    CiRoot= 5 5 1
    RFRoot= 5
    NONEquilibrium
  &CASPT2
    Multistate= 1 5
    RFPert

A number of simple examples as how to proceed with the most frequent
quantum chemical problems computed with |molcas| have been given above. Certainly there are many more
possibilities in |molcas|, such as calculation of 3D band
systems in solids at a semiempirical level, obtaining valence-bond structures,
the use of QM/MM methods in combination with a external MM code, the introduction
of external homogeneous or non homogeneous perturbations, generation of atomic
basis sets, application of different localization schemes, analysis of first
order polarizabilities, calculation of vibrational intensities, analysis, generation,
and fitting of potentials, computation of vibro-rotational spectra for diatomic
molecules, introduction of relativistic effects, etc. All those aspects are
explained in the manual and are much more specific. Next section :ref:`TUT:sec:pg-based-tut`
details the basic structure of the inputs, program by program, while easy examples
can also be found. Later, another chapter includes a number of extremely detailed
examples with more elaborated quantum chemical examples, in which also scientific
comments are included. Examples include calculations on high symmetry molecules,
geometry optimizations and Hessians, computing reaction paths, high quality wave
functions, excited states, solvent models, and computation of relativistic effects.
