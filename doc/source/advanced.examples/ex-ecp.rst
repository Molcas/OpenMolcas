.. index::
   single: ECP
   single: Pseudopotentials
   single: Embedding potentials
   single: SEWARD; ECP
   single: SEWARD; Pseudopotentials
   single: SEWARD; Embedding potentials
   single: SEWARD; Core potentials
   single: Core; Core potentials
   single: Integrals; Core potentials

.. _TUT\:sec\:ecp:

Core and Embedding Potentials within the :program:`SEWARD` Program
==================================================================

.. only:: html

  .. contents::
     :local:
     :backlinks: none

|molcas| is able to perform *effective core potential* (ECP)
and *embedded cluster* (EC) calculations.
In ECP calculations :cite:`Wahlgren:92,Seijo:99`
the *core* electrons of a molecule are kept frozen and represented by a set of atomic
effective potentials, while only the valence electrons are explicitly handled
in the quantum mechanical calculation. In EC calculations only the electrons
assigned to a piece of the whole system, the *cluster*, are explicitly
treated in a quantum mechanical calculation, while the rest of the whole
system, the *environment*, is kept frozen and represented by embedding
potentials which act onto the *cluster*. For an explanation of the
type of potentials and approaches used in |molcas| the reader is referred
to the section :ref:`UG:sec:the_ecp_libraries` of the user's guide.

To use such type of effective potentials implies to compute a set
of atomic integrals and therefore involves only the :program:`SEWARD` program.
The remaining |molcas| programs will simply use the integrals in the
standard way and no indication of the use of ECP will appear in the
outputs further on; the difference is of course that the absolute energies
obtained for the different methods are not comparable to those obtained
in an all-electron calculation. Therefore, the only input required to
use ECP or EC is the :program:`SEWARD` input, according to the examples
given below. In the input files of the subsequent |molcas| programs the
orbitals corresponding to the excluded core orbitals should of course not be
included, and not the excluded electrons.

:program:`seward` input for Effective Core Potential calculations
-----------------------------------------------------------------

Astatine (:math:`\ce{At}`) is the atomic element number 85 which has the main configuration
in its electronic ground state: [*core*] 6s\ :math:`^2` 5d\ :math:`^{10}` 6p\ :math:`^5`. In the
*core* 68 electrons are included, corresponding to the xenon configuration
plus the 4f\ :math:`^{14}` lantanide shell. To perform an ECP calculation in a
molecular system containing :math:`\ce{At}` it is necessary to specify which type of
effective potential will substitute the *core* electrons and which valence
basis set will complement it. Although the core ECP's (strictly AIMP's, see
section :ref:`UG:sec:the_ecp_libraries` of the user's guide) can be safely
mixed together with all-electron basis set, the valence basis sets included
in the |molcas| AIMP library have been explicitly optimized to complement the
AIMP potentials.

.. index::
   single: Relativistic effect; Core potentials

The file :file:`ECP` in the |molcas| directory :file:`$MOLCAS/basis_library` contains the
list of available core potentials and valence basis sets. Both the relativistic
(CG-AIMP's) and the nonrelativistic (NR-AIMP's) potentials are included. As an
example, this is the head of the entry corresponding to the relativistic ECP
for :math:`\ce{At}`: ::

  /At.ECP.Barandiaran.13s12p8d5f.1s1p2d1f.17e-CG-AIMP.
  Z.Barandiaran, L.Seijo, J.Chem.Phys. 101(1994)4049; L.S. JCP 102(1995)8078.
  core[Xe,4f] val[5d,6s,6p]  SO-corr  (11,1,1/9111/611*/4o1)=3s4p3d2f recommended
  *
  * - spin-orbit basis set correction from
  *   L.Seijo, JCP 102(1995)8078.
  *
  * - (5o) f orthogonality function is the 4f core orbital
  *
  *ATQR-DSP(A3/A2/71/5)-SO       (A111/9111/611/41)

The first line is the label line written in the usual :program:`SEWARD` format:
element symbol, basis label, first author, size of the primitive set, size
of the contracted set (in both cases referred to the valence basis set), and
type of ECP used. In this case there are 17 valence electrons and the
effective potential is a Cowan--Griffin-relativistic core AIMP. The number of
primitive functions for the valence basis set (13s12p8d5f here) will split
into different subsets (within a segmented contraction scheme) according to
the number of contracted functions. In the library, the contracted
basis functions have been set to the minimal basis size: 1s1p2d1f for the
valence electrons in :math:`\ce{At}`. This means the following partition: 1s contracted
function including 13 primitive functions; 1p contracted function including
12 primitive functions; 2d contracted functions, the first one containing
seven primitive functions and the second one primitive function
(see the library), and finally 1f contracted function containing five
primitive functions.

In the :program:`SEWARD` input the user can modify the contraction scheme
simply varying the number of contracted functions. There is a recommended size
for the valence basis set which is printed in the third line for each atom entry
on the library: 3s4p3d2f for :math:`\ce{At}`. For example, the simplest way to include the
atom core potential and valence basis set in the :program:`SEWARD` input would
be: ::

  At.ECP...3s4p3d2f.17e-CG-AIMP.

This means a partition for the valence basis set as showed in
:numref:`block:valbas_ecp`.

.. code-block:: none
   :caption: Partition of a valence basis set using the ECP's library
   :name: block:valbas_ecp

   Basis set:AT.ECP...3S4P3D2F.17E-CG-AIMP.

                    Type
                     s
             No.      Exponent    Contraction Coefficients
              1   .133037396D+07  -.000154   .000000   .000000
              2   .993126141D+05  -.001030   .000000   .000000
              3   .128814005D+05  -.005278   .000000   .000000
              4   .247485916D+04  -.014124   .000000   .000000
              5   .214733934D+03   .069168   .000000   .000000
              6   .111579706D+03   .020375   .000000   .000000
              7   .370830653D+02  -.259246   .000000   .000000
              8   .113961072D+02   .055751   .000000   .000000
              9   .709430236D+01   .649870   .000000   .000000
             10   .448517638D+01  -.204733   .000000   .000000
             11   .157439587D+01  -.924035   .000000   .000000
             12   .276339384D+00   .000000  1.000000   .000000
             13   .108928284D+00   .000000   .000000  1.000000

                    Type
                     p
             No.      Exponent    Contraction Coefficients
             14   .608157825D+04   .000747   .000000   .000000   .000000
             15   .128559298D+04   .009304   .000000   .000000   .000000
             16   .377428675D+03   .026201   .000000   .000000   .000000
             17   .552551834D+02  -.087130   .000000   .000000   .000000
             18   .233740022D+02  -.044778   .000000   .000000   .000000
             19   .152762905D+02   .108761   .000000   .000000   .000000
             20   .838467359D+01   .167650   .000000   .000000   .000000
             21   .234820847D+01  -.290968   .000000   .000000   .000000
             22   .119926577D+01  -.237719   .000000   .000000   .000000
             23   .389521915D+00   .000000  1.000000   .000000   .000000
             24   .170352883D+00   .000000   .000000  1.000000   .000000
             25   .680660800D-01   .000000   .000000   .000000  1.000000

                    Type
                     d
             No.      Exponent    Contraction Coefficients
             26   .782389711D+03   .007926   .000000   .000000
             27   .225872717D+03   .048785   .000000   .000000
             28   .821302011D+02   .109617   .000000   .000000
             29   .173902999D+02  -.139021   .000000   .000000
             30   .104111329D+02  -.241043   .000000   .000000
             31   .195037661D+01   .646388   .000000   .000000
             32   .689437556D+00   .000000  1.000000   .000000
             33   .225000000D+00   .000000   .000000  1.000000

                    Type
                     f
             No.      Exponent    Contraction Coefficients
             34   .115100000D+03   .065463   .000000
             35   .383200000D+02   .270118   .000000
             36   .151600000D+02   .468472   .000000
             37   .622900000D+01   .387073   .000000
             38   .242100000D+01   .000000  1.000000

Therefore, the primitive set will always be split following the scheme:
the first contracted function will contain the total number of primitives
minus the number of remaining contracted functions and each of the
remaining contracted functions will contain one single uncontracted
primitive function. In the present example possible contraction patterns
are: contracted 1s1p2d1f (13/12/8,1/5 primitives per contracted function, respectively),
2s2p3d2f (12,1/11,1/7,1,1/4,1), 3s3p4d2f (11,1,1/10,1,1/6,1,1,1/4,1), etc.
Any other scheme which cannot be generated in this way must be included in
the input using the Inline format for basis sets or an additional user's library.
When the Inline option is
used both the valence basis set and the AIMP potential must be included in
the input, as it will be shown in the next section.

For an explanation of the remaining items in the library the reader is referred
to the section :ref:`UG:sec:the_ecp_libraries` of the user's guide.

:numref:`block:hat_scf` contains the sample input required to compute the
SCF wave function for the astatine hydride molecule at an internuclear
distance of 3.2 au.
The Cowan--Griffin-relativistic core-AIMP has been
used for the :math:`\ce{At}` atom with a size for the valence basis set recommended in the
:file:`ECP` library: 3s4p3d2f.

.. extractcode-block:: none
   :filename: advanced/ECP.HAt.input
   :caption: Sample input required by SEWARD and SCF programs to compute the SCF
             wave function of :math:`\ce{HAt}` using a relativistic ECP
   :name: block:hat_scf

   &GATEWAY
   Title
   HAt molecule using 17e-Cowan-Griffin-relativistic core-AIMP
   coord
   2
   coordinates in bohr
   At 0 0 0
   H  0 0 3.2
   group
   X Y
   Basis set
   H.ano-l-vtzp
   Basis set
   At.ECP...3s4p3d2f.17e-CG-AIMP.
   &SEWARD
   &SCF
   Title
    HAt g.s. (At-val=5d,6s,6p)
   Occupied
    4 2 2 1

.. index::
   single: Embedded clusters
   single: Lattice

:program:`seward` input for Embedded Cluster calculations
---------------------------------------------------------

To perform embedded cluster (EC) calculations requires certain degree
of experience and therefore the reader is referred to the literature
quoted in section :ref:`UG:sec:the_ecp_libraries` of the user's guide.
On the following a detailed example is however presented.
It corresponds to EC calculations useful for local properties
associated to a :math:`\ce{Tl^+}` impurity in :math:`\ce{KMgF3}`. First, a cluster must be
specified. This is the piece of the system which is explicitly treated by the
quantum mechanical calculation. In the present example the cluster will be
formed by the unit :math:`\ce{(TlF_{12})^{11-}}`. A flexible basis for the cluster must be
determined. :numref:`block:tlf_input` contains the basis set selection
for the thallium and fluorine atoms. In this case ECP-type basis sets
have been selected. For :math:`\ce{Tl}` a valence basis set of size 3s4p4d2f has
been used combined with the relativistic core-AIMP potentials as they
appear in the :file:`ECP` library. For the :math:`\ce{F}` atom the valence
basis set has been modified from that appearing in the :file:`ECP`
library. In this case the exponent of the p-diffuse function and the p
contraction coefficients
of the :math:`\ce{F}` basis set have been optimized in calculations on the fluorine
anion included in the specific lattice in order to obtain a more
flexible description of the anion. This
basis set must be introduced Inline, and then also the ECP potential
must be added to the input. The user can compare the basis set
and ECP for :math:`\ce{F}` in :numref:`block:tlf_input` with the entry of :file:`ECP`
under /F.ECP.Huzinaga.5s6p1d.1s2p1d.7e-NR-AIMP. The entry for the
Inline format must finish with the line End of Spectral Representation Operator.

Once the cluster has been defined it is necessary to represent the embedding
lattice. Presently, |molcas| includes embedding potentials for ions of
several elpasolites, fluoro-perovskites, rocksalt structure oxides and halides,
and fluorites. The embedding potentials for any other structure can be included
in the input using the Inline format
or included in a private user library.
In the selected example a fluoro-perovskite lattice has
been selected: :math:`\ce{KMgF3}`.
Here, the :math:`\ce{Tl^+}` impurity substitutes a :math:`\ce{K^+}` ion in an :math:`O_h` site with
12 coordination.
The first coordination shell of fluorine ions has been included into the cluster
structure and the interactions to the :math:`\ce{Tl}` atom will be computed by quantum
mechanical methods. The rest of the lattice will be represented by the
structure :math:`\ce{KMgF3}` with five shells of ions at experimental sites.
The shells have been divided in two types. Those shells closer to the
cluster are included as embedding potentials from the library :file:`ECP`.
For example the potassium centers will use the entry on :numref:`block:tlf_k`.

.. code-block:: none
   :caption: Sample input for an embedded core potential for a shell of potassium cations
   :name: block:tlf_k

   Basis set
   K.ECP..0s.0s.0e-AIMP-KMgF3.
   PSEUdocharge
   K2-1    0.0000000000   0.0000000000   7.5078420000
   K2-2    0.0000000000   7.5078420000   0.0000000000
   K2-3    0.0000000000   7.5078420000   7.5078420000
   K2-4    7.5078420000   0.0000000000   0.0000000000
   K2-5    7.5078420000   0.0000000000   7.5078420000
   K2-6    7.5078420000   7.5078420000   0.0000000000
   K2-7    7.5078420000   7.5078420000   7.5078420000
   End Of Basis

No basis set is employed to represent the potassium centers on :numref:`block:tlf_k`,
which just act as potentials embedding the cluster. The keyword
:kword:`PSEUdocharge` ensures that the interaction energy between the embedding
potentials is not included in the "Nuclear repulsion energy"
and that their location is not varied in a geometry optimization (:program:`SLAPAF`).
The first shells of :math:`\ce{Mg^{+2}}` and :math:`\ce{F^-}` will be introduced in the same way.

The remaining ions of the lattice will be treated as point charges.
To add a point charge on the :program:`SEWARD` input it is possible to proceed
in two ways. One possibility is to employ the usual label to introduce an atom
with its basis functions set to zero and the keyword :kword:`CHARge` set to the
value desired for the charge of the center. This way of introducing point charges must not be
used when geometry optimizations with the :program:`SLAPAF` program is going to
be performed because :program:`SLAPAF` will recognize the point charges as atoms
whose positions should be optimized. Instead the keyword :kword:`XFIEld` can be
used as it is illustrated in :numref:`block:tlf_input`. :kword:`XFIEld` must
be followed by a line containing the number of point charges, and by subsequent
lines containing the cartesian coordinates and the introduced charge or the
three components of the dipole moment at the specified geometry. In any case
the seven positions in each line must be fulfilled. To ensure the neutral
character of the whole system the point charges placed on the terminal edges,
corners or faces of the lattice must have the proper fractional values.

:numref:`block:tlf_input` contains the complete sample input to perform a
SCF energy calculation on the system :math:`\ce{(TlF_{12})^{11-}{:}KMgF3}`.

.. extractcode-block:: none
   :filename: advanced/ECP.TlF12.input
   :caption: Sample input for a SCF geometry optimization of the :math:`\ce{(TlF_{12})^{11-}{:}KMgF3}` system
   :name: block:tlf_input

   &GATEWAY
   Title
   |                          Test run TlF12:KMgF3.1                              |
   |** Molecule **   (TlF12)11- cluster embedded in a lattice of KMgF3            |
   |** Basis set and ECP **                                                       |
   |  * Tl * (11,1,1/9,1,1,1/5,1,1,1/4,1)                             from ECP    |
   |         13e-Cowan-Griffin-relativistic core-AIMP                 from ECP    |
   |  * F *  (4,1/4,1,1) diffuse-p optimized in KMgF3:F(-)                  inline|
   |          7e-nonrelativistic core-AIMP                                  inline|
   |  KMgF3 embedding-AIMPs                                           from ECP    |
   |** cluster geometry **   r(Tl-F)/b= 5.444 = 3.84948932 * sqrt(2)              |
   |** lattice **  (perovskite structure) 5 shells of ions at experimental sites  |
   Symmetry
   X Y Z

   Basis set
   Tl.ECP.Barandiaran.13s12p8d5f.3s4p4d2f.13e-CG-AIMP.
   Tl     0.00000   0.00000   0.00000
   End Of Basis

   Basis set
   F.ECP.... / Inline
   *    basis set and core-AIMP as in: F.ECP.Huzinaga.5s6p1d.2s4p1d.7e-NR-AIMP.
   *    except that the p-diffuse and the p contraction coeffs. have been
   *    optimized in KMgF3-embedded F(-) scf calculations.
     7.000000         1
       5    2
      405.4771610
      61.23686380
      13.47117730
      1.095173720
      .3400847530
     -.013805187800   .000000000000
     -.089245064800   .000000000000
     -.247937861000   .000000000000
      .632895340000   .000000000000
      .000000000000   .465026336000
       6    3
      44.13600920
      9.982597110
      2.947082680
      .9185111850
      .2685213550
      .142
      .015323038700   .000000000000   .000000000000
      .095384703000   .000000000000   .000000000000
      .291214218000   .000000000000   .000000000000
      .441351868000   .000000000000   .000000000000
      .000000000000   .427012588000   .000000000000
      .000000000000   .000000000000  1.000000000000
   *
   * Core AIMP: F-1S
   *
   * Local Potential Paramenters : (ECP convention)
   *                               A(AIMP)=-Zeff*A(ECP)
   M1
       7
      279347.4000
      31889.74900
      5649.977600
      1169.273000
      269.0513200
      71.29884600
      22.12150700

      .004654725000
      .007196816857
      .015371258571
      .032771900000
      .070383742857
      .108683807143
      .046652035714
   M2
       0
   COREREP
      1.0
   PROJOP
       0
      14    1
     52.7654040
      210965.4100
      31872.59200
      7315.837400
      2077.215300
      669.9991000
      232.1363900
      84.99573000
      32.90124100
      13.36331800
      5.588141500
      2.319058700
      .9500928100
      .3825419200
      .1478404000
      .000025861368
      .000198149380
      .001031418900
      .004341016600
      .016073698000
      .053856655000
      .151324390000
      .318558040000
      .404070310000
      .190635320000
      .011728993000
      .002954046500
     -.000536098280
      .000278474090
   *
   Spectral Representation Operator
   Valence primitive basis
   Exchange
   End of Spectral Representation Operator
   F_1        3.849489320       3.849489320        .000000000
   F_2         .000000000       3.849489320       3.849489320
   F_3        3.849489320        .000000000       3.849489320
   * 3*4 = 12
   End Of Basis

   * end of cluster data: TlF12

   * beginning of lattice embedding data: KMgF3

   Basis set
   K.ECP.Lopez-Moraza.0s.0s.0e-AIMP-KMgF3.
   pseudocharge
   * K(+) ions as embedding AIMPs
   K2-1    0.0000000000   0.0000000000   7.5078420000
   K2-2    0.0000000000   7.5078420000   0.0000000000
   K2-3    0.0000000000   7.5078420000   7.5078420000
   K2-4    7.5078420000   0.0000000000   0.0000000000
   K2-5    7.5078420000   0.0000000000   7.5078420000
   K2-6    7.5078420000   7.5078420000   0.0000000000
   K2-7    7.5078420000   7.5078420000   7.5078420000
   * 3*2 + 3*4 + 1*8 = 26
   End Of Basis

   Basis set
   Mg.ECP.Lopez-Moraza.0s.0s.0e-AIMP-KMgF3.
   pseudocharge
   * Mg(2+) ions as embedding AIMPs
   MG1-1   3.7539210000   3.7539210000   3.7539210000
   MG3-1   3.7539210000   3.7539210000  11.2617630000
   MG3-2   3.7539210000  11.2617630000   3.7539210000
   MG3-3   3.7539210000  11.2617630000  11.2617630000
   MG3-4  11.2617630000   3.7539210000   3.7539210000
   MG3-5  11.2617630000   3.7539210000  11.2617630000
   MG3-6  11.2617630000  11.2617630000   3.7539210000
   MG3-7  11.2617630000  11.2617630000  11.2617630000
   * 8*8 = 64
   End Of Basis

   Basis set
   F.ECP.Lopez-Moraza.0s.0s.0e-AIMP-KMgF3.
   pseudocharge
   * F(-) ions as embedding AIMPs
   F2-1    3.7539210000   3.7539210000   7.5078420000
   F2-2    3.7539210000   7.5078420000   3.7539210000
   F2-3    7.5078420000   3.7539210000   3.7539210000
   F3-1    0.0000000000   3.7539210000  11.2617630000
   F3-2    3.7539210000   0.0000000000  11.2617630000
   F3-3    3.7539210000  11.2617630000   0.0000000000
   F3-4    0.0000000000  11.2617630000   3.7539210000
   F3-5    3.7539210000  11.2617630000   7.5078420000
   F3-6    0.0000000000  11.2617630000  11.2617630000
   F3-7    3.7539210000   7.5078420000  11.2617630000
   F3-8   11.2617630000   3.7539210000   0.0000000000
   F3-9   11.2617630000   0.0000000000   3.7539210000
   F3-10   11.2617630000   3.7539210000   7.5078420000
   F3-11    7.5078420000   3.7539210000  11.2617630000
   F3-12   11.2617630000   0.0000000000  11.2617630000
   F3-13   11.2617630000  11.2617630000   0.0000000000
   F3-14    7.5078420000  11.2617630000   3.7539210000
   F3-15   11.2617630000   7.5078420000   3.7539210000
   F3-16   11.2617630000  11.2617630000   7.5078420000
   F3-17    7.5078420000  11.2617630000  11.2617630000
   F3-18   11.2617630000   7.5078420000  11.2617630000
   * 9*4 +  12*8 = 132
   End Of Basis

   * The rest of the embedding lattice will be represented by point charges,
   * which enter into the calculation in the form of a XField.
   *
   XField
    95
   *
   * K(+) ions as point charges
       0.0000000000   0.0000000000  15.0156840000       +1.0  0.  0.  0.
       0.0000000000   7.5078420000  15.0156840000       +1.0  0.  0.  0.
       0.0000000000  15.0156840000   0.0000000000       +1.0  0.  0.  0.
       0.0000000000  15.0156840000   7.5078420000       +1.0  0.  0.  0.
       0.0000000000  15.0156840000  15.0156840000       +1.0  0.  0.  0.
       7.5078420000   0.0000000000  15.0156840000       +1.0  0.  0.  0.
       7.5078420000   7.5078420000  15.0156840000       +1.0  0.  0.  0.
       7.5078420000  15.0156840000   0.0000000000       +1.0  0.  0.  0.
       7.5078420000  15.0156840000   7.5078420000       +1.0  0.  0.  0.
       7.5078420000  15.0156840000  15.0156840000       +1.0  0.  0.  0.
      15.0156840000   0.0000000000   0.0000000000       +1.0  0.  0.  0.
      15.0156840000   0.0000000000   7.5078420000       +1.0  0.  0.  0.
      15.0156840000   0.0000000000  15.0156840000       +1.0  0.  0.  0.
      15.0156840000   7.5078420000   0.0000000000       +1.0  0.  0.  0.
      15.0156840000   7.5078420000   7.5078420000       +1.0  0.  0.  0.
      15.0156840000   7.5078420000  15.0156840000       +1.0  0.  0.  0.
      15.0156840000  15.0156840000   0.0000000000       +1.0  0.  0.  0.
      15.0156840000  15.0156840000   7.5078420000       +1.0  0.  0.  0.
      15.0156840000  15.0156840000  15.0156840000       +1.0  0.  0.  0.
   *
   * F(-) ions as point charges
       3.7539210000   3.7539210000  15.0156840000       -1.0  0.  0.  0.
       3.7539210000  11.2617630000  15.0156840000       -1.0  0.  0.  0.
       3.7539210000  15.0156840000   3.7539210000       -1.0  0.  0.  0.
       3.7539210000  15.0156840000  11.2617630000       -1.0  0.  0.  0.
      11.2617630000   3.7539210000  15.0156840000       -1.0  0.  0.  0.
      11.2617630000  11.2617630000  15.0156840000       -1.0  0.  0.  0.
      11.2617630000  15.0156840000   3.7539210000       -1.0  0.  0.  0.
      11.2617630000  15.0156840000  11.2617630000       -1.0  0.  0.  0.
      15.0156840000   3.7539210000   3.7539210000       -1.0  0.  0.  0.
      15.0156840000   3.7539210000  11.2617630000       -1.0  0.  0.  0.
      15.0156840000  11.2617630000   3.7539210000       -1.0  0.  0.  0.
      15.0156840000  11.2617630000  11.2617630000       -1.0  0.  0.  0.
   *
   * Mg(2+) ions in face, as fractional point charges
      3.7539210000   3.7539210000  18.7696050000        +1.0  0.  0.  0.
      3.7539210000  11.2617630000  18.7696050000        +1.0  0.  0.  0.
      3.7539210000  18.7696050000   3.7539210000        +1.0  0.  0.  0.
      3.7539210000  18.7696050000  11.2617630000        +1.0  0.  0.  0.
     11.2617630000   3.7539210000  18.7696050000        +1.0  0.  0.  0.
     11.2617630000  11.2617630000  18.7696050000        +1.0  0.  0.  0.
     11.2617630000  18.7696050000   3.7539210000        +1.0  0.  0.  0.
     11.2617630000  18.7696050000  11.2617630000        +1.0  0.  0.  0.
     18.7696050000   3.7539210000   3.7539210000        +1.0  0.  0.  0.
     18.7696050000   3.7539210000  11.2617630000        +1.0  0.  0.  0.
     18.7696050000  11.2617630000   3.7539210000        +1.0  0.  0.  0.
     18.7696050000  11.2617630000  11.2617630000        +1.0  0.  0.  0.
   *
   * Mg(2+) ions in edge, as fractional point charges
      3.7539210000  18.7696050000  18.7696050000     +0.5  0.  0.  0.
     11.2617630000  18.7696050000  18.7696050000     +0.5  0.  0.  0.
     18.7696050000   3.7539210000  18.7696050000     +0.5  0.  0.  0.
     18.7696050000  11.2617630000  18.7696050000     +0.5  0.  0.  0.
     18.7696050000  18.7696050000   3.7539210000     +0.5  0.  0.  0.
     18.7696050000  18.7696050000  11.2617630000     +0.5  0.  0.  0.
   *
   * Mg(2+) ions in corner, as fractional point charges
     18.7696050000  18.7696050000  18.7696050000      +0.25  0. 0. 0.
   *
   * F(-) ions in face, as fractional point charges
      0.0000000000   3.7539210000  18.7696050000       -0.5  0. 0. 0.
      3.7539210000   0.0000000000  18.7696050000       -0.5  0. 0. 0.
      0.0000000000  11.2617630000  18.7696050000       -0.5  0. 0. 0.
      3.7539210000   7.5078420000  18.7696050000       -0.5  0. 0. 0.
      3.7539210000  18.7696050000   0.0000000000       -0.5  0. 0. 0.
      0.0000000000  18.7696050000   3.7539210000       -0.5  0. 0. 0.
      3.7539210000  18.7696050000   7.5078420000       -0.5  0. 0. 0.
      0.0000000000  18.7696050000  11.2617630000       -0.5  0. 0. 0.
      3.7539210000  18.7696050000  15.0156840000       -0.5  0. 0. 0.
      3.7539210000  15.0156840000  18.7696050000       -0.5  0. 0. 0.
      7.5078420000   3.7539210000  18.7696050000       -0.5  0. 0. 0.
     11.2617630000   0.0000000000  18.7696050000       -0.5  0. 0. 0.
      7.5078420000  11.2617630000  18.7696050000       -0.5  0. 0. 0.
     11.2617630000   7.5078420000  18.7696050000       -0.5  0. 0. 0.
     11.2617630000  18.7696050000   0.0000000000       -0.5  0. 0. 0.
      7.5078420000  18.7696050000   3.7539210000       -0.5  0. 0. 0.
     11.2617630000  18.7696050000   7.5078420000       -0.5  0. 0. 0.
      7.5078420000  18.7696050000  11.2617630000       -0.5  0. 0. 0.
     11.2617630000  18.7696050000  15.0156840000       -0.5  0. 0. 0.
     11.2617630000  15.0156840000  18.7696050000       -0.5  0. 0. 0.
     18.7696050000   3.7539210000   0.0000000000       -0.5  0. 0. 0.
     18.7696050000   0.0000000000   3.7539210000       -0.5  0. 0. 0.
     18.7696050000   3.7539210000   7.5078420000       -0.5  0. 0. 0.
     18.7696050000   0.0000000000  11.2617630000       -0.5  0. 0. 0.
     18.7696050000   3.7539210000  15.0156840000       -0.5  0. 0. 0.
     15.0156840000   3.7539210000  18.7696050000       -0.5  0. 0. 0.
     18.7696050000  11.2617630000   0.0000000000       -0.5  0. 0. 0.
     18.7696050000   7.5078420000   3.7539210000       -0.5  0. 0. 0.
     18.7696050000  11.2617630000   7.5078420000       -0.5  0. 0. 0.
     18.7696050000   7.5078420000  11.2617630000       -0.5  0. 0. 0.
     18.7696050000  11.2617630000  15.0156840000       -0.5  0. 0. 0.
     15.0156840000  11.2617630000  18.7696050000       -0.5  0. 0. 0.
     15.0156840000  18.7696050000   3.7539210000       -0.5  0. 0. 0.
     18.7696050000  15.0156840000   3.7539210000       -0.5  0. 0. 0.
     15.0156840000  18.7696050000  11.2617630000       -0.5  0. 0. 0.
     18.7696050000  15.0156840000  11.2617630000       -0.5  0. 0. 0.
   *
   * F(-) ions in edge, as fractional point charges
      0.0000000000  18.7696050000  18.7696050000       -0.25  0. 0. 0.
      7.5078420000  18.7696050000  18.7696050000       -0.25  0. 0. 0.
     18.7696050000   0.0000000000  18.7696050000       -0.25  0. 0. 0.
     18.7696050000   7.5078420000  18.7696050000       -0.25  0. 0. 0.
     18.7696050000  18.7696050000   0.0000000000       -0.25  0. 0. 0.
     18.7696050000  18.7696050000   7.5078420000       -0.25  0. 0. 0.
     18.7696050000  18.7696050000  15.0156840000       -0.25  0. 0. 0.
     15.0156840000  18.7696050000  18.7696050000       -0.25  0. 0. 0.
     18.7696050000  15.0156840000  18.7696050000       -0.25  0. 0. 0.

   *  end of lattice embedding data: KMgF3

   * 13 cluster components  and 881 lattice components

   &SEWARD
   &SCF
   Title
    (TlF12)11- run as D2h
   Occupied
    12    7    7    6    7    6    6    3
