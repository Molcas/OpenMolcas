.. index::
   single: Program; Poly_aniso
   single: Poly_aniso

.. _UG\:sec\:poly_aniso:

:program:`poly_aniso`
=====================

.. only:: html

  .. contents::
     :local:
     :backlinks: none

.. xmldoc:: <MODULE NAME="POLY_ANISO" APPEAR="Poly_Aniso">
            %%Description:
            <HELP>
            The POLY_ANISO program allows the non-perturbative calculation of
            effective spin (pseudospin) Hamiltonians and static magnetic properties
            of mononuclear complexes and fragments completely ab initio,including
            the spin-orbit interaction. As a starting point it uses the results
            of a RASSI calculation for the ground and several excited spin-orbital
            multiplets.
            The following quantities can be computed:
            ||
            ||1. Parameters of pseudospin magnetic Hamiltonians:
            ||   a) First order (linear after pseudospin) Zeeman splitting tensor (g tensor),
            ||      including the determination of the sign of the product gX*gY*gZ
            ||   b) Second order (bilinear after pseudospin) zero-field splitting tensor (D tensor)
            ||   c) Higher order zero-field splitting tensors (D^2, D^4, D^6, ..., etc.)
            ||   d) Higher order Zeeman splitting tensors (G^1, G^3, G^5, ..., etc.)
            ||   e) Angular Moments along the main magnetic axes
            ||
            ||2. Crystal-Field parameters for the ground atomic multiplet for lanthanides
            ||
            ||3. Static magnetic properties:
            ||   a) Van Vleck susceptibility tensor
            ||   b) Powder magnetic susceptibility function
            ||   c) Magnetization vector for specified directions of the applied magnetic field
            ||   d) Powder magnetization
            </HELP>

The :program:`POLY_ANISO` program is a routine which allows a semi-ab initio
description of the (low-lying) electronic structure and magnetic properties
of polynuclear compounds. It is based on the localized nature of the magnetic
orbitals (i.e. the |d.| or |f.| orbitals containing unpaired electrons :cite:`Anderson1959,Anderson1963`).
For many compounds of interest, the localized character of magnetic orbitals leads
to very weak character of the interactions between magnetic centers. Due to this weakness of the
interaction, the metals' orbitals and corresponding localized ground and excited states
may be optimized in the absence of the magnetic interaction at all. For this purpose, various fragmentation
models may be applied. The most commonly used fragmentation model is exemplified in :numref:`fig:fragment`.

.. figure:: fragment.*
   :name: fig:fragment
   :width: 100%
   :align: center

   Fragmentation model of a polynuclear compound. The upper scheme shows a schematic overview of a tetranuclear compound and the resulting four mononuclear fragments obtained by *diamagnetic atom substitution* method. By this scheme, the neighboring magnetic centers, containing unpaired electrons are computationally replaced by their diamagnetic equivalents. As example, transition metal sites TM(II) are best replaced by either diamagnetic :math:`\ce{Zn(II)}` or :math:`\ce{Sc(III)}`, in function which one is the closest. For lanthanides :math:`\ce{Ln(III)}` the same principle is applicable, :math:`\ce{La(III)}` or :math:`\ce{Lu(III)}` are best suited to replace a given magnetic lanthanide. Individual mononuclear metal framgents are then investigated by common CASSCF/CASPT2/RASSI/SINGLE_ANISO computational method. A single file for each magnetic site, produced by the :program:`SINGLE_ANISO` run, is needed by the :program:`POLY_ANISO` code as input.

Magnetic interaction between metal sites is very important for accurate description of low-lying states and their properties.
It can be considered as a sum of various interaction mechanisms: magnetic exchange, dipole-dipole interaction, antisymmetric exchange, etc.
In the :program:`POLY_ANISO` code we have implemented several mechanisms.

The description of the magnetic exchange interaction is done within the Lines model :cite:`Lines1971`.
This model is exact in three cases:

a) interaction between two isotropic spins (Heisenberg),
b) interaction between one Ising spin (only :math:`S_z` component) and one isotropic (i.e. usual) spin, and
c) interaction between two Ising spins.

In all other cases of interaction between magnetic sites with intermediate anisotropy, the Lines model represents an
approximation. However, it was succesfully applied for a wide variety of polynuclear compounds so far.

.. compound::

  In addition to the magnetic exchange, magnetic dipole-dipole interaction can be accounted exactly, by
  using the information about each metal site already computed *ab initio*. In the case of
  strongly anisotropic lanthanide compounds, the dipole-dipole interaction is usualy the dominant
  one. Dipolar magnetic coupling is one kind of long-range interaction between magnetic moments.
  For example, a system containing two magnetic dipoles :math:`\mu_1` and :math:`\mu_2`, separated by
  distance :math:`\vec{r}` have a total energy:

  .. math:: E_{\text{dip}} = \frac{\mu_{\text{Bohr}}^{2}}{r^3} [\vec{\mu}_1 \cdot \vec{\mu}_2 - 3(\vec{\mu}_1 \vec{n}_{12}) \cdot (\vec{\mu}_2 \vec{n}_{12})],

  where :math:`\vec{\mu}_{1,2}` are the magnetic moments of sites 1 and 2, respectively; :math:`r` is the distance between
  the two magnetic dipoles, :math:`\vec{n}_{12}` is the directional vector connecting the two magnetic dipoles (of unit length).
  :math:`\mu_{\text{Bohr}}^2` is the square of the Bohr magneton; with an approximative value of 0.43297 in :math:`\text{cm}^{-1}/\text{T}`.
  As inferred from the above Equation, the dipolar magnetic interaction depends on the distance and on the angle between the magnetic moments on magnetic
  centers. Therefore, the cartesian coordinates of all non-equivalent magnetic centers must be provided in the input (see the keyword :kword:`COOR`).

Files
-----

Input files
...........

The program :program:`Poly_Aniso` needs the following files:

.. class:: filelist

:file:`aniso_XX.input`
  This is an ASCII text file generated by the |molcas|/SINGLE_ANISO program.
  It should be provided for :program:`POLY_ANISO` :file:`aniso_i.input` (:math:`i=1, 2, 3`, etc.): one file for each magnetic center.
  In cases when the entire polynuclear cluster or molecule has exact point group symmetry, only
  :file:`aniso_i.input` files for crystallographically non-equivalent centers should be given.

:file:`chitexp.input`
  set directly in the standard input (key :kword:`TEXP`)

:file:`magnexp.input`
  set directly in the standard input (key :kword:`HEXP`)

Output files
............

.. class:: filelist

:file:`zeeman_energy_xxx.txt`
  A series of files named :file:`zeeman_energy_xxx.txt` is produced in the :file:`$WorkDir` only in case keyword :kword:`ZEEM` is
  employed (see below). Each file is an ASCII text formated and contains Zeeman spectra of the investigated
  compound for each value of the applied magnetic field.

:file:`chit_compare.txt`
  A text file contining the experimental and calculated magnetic susceptibility data.

:file:`magn_compare.txt`
  A text file contining the experimental and calculated powder magnetisation data.

Files :file:`chit_compare.txt` and :file:`chit_compare.txt` may be used in connection with a simple gnuplot script
in order to plot the comparison between experimental and calculated data.

.. index::
   pair: Input; Poly_aniso

.. _UG\:sec\:poly_aniso_input:

Input
-----

This section describes the keywords used to control the standard input file.
Only two keywords :kword:`NNEQ`, :kword:`PAIR` (and :kword:`SYMM` if the polynuclear cluster has symmetry) are
mandatory for a minimal execution of the program, while the other keywords allow
customization of the execution of the :program:`POLY_ANISO`.

Mandatory keywords defining the calculation
...........................................

*Keywords defining the polynuclear cluster*

.. class:: keywordlist

:kword:`NNEQ`
  This keyword defines several important parameters of the calculation. On the
  first line after the keyword the program reads 2 values:
  1) the number of types of different magnetic centers (NON-EQ) of the cluster and
  2) a letter ``T`` or ``F`` in the second position of the same line.
  The number of NON-EQ is the total number of magnetic centers of the cluster
  which cannot be related by point group symmetry.
  In the second position the answer to the question: *Have all NON-EQ centers been computed ab initio?*
  is given: ``T`` for *True* and ``F`` for *False*.
  On the third position, the answer to the question: *Are the rassi.h5 files to be read for input?* 
  is given. For the current status, the letter ``F`` is the only option.
  On the following line the program will read NON-EQ values specifying the
  number of equivalent centers of each type.
  On the following line the program will read NON-EQ integer numbers specifying
  the number of low-lying spin-orbit functions from each center forming the local
  exchange basis.

  Some examples valid for situations where all sites have been
  computed *ab initio* (case ``T``, *True*):

  .. class:: poly_aniso

  +----------------------------------------+----------------------------------------+----------------------------------------+
  | ::                                     | ::                                     | ::                                     |
  |                                        |                                        |                                        |
  |   NNEQ                                 |   NNEQ                                 |   NNEQ                                 |
  |   2  T  F                              |   3  T  F                              |   6  T  F                              |
  |   1  2                                 |   2  1  1                              |   1  1  1  1  1  1                     |
  |   2  2                                 |   4  2  3                              |   2  4  3  5  2  2                     |
  +----------------------------------------+----------------------------------------+----------------------------------------+
  | There are two kinds of magnetic centers| There are three kinds of magnetic      | There are 6 kinds of magnetic centers  |
  | in the cluster; both have been computed| centers in the cluster; all three have | in the cluster; all six have been      |
  | ab initio; the cluster consists of 3   | been computed ab initio; the cluster   | computed ab initio; the cluster        |
  | magnetic centers: one center of the    | consists of four magnetic centers: two | consists of 6 magnetic centers: one    |
  | first kind and two centers of the      | centers of the first kind, one center  | center of each kind. From the center of|
  | second kind. From each center we take  | of the second kind and one center of   | the first kind we take into exchange   |
  | into the exchange coupling only the    | the third kind. From each of the       | coupling two spin-orbit states, four   |
  | ground doublet. As a result the        | centers of the first kind we take into | states from the second center, three   |
  | :math:`N_{\text{exch}}=2^1 \times      | exchange coupling four spin-orbit      | states from the third center, five     |
  | 2^2=8` :file:`aniso_1.input` (for ---  | states, two states from the second kind| states from the fourth center and two  |
  | type 1) and :file:`aniso_2.input`      | and three states from the third center.| states from the fifth and sixth        |
  | (for --- type 2) files must be present.| As a result the                        | centers. As a result the               |
  |                                        | :math:`N_{\text{exch}}=4^2 \times 2^1  | :math:`N_{\text{exch}}=2^1 \times 4^1  |
  |                                        | \times 3^1=96`. Three files            | \times 3^1 \times 5^1 \times 2^1 \times|
  |                                        | :file:`aniso_i.input` for each center  | 2^1=480`. Six files                    |
  |                                        | (:math:`i=1,2,3`) must be present.     | :file:`aniso_i.input` for each center  |
  |                                        |                                        | (:math:`i=1,2,\ldots,6`) must be       |
  |                                        |                                        | present.                               |
  +----------------------------------------+----------------------------------------+----------------------------------------+

  Only in cases when some centers have NOT been computed ab initio (i.e. for which no :file:`aniso_i.input` file exists),
  the program will read an additional line consisting of NON-EQ letters (``A`` or ``B``) specifying the type of each of
  the NON-EQ centers:
  ``A`` --- the center is computed ab initio and ``B`` --- the center is considered isotropic.
  On the following ``number-of-B-centers`` line(s) the isotropic :math:`g` factors of the
  center(s) defined as ``B`` are read. The spin of the ``B`` center(s) is defined: :math:`S=(N-1)/2`,
  where :math:`N` is the corresponding number of states to be taken into the exchange coupling
  for this particular center.

  Some examples valid for mixed situations: the system consists of centers computed *ab initio* and
  isotropic centers (case ``F``, *False*):

.. class:: poly_aniso

  +----------------------------------------+----------------------------------------+----------------------------------------+
  | ::                                     | ::                                     | ::                                     |
  |                                        |                                        |                                        |
  |   NNEQ                                 |   NNEQ                                 |   NNEQ                                 |
  |   2  F  F                              |   3  F  F                              |   6  T  F                              |
  |   1  2                                 |   2  1  1                              |   1  1  1  1  1  1                     |
  |   2  2                                 |   4  2  3                              |   2  4  3  5  2  2                     |
  |   A  B                                 |   A  B  B                              |   B  B  A  A  B  A                     |
  |   2.3                                  |   2.1                                  |   2.12                                 |
  |                                        |   2.0                                  |   2.43                                 |
  |                                        |                                        |   2.00                                 |
  +----------------------------------------+----------------------------------------+----------------------------------------+
  | There are two kinds of magnetic centers| There are three kinds of magnetic      | There are six kinds of magnetic centers|
  | in the cluster; the center of the first| centers in the cluster; the first      | in the cluster; only three centers have|
  | type has been computed *ab initio*,    | center type has been computed *ab*     | *been computed *ab initio*, while the  |
  | while the centers of the second type   | initio*, while the centers of the      | other three centers are considered     |
  | are considered isotropic with          | second and third types are considered  | isotropic; the :math:`g` factor of the |
  | :math:`g=2.3`; the cluster consists of | isotropic with :math:`g=2.1` (second   | first center is 2.12 (:math:`S=1/2`);  |
  | three magnetic centers: one center of  | type) and :math:`g=2.0` (third type);  | of the second center 2.43              |
  | the first kind and two centers of the  | the cluster consists of four magnetic  | (:math:`S=3/2`); of the fifth center   |
  | second kind. Only the ground doublet   | centers: two centers of the first kind,| 2.00 (:math:`S=1/2`); the entire       |
  | state from each center is considered   | one center of the second kind and one  | cluster consists of six magnetic       |
  | for the exchange coupling. As a result | center of the third kind. From each of | centers: one center of each kind. From |
  | the :math:`N_{\text{exch}}=2^1 \times  | the centers of the first kind, four    | the center of the first kind, two      |
  | 2^2=8`. File :file:`aniso_1.input` (for| spin-orbit states are considered for   | spin-orbit states are considered in the|
  | --- type 1) must be present.           | the exchange coupling, two states from | exchange coupling, four states from the|
  |                                        | the second kind and three states from  | second center, three states from the   |
  |                                        | the center of the third kind. As a     | third center, five states from the     |
  |                                        | result the :math:`N_{\text{exch}}=4^2  | fourth center and two states from the  |
  |                                        | \times 2^1 \times 3^1=96`. The file    | fifth and sixth centers. As a result   |
  |                                        | :file:`aniso_1.input` must be present. | the :math:`N_{\text{exch}}=2^1 \times  |
  |                                        |                                        | 4^1 \times 3^1 \times 5^1 \times 2^1   |
  |                                        |                                        | \times 2^1=480`. Three files           |
  |                                        |                                        | :file:`aniso_3.input` and              |
  |                                        |                                        | :file:`aniso_4.input` and              |
  |                                        |                                        | :file:`aniso_6.input` must be present. |
  +----------------------------------------+----------------------------------------+----------------------------------------+

  There is no maximal value for :kword:`NNEQ`, although the calculation becomes quite heavy in case the number of
  exchange functions is large.

.. xmldoc:: <KEYWORD MODULE="POLY_ANISO" NAME="NNEQ" APPEAR="Definition of input magnetic sites" KIND="CUSTOM" LEVEL="BASIC">
              %%Keyword: NNEQ <basic>
              <HELP>
                 This keyword defines several important parameters of the calculation. On the
                 first line after the keyword the program reads 2 values:
                 1) the number of types of different magnetic centers (NON-EQ) of the cluster and
                 2) a letter ``T`` or ``F`` in the second position of the same line.
                 The number of NON-EQ is the total number of magnetic centers of the cluster
                 which cannot be related by point group symmetry.
                 In the second position the answer to the question: *Have all NON-EQ centers been computed ab initio?*
                 is given: ``T`` for *True* and ``F`` for *False*.
                 On the third position, the answer to the question: *Are the rassi.h5 files to be read for input?* 
                 is given. For the current status, the letter ``F`` is the only option.
                 On the following line the program will read NON-EQ values specifying the
                 number of equivalent centers of each type.
                 On the following line the program will read NON-EQ integer numbers specifying
                 the number of low-lying spin-orbit functions from each center forming the local
                 exchange basis.
              </HELP>
              </KEYWORD>


:kword:`SYMM`
  Specifies rotation matrices to symmetry equivalent sites. This keyword is mandatory in the case more centers of a given type are present in the calculation.
  This keyword is mandatory when the calculated polynuclear compound has exact crystallographic point group symmetry. In other words, when the number of
  equivalent centers of any kind :math:`i` is larger than 1, this keyword must be employed. Here the rotation matrices from the one
  center to all the other of the same type are declared.
  On the following line the program will read the number ``1`` followed on the next lines by as many :math:`3\times3` rotation matrices as the total number of
  equivalent centers of type ``1``. Then the rotation matrices of centers of type ``2``, ``3`` and so on, follow in the same format.
  When the rotation matrices contain irrational numbers (e.g. :math:`\sin{\frac{\pi}{6}}=\frac{\sqrt{3}}{2}`), then more digits than presented in the examples
  below are advised to be given: :math:`\frac{\sqrt{3}}{2}=0.86602540378`.

  Examples:

  .. class:: poly_aniso

  +----------------------------------------+----------------------------------------+----------------------------------------+
  | ::                                     | ::                                     | ::                                     |
  |                                        |                                        |                                        |
  |     NNEQ                               |   NNEQ                                 |   NNEQ                                 |
  |     2  F  F                            |   3  F  F                              |   6  F  F                              |
  |     1  2                               |   2  1  1                              |   1  1  1  1  1  1                     |
  |     2  2                               |   4  2  3                              |   2  4  3  5  2  2                     |
  |     A  B                               |   A  B  B                              |   B  B  A  A  B  A                     |
  |     2.3                                |   2.1                                  |   2.12                                 |
  |                                        |   2.0                                  |   2.43                                 |
  |     SYMM                               |   2.0                                  |   2.00                                 |
  |     1                                  |                                        |                                        |
  |     1.0 0.0 0.0                        |   SYMM                                 |                                        |
  |     0.0 1.0 0.0                        |   1                                    |                                        |
  |     0.0 0.0 1.0                        |   1.0 0.0 0.0                          |                                        |
  |     2                                  |   0.0 1.0 0.0                          |                                        |
  |     1.0 0.0 0.0                        |   0.0 0.0 1.0                          |                                        |
  |     0.0 1.0 0.0                        |   0.0 -1.0 0.0                         |                                        |
  |     0.0 0.0 1.0                        |   1.0 0.0  0.0                         |                                        |
  |     -1.0 0.0 0.0                       |   0.0 0.0  1.0                         |                                        |
  |     0.0 -1.0 0.0                       |   2                                    |                                        |
  |     0.0 0.0 -1.0                       |   1.0 0.0 0.0                          |                                        |
  |                                        |   0.0 1.0 0.0                          |                                        |
  |                                        |   0.0 0.0 1.0                          |                                        |
  |                                        |   3                                    |                                        |
  |                                        |   1.0 0.0 0.0                          |                                        |
  |                                        |   0.0 1.0 0.0                          |                                        |
  |                                        |   0.0 0.0 1.0                          |                                        |
  +----------------------------------------+----------------------------------------+----------------------------------------+
  | The cluster computed here is a         | In this input a tetranuclear compound  | In this case the computed system has no|
  | trinuclear compound, with one center   | is defined, all centers are computed ab| symmetry. Therefore, the :kword:`SYMM` |
  | computed ab initio, while the other two| initio. There are two centers of type  | keyword may be skipped.                |
  | centers, related to each other by      | "1", related one to each other by      |                                        |
  | inversion, are considered isotropic    | :math:`C_2` symmetry around the        |                                        |
  | with :math:`g_x=g_y=g_z=2.3`. The      | Cartesian Z axis. Therefore the        |                                        |
  | rotation matrix for the first center is| :kword:`SYMM` keyword is mandatory.    |                                        |
  | :math:`I` (identity, unity) since the  | There are two matrices for centers of  |                                        |
  | center is unique. For the centers of   | type 1, and one matrix (identity) for  |                                        |
  | type 2, there are two matrices         | the centers of type 2 and type 3.      |                                        |
  | :math:`3\times3` since we have two     |                                        |                                        |
  | centers in the cluster. The rotation   |                                        |                                        |
  | matrix of the first center of type 2 is|                                        |                                        |
  | Identity while the rotation matrix for |                                        |                                        |
  | the equivalent center of type 2 is the |                                        |                                        |
  | inversion matrix.                      |                                        |                                        |
  +----------------------------------------+----------------------------------------+----------------------------------------+

  More examples are given in the *Tutorial* section.


.. xmldoc:: <KEYWORD MODULE="POLY_ANISO" NAME="SYMM" APPEAR="Definition of symmetry of the polynuclear cluster, if any" KIND="CUSTOM" LEVEL="BASIC">
              %%Keyword: SYMM <basic>
              <HELP>
                 Specifies rotation matrices to symmetry equivalent sites. This keyword is mandatory in the case more centers of a given type are present in the calculation.
                 This keyword is mandatory when the calculated polynuclear compound has exact crystallographic point group symmetry. In other words, when the number of
                 equivalent centers of any kind :math:`i` is larger than 1, this keyword must be employed. Here the rotation matrices from the one
                 center to all the other of the same type are declared.
                 On the following line the program will read the number ``1`` followed on the next lines by as many :math:`3\times3` rotation matrices as the total number of
                 equivalent centers of type ``1``. Then the rotation matrices of centers of type ``2``, ``3`` and so on, follow in the same format.
                 When the rotation matrices contain irrational numbers (e.g. :math:`\sin{\frac{\pi}{6}}=\frac{\sqrt{3}}{2}`), then more digits than presented in the examples
                 below are advised to be given: :math:`\frac{\sqrt{3}}{2}=0.86602540378`.
              </HELP>
              </KEYWORD>









*Keywords defining the magnetic exchange interactions*

This section defines the keywords used to set up the interacting pairs of magnetic centers
and the corresponding exchange interactions.

A few words about the numbering of the magnetic centers of the
cluster in the :program:`POLY_ANISO`. First all equivalent centers of the type 1 are
numbered, then all equivalent centers of the type 2, etc. These labels of the magnetic
centers are used further for the declaration of the magnetic coupling.
The pseudo-code is: ::

  k=0
  Do i=1, number-of-non-equivalent-sites
    Do j=1, number-of-equivalent-sites-of-type(i)
       k=k+1
       site-number(i,j)=k
    End Do
  End Do

.. class:: keywordlist

:kword:`PAIR` or :kword:`LIN1`
  Specifies the Lines interaction(s) between metal pairs. One parameter per interacting pair is required.

  ::

    LIN1
       READ number-of-interacting-pairs
       Do i=1, number-of-interacting-pairs
          READ site-1, site-2,   J
       End Do

.. xmldoc:: <KEYWORD MODULE="POLY_ANISO" NAME="PAIR" KIND="REALS_COMPUTED" SIZE="3" LEVEL="UNDOCUMENTED" />
.. xmldoc:: <KEYWORD MODULE="POLY_ANISO" NAME="LIN1" KIND="REALS_COMPUTED" SIZE="3" LEVEL="UNDOCUMENTED" />

:kword:`ALIN` :kword:`LIN3`
  Specifies the anisotropic interactions between metal pairs. Three parameters per interacting pair are required.

  ::

    LIN3
       READ number-of-interacting-pairs
       Do i=1, number-of-interacting-pairs
          READ site-1, site-2,   Jxx, Jyy, Jzz
       End Do

  :math:`J_{\alpha\beta}`, where :math:`\alpha` and :math:`\beta` are main values of the Cartesian components of the (:math:`3\times3`) matrix defining the exchange interaction between site-1 and site-2.

.. xmldoc:: <KEYWORD MODULE="POLY_ANISO" NAME="LIN3" KIND="REALS_COMPUTED" SIZE="5" LEVEL="UNDOCUMENTED" />




:kword:`LIN9`
  Specifies the full anisotropic interaction matrices between metal pairs. Nine parameters per interacting pair is required.

  ::

    LIN9
       READ number-of-interacting-pairs
       Do i=1, number-of-interacting-pairs
          READ site-1, site-2,   Jxx, Jxy, Jxz,   Jyx, Jyy, Jyz,  Jzx, Jzy, Jzz
       End Do

 :math:`J_{\alpha\beta}`, where :math:`\alpha` and :math:`\beta` are main values of the Cartesian components of the (:math:`3\times3`) matrix defining the exchange interaction between site-1 and site-2.

.. xmldoc:: <KEYWORD MODULE="POLY_ANISO" NAME="LIN9" KIND="REALS_COMPUTED" SIZE="11" LEVEL="UNDOCUMENTED" />



:kword:`COOR`
  Specifies the symmetrized coordinates of the metal sites. This keyword enables computation of dipole-dipole
  magnetic interaction between metal sites defined in the keywords :kword:`PAIR`, :kword:`ALIN`, :kword:`LIN1`, :kword:`LIN3` or :kword:`LIN9`.

  ::

    COOR
       Do i=1, number-of-non-equivalent-sites
          READ coordinates of center 1
          READ coordinates of center 2
          ...
       End Do

.. xmldoc:: <KEYWORD MODULE="POLY_ANISO" NAME="COOR" KIND="REALS_LOOKUP" SIZE="3NNEQ" LEVEL="BASIC">
              %%Keyword: COOR <basic>
              <HELP>
                 Specifies the symmetrized coordinates of the metal sites. This keyword enables computation of dipole-dipole interaction.
              </HELP>
              </KEYWORD>






*Other keywords*

Normally :program:`POLY_ANISO` runs without specifying any of the following keywords.

Argument(s) to a keyword are always supplied on the next line of the input file.

Optional general keywords to control the input
..............................................

.. class:: keywordlist

:kword:`MLTP`
  The number of molecular multiplets (i.e. groups of spin-orbital eigenstates) for which
  :math:`g`, :math:`D` and higher magnetic tensors will be calculated (default :kword:`MLTP`\=1).
  The program reads two lines: the first is the number of multiplets (:math:`N_{\text{MULT}}`) and the
  second the array of :math:`N_{\text{MULT}}` numbers specifying the dimension (multiplicity) of each multiplet.

  Example: ::

      MLTP
      10
      2 4 4 2 2   2 2 2 2 2

    :program:`POLY_ANISO` will compute the :math:`g` and :math:`D{-}` tensors for 10 groups of states.
    The groups 1 and 4--10 are doublets (:math:`\tilde{S}=\ket{1/2}`), while the groups 2 and 3 are quadruplets,
    having the effective spin :math:`\tilde{S}=\ket{3/2}`. For the latter cases, the ZFS (:math:`D{-}`) tensors will be computed.

  .. xmldoc:: <KEYWORD MODULE="POLY_ANISO" NAME="MLTP" KIND="INTS_COMPUTED" SIZE="1" LEVEL="BASIC" DEFAULT_VALUE="1">
              %%Keyword: MLTP <basic>
              <HELP>
              The number of molecular multiplets (i.e. groups of spin-orbital eigenstates) for
              which g, D and higher magnetic tensors will be calculated.
              The program reads two lines: the first is the number of multiplets (NMULT) and
              on the second line the array of NMULT numbers specifying the dimension of each multiplet.
              By default, the code will first analyze the energy spectra by itself and will
              compute the g and D tensors for ten low-lying groups of states. By using this
              keyword the user overwrites the default.
              </HELP>
              </KEYWORD>

:kword:`TINT`
  Specifies the temperature points for the evaluation of the magnetic susceptibility. The program will read four numbers: :math:`T_{\text{min}}`, :math:`T_{\text{max}}`, :math:`n_T`.

  .. container:: list

    :math:`T_{\text{min}}` --- the minimal temperature (Default 0.0 K)

    :math:`T_{\text{max}}` --- the maximal temperature (Default 300.0 K)

    :math:`n_T` --- number of temperature points (Default 101)

  Example: ::

    TINT
    0.0  330.0  331

  :program:`POLY_ANISO` will compute temperature dependence of the magnetic susceptibility in 331 points evenly distributed in temperature interval: 0.0 K -- 330.0 K.

  .. xmldoc:: <KEYWORD MODULE="POLY_ANISO" NAME="TINT" KIND="REAL" LEVEL="BASIC">
              %%Keyword: TINT <basic>
              <HELP>
              Specifies the temperature points for the evaluation of the magnetic susceptibility.
              The program will read four numbers: Tmin, Tmax, nT, and dltT0. Units of temperature = Kelvin (K).
              ||Tmin  -- the minimal temperature (Default 0.0 K)
              ||Tmax  -- the maximal temperature (Default 300.0 K)
              ||nT    -- number of temperature points (Default 101)
              </HELP>
              </KEYWORD>

:kword:`HINT`
  Specifies the field points for the evaluation of the magnetization in a certain direction. The program will read four numbers: :math:`H_{\text{min}}`, :math:`H_{\text{max}}` and :math:`n_H`

  .. container:: list

    :math:`H_{\text{min}}` --- the minimal field (Default 0.0 T)

    :math:`H_{\text{max}}` --- the maximal filed (Default 10.0 T)

    :math:`n_H` --- number of field points (Default 101)

  .. compound::

    Example: ::

      HINT
      0.0  20.0  201

    :program:`POLY_ANISO` will compute the molar magnetization in 201 points evenly distributed in field interval: 0.0 T -- 20.0 T.

  .. xmldoc:: <KEYWORD MODULE="POLY_ANISO" NAME="HINT" KIND="REAL" LEVEL="BASIC">
              %%Keyword: HINT <basic>
              <HELP>
              Specifies the field points for the evaluation of the molar magnetization.
              The program will read four numbers: Hmin, Hmax, nH, and dltH0. Units of magnetic field = tesla (T).
              ||
              ||Hmin -- the minimal field (Default 0.0 T)
              ||Hmax -- the maximal field (Default 300.0 T)
              ||nH   -- number of field points (Default 101)
              </HELP>
              </KEYWORD>

:kword:`TMAG`
  Specifies the temperature(s) at which the field-dependent magnetization is calculated. Default is one temperature point, :math:`T`\=2.0 K.
  Example: ::

    TMAG
    6   1.8 2.0 2.4  2.8 3.2 4.5

  .. xmldoc:: <KEYWORD MODULE="POLY_ANISO" NAME="TMAG" KIND="REAL" LEVEL="BASIC">
              %%Keyword: TMAG <basic>
              <HELP>
              Specifies the temperature at which the field-dependent magnetization is calculated. Default is 2.0 K
              </HELP>
              </KEYWORD>

:kword:`ENCU`
  .. compound::

    This flag is used to define the cut-off energy for the lowest states for which Zeeman interaction is taken into account exactly. The contribution to the magnetization coming from states that are higher in energy than :math:`E` (see below) is done by second order perturbation theory. The program will read two integer numbers: :math:`N_K` and :math:`M_G`. Default values are: :math:`N_K`\=100, :math:`M_G`\=100.

    .. math:: E=N_K \cdot k_{\text{B}} \cdot T_{\text{MAG}} + M_G \cdot \mu_{\text{Bohr}} \cdot H_{\text{max}}

    The field-dependent magnetization is calculated at the (highest) temperature value defined in either :kword:`TMAG` or :kword:`HEXP`.
    Example: ::

      ENCU
      250  150

    If :math:`H_{\text{max}}` = 10 T and :kword:`TMAG` = 1.8 K, then the cut-off energy is:

    .. math:: E=100 \cdot 250 \cdot k_{\text{B}} \cdot 1.8 + 150 \cdot \mu_{\text{Bohr}} \cdot 10 = 1013.06258\,\text{cm}^{-1}

    This means that the magnetization coming from all spin-orbit states with energy lower than :math:`E=1013.06258\,\text{cm}^{-1}` will be computed exactly.
    :kword:`ERAT`, :kword:`NCUT` and :kword:`ENCU` are mutually exclusive.

  .. xmldoc:: <KEYWORD MODULE="POLY_ANISO" NAME="ENCU" KIND="INT" LEVEL="BASIC">
              %%Keyword: ENCU <basic>
              <HELP>
              This keyword is used to define the cut-off energy for the lowest states for which
              Zeeman interaction is taken into account exactly. The contribution to the
              magnetization coming from states that are higher in energy than E (see below)
              is done by second order perturbation theory. The program will read two integer
              numbers: NK and MG. Default values are: NK=100, MG=100. The field-dependent magnetization
              is calculated at the temperature value TMAG.
              </HELP>
              </KEYWORD>

:kword:`ERAT`
  .. compound::

    This flag is used to define the cut-off energy for the lowest states for which Zeeman interaction
    is taken into account exactly. The contribution to the molar magnetization coming from states that
    are higher in energy than :math:`E` (see below) is done by second order perturbation theory.
    The program reads one real number in the domain (0.0--1.0). Default is 1.0 (all exchange states are
    included in the Zeeman interaction).

    .. math:: E = \text{ERAT} \cdot \text{Maximal-spread-of-exchange-splitting}

    The field-dependent magnetization is calculated at all temperature points defined in either :kword:`TMAG` or :kword:`HEXT`.
    Example: ::

      ERAT
      0.75

    :kword:`ERAT`, :kword:`NCUT` and :kword:`ENCU` are mutually exclusive.

  .. xmldoc:: <KEYWORD MODULE="POLY_ANISO" NAME="ERAT" KIND="INT" LEVEL="BASIC">
              %%Keyword: ERAT <basic>
              <HELP>
              This keyword is used to define the cut-off energy for the lowest states for which
              Zeeman interaction is taken into account exactly. The contribution to the
              magnetization coming from states that are higher in energy than E (see below)
              is done by second order perturbation theory. The program will read one real number in the domain 0.0-1.0.
              The field-dependent magnetization
              is calculated at the temperature value TMAG.
              </HELP>
              </KEYWORD>

:kword:`NCUT`
  .. compound::

    This flag is used to define the number of low-lying exchange states for which Zeeman interaction is taken into
    account exactly. The contribution to the magnetization coming from the remaining exchange states is done by second
    order perturbation theory. The program will read one integer number. The field-dependent magnetization is calculated at all temperature points defined in either :kword:`TMAG` or :kword:`HEXT`.
    Example: ::

      NCUT
      125

    In case the defined number is larger than the total number of exchange states in the calculation (:math:`N_{\text{exch}}`), then :math:`n_{\text{Cut}}` is set to be equal to :math:`N_{\text{exch}}`.
    :kword:`ERAT`, :kword:`NCUT` and :kword:`ENCU` are mutually exclusive.

:kword:`MVEC`

    Defines the number of directions for which the magnetization vector will be computed.
    On the first line below the keyword, the number of directions should be mentioned (:math:`N_{\text{DIR}}`. Default 0).
    The program will read :math:`N_{\text{DIR}}` lines for cartesian coordinates specifying the direction :math:`i` of the
    applied magnetic field (:math:`\theta_i` and :math:`\phi_i`). These values may be arbitrary real numbers.
    The direction(s) of applied magnetic field are obtained by normalizing the length of each vector to one.
    Example: ::

      MVEC
      4
      0.0000  0.0000   0.1000
      1.5707  0.0000   2.5000
      1.5707  1.5707   1.0000
      0.4257  0.4187   0.0000

    The above input requests computation of the magnetization vector in four directions of applied field.
    The actual directions on the unit sphere are: ::

      4
      0.00000  0.00000  1.00000
      0.53199  0.00000  0.84675
      0.53199  0.53199  0.33870
      0.17475  0.17188  0.00000

  .. xmldoc:: <KEYWORD MODULE="POLY_ANISO" NAME="MVEC" KIND="REALS_COMPUTED" SIZE="3" LEVEL="BASIC">
              %%Keyword: MVEC <basic>
              <HELP>
              Defines the number of directions for which the magnetization vector will be computed.
              On the first line below the keyword, the number of directions should be mentioned (NDIR. Default 0).
              The program will read NDIR lines for spherical coordinates specifying the direction
              "i" of the magnetic field (theta_i and phi_i). These values should be in radians.
              </HELP>
              </KEYWORD>


:kword:`ZEEM`

    Defines the number of directions for which the Zeeman energy will be computed /saved / plotted.
    On the first line below the keyword, the number of directions should be mentioned (:math:`N_{\text{DIR}}`. Default 0).
    The program will read :math:`N_{\text{DIR}}` lines for cartesian coordinates specifying the direction :math:`i` of the
    applied magnetic field (:math:`\theta_i` and :math:`\phi_i`). These values may be arbitrary real numbers.
    The direction(s) of applied magnetic field are obtained by normalizing the length of each vector to one.
    Example: ::

      MVEC
      4
      0.0000  0.0000   0.1000
      1.5707  0.0000   2.5000
      1.5707  1.5707   1.0000
      0.4257  0.4187   0.0000

    The above input requests computation of the magnetization vector in four directions of applied field.
    The actual directions on the unit sphere are: ::

      4
      0.00000  0.00000  1.00000
      0.53199  0.00000  0.84675
      0.53199  0.53199  0.33870
      0.17475  0.17188  0.00000

  .. xmldoc:: <KEYWORD MODULE="POLY_ANISO" NAME="ZEEM" KIND="REALS_COMPUTED" SIZE="3" LEVEL="BASIC">
              %%Keyword: ZEEM <basic>
              <HELP>
              Defines the number of directions for which the magnetization vector will be computed.
              On the first line below the keyword, the number of directions should be mentioned (NDIR. Default 0).
              The program will read NDIR lines for spherical coordinates specifying the direction
              "i" of the magnetic field (theta_i and phi_i). These values should be in radians.
              </HELP>
              </KEYWORD>




:kword:`MAVE`
  This keyword specifies the grid density used for the computation of powder molar
  magnetization. The program uses Lebedev--Laikov distribution of points on the unit sphere.
  The program reads two integer numbers: :math:`n_{\text{sym}}` and :math:`n_{\text{grid}}`. The :math:`n_{\text{sym}}` defines which
  part of the sphere is used for averaging. It takes one of the three values: 1 (half-sphere),
  2 (a quater of a sphere) or 3 (an octant of the sphere). :math:`n_{\text{grid}}` takes values from 1
  (the smallest grid) till 32 (the largest grid, i.e. the densiest). The default is to
  consider integration over a half-sphere (since :math:`M(H)=-M(-H)`): :math:`n_{\text{sym}}=1` and :math:`n_{\text{grid}}=15`
  (i.e. 185 points distributed over half-sphere). In case of symmetric compounds, powder
  magnetization may be averaged over a smaller part of the sphere, reducing thus the number
  of points for the integration. The user is responsible to choose the appropriate integration scheme.
  Default value for :math:`n_{\text{grid}}=15` (185 directions equally distributed in the given area).
  Note that the program's default is rather conservative.

  .. xmldoc:: <KEYWORD MODULE="POLY_ANISO" NAME="MAVE" KIND="INT" LEVEL="BASIC">
              %%Keyword: MAVE <basic>
              <HELP>
              Specifies the number of directions of the applied magnetic field for the computation
              of the powder molar magnetization. The program will read two numbers: N_theta and N_phi.
              ||N_theta -- number of "theta" points in the interval (0, pi/2) (i.e. on the Z axis ) (Default 12)
              ||N_phi   -- number of  "phi"  points in the interval (0, 2*pi).(i.e. on the equator) (Default 24)
              The number of directions over which the actual averaging will take place is roughly the product of N_theta and N_phi.
              </HELP>
              </KEYWORD>

:kword:`TEXP`
  This keyword allows computation of the magnetic susceptibility :math:`\chi T(T)` at experimental points.
  On the line below the keyword, the number of experimental points :math:`N_T` is defined, and on the next :math:`N_T` lines the program reads the experimental temperature (in K) and the experimental magnetic susceptibility (in :math:`\text{cm}^3\text{K}\text{mol}^{-1}`). :kword:`TEXP` and :kword:`TINT` keywords are mutually exclusive. The magnetic susceptibility routine will also print the standard deviation from the experiment.

  ::

    TEXP
       READ  number-of-T-points
       Do i=1, number-of-T-points
          READ ( susceptibility(i, Temp), TEMP = 1, number-of-T-points )
       End Do

  .. xmldoc:: <KEYWORD MODULE="POLY_ANISO" NAME="TEXP" KIND="REALS_COMPUTED" SIZE="2" LEVEL="BASIC">
              %%Keyword: TEXP <basic>
              <HELP>
              This keyword allows computation of the magnetic susceptibility at experimental
              temperature points. On the line below the keyword, the number of experimental
              points NT is defined, and on the next NT lines the program reads the experimental
              temperature (in K) and the experimental magnetic susceptibility (in cm^3Kmol^{-1}).
              TEXP and TINT keywords are mutually exclusive. The POLY_ANISO will also print the
              standard deviation from the experiment.
              </HELP>
              </KEYWORD>

:kword:`HEXP`
  This keyword allows computation of the molar magnetization :math:`M_{\text{mol}} (H)` at experimental points.
  On the line below the keyword,the number of experimental points :math:`N_H` is defined, and on the next :math:`N_H` lines the program reads the experimental field intensity (tesla) and the experimental magnetization (in :math:`\mu_{\text{Bohr}}`). :kword:`HEXP` and :kword:`HINT` are mutually exclusive. The magnetization routine will print the standard deviation from the experiment.

  ::

    HEXP
       READ  number-of-T-points-for-M,  all-T-points-for-M-in-K
       READ  number-of-field-points
       Do i=1, number-of-field-points
          READ ( Magn(i, iT), iT=1, number-of-T-points-for-M )
       End Do

  .. xmldoc:: <KEYWORD MODULE="POLY_ANISO" NAME="HEXP" KIND="CUSTOM" LEVEL="BASIC">
              %%Keyword: HEXP <basic>
              <HELP>
              This keyword allows computation of the molar magnetization at experimental field points.
              On the line below the keyword,the number of experimental points NH is defined, and on
              the next NH lines the program reads the experimental field strength (tesla) and the
              experimental magnetization (in Bohr magnetons). HEXP and HINT are mutually exclusive.
              The POLY_ANISO will print the standard deviation from the experiment.
              </HELP>
              </KEYWORD>

:kword:`ZJPR`
  This keyword specifies the value (in :math:`\text{cm}^{-1}`) of a phenomenological parameter of a mean molecular field acting on the spin of the complex (the average intermolecular exchange constant). It is used in the calculation of all magnetic properties (not for spin Hamiltonians) (Default is 0.0)

  .. xmldoc:: <KEYWORD MODULE="POLY_ANISO" NAME="ZJPR" KIND="REAL" LEVEL="BASIC">
              %%Keyword: ZJPR <basic>
              <HELP>
              This keyword specifies the value (in cm^-1) of a phenomenological parameter of a
              mean molecular field acting on the spin of the complex (the average intermolecular
              exchange constant). It is used in the calculation of all magnetic properties (not for
              spin Hamiltonians) (Default is 0.0)
              </HELP>
              </KEYWORD>


:kword:`ABCC`
  This keyword will enable computation of magnetic and anisotropy axes in the
  crystallographic :math:`abc` system. On the next line, the program will read six real
  values, namely :math:`a`, :math:`b`, :math:`c`, :math:`\alpha`, :math:`\beta`, and :math:`\gamma`, defining the
  crystal lattice. On the second line, the program will read the Cartesian coordinates
  of the magnetic center. The computed values in the output correspond to the
  crystallographic position of three "dummy atoms" located on the corresponding anisotropy axes, at the distance of 1 ångstrom from the metal site. ::

    ABCC
    20.17   19.83   18.76    90  120.32  90
    12.329  13.872  1.234

  .. xmldoc:: <KEYWORD MODULE="POLY_ANISO" NAME="ABCC" KIND="STRING" LEVEL="BASIC">
              %%Keyword: ABCC <basic>
              <HELP>
              This keyword will enable computation of magnetic and anisotropy axes in the
              crystallographic abc system. On the next line, the program will read six real
              values, namely (a, b, c, alpha, beta, and gamma), defining the crystal lattice.
              On the second line, the program will read the Cartesian coordinates of the
              magnetic center. The computed values in the output correspond to the crystallographic
              position of three "dummy atoms" located on the corresponding anisotropy axes, at the
              distance of 1.0 angstrom from the metal site.
              </HELP>
              </KEYWORD>


:kword:`XFIE`
  This keyword specifies the value (in :math:`\text{T}`) of applied magnetic field
  for the computation of magnetic susceptibility by :math:`dM/dH` and :math:`M/H` formulas.
  A comparison with the usual formula (in the limit of zero applied field) is provided.
  (Default is 0.0)

  .. xmldoc:: <KEYWORD MODULE="POLY_ANISO" NAME="XFIE" KIND="REAL" LEVEL="BASIC">
              %%Keyword: XFIE <basic>
              <HELP>
              This keyword specifies the value (in Tesla) of applied magnetic field
              for the computation of magnetic susceptibility by: dM/dH and M/H formulas.
              A comparison with the usual formula (in the limit of zero applied field) is provided.
              (Default is 0.0)
              </HELP>
              </KEYWORD>

:kword:`PRLV`
  This keyword controls the print level.

  .. container:: list

    2 --- normal. (Default)

    3 or larger (debug)

  .. xmldoc:: <KEYWORD MODULE="POLY_ANISO" NAME="PRLV" KIND="INT" LEVEL="BASIC">
              %%Keyword: PRLV <basic>
              <HELP>
              This keyword controls the print level.
              ||2 -- normal. (Default)
              ||3 or larger (debug)
              </HELP>
              </KEYWORD>


:kword:`PLOT`
  This keyword will generate a few plots (png or eps format) via an interface to the linux program *gnuplot*. 
  The interface generates a datafile, a gnuplot script and attempts execution of the script for generation of the image. 
  The plots are generated only if the respective function is invoked. The magnetic susceptibility, molar magnetisation and blocking barrier (UBAR) plots are generated.
  The files are named: file:`XT.dat`, file:`XT.plt`, file:`XT.png`, file:`MH.dat`, file:`MH.plt`, file:`MH.png`, file:`BARRIER_TME.dat`, file:`BARRIER_ENE.dat`, file:`BARRIER.plt` and file:`BARRIER.png`.


  .. xmldoc:: <KEYWORD MODULE="SINGLE_ANISO" NAME="PLOT" KIND="SINGLE" LEVEL="BASIC">
              %%Keyword: PLOT <basic>
              <HELP>
              This keyword will generate a few plots (png or eps format) via an interface to the linux program "gnuplot".
              The interface generates a datafile, a gnuplot script and attempts execution of the script for generation of the image.
              The plots are generated only if the respective function is invoked. The magnetic susceptibility, molar magnetisation and blocking barrier (UBAR) plots are generated.
              The files are named: `XT.dat`, `XT.plt`, `XT.png`, `MH.dat`, `MH.plt`, `MH.png`, `BARRIER_TME.dat`, `BARRIER_ENE.dat`, `BARRIER.plt` and `BARRIER.png`.
              </HELP>
              </KEYWORD>



.. xmldoc:: </MODULE>
