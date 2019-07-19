.. _UG\:sec\:the_basis_set_libraries:

The Basis Set Libraries
=======================

.. only:: html

  .. contents::
     :local:
     :backlinks: none

The basis sets library contains both all-electron and effective core
potentials. They will be briefly described below and we
refer to the publications for more details. The user can also add new
basis sets to the basis directory and the structure of the file will therefore
be described below.

Dummy atoms
-----------

Note that to use dummy atoms the user should employ the basis set label "``X....``". This
will signify centers associated with no charge and no basis functions.

The All Electron Basis Set Library
----------------------------------

The basis set library of |molcas| contains an extensive set of basis sets both
segmented and generally contracted. The files in the basis directory are named
in upper case after the basis type
label (see below). Three sets of generally contracted basis sets have been
especially designed for |molcas|. They are based on the Atomic Natural Orbital
(ANO) concept and are labeled ANO-X (X=S, L, or RCC). They
have been designed to give a balanced description of the atoms in ground,
excited, and ionized states. A more detailed description of these basis sets is
given below. A fourth basis set, which is especially designed for the
calculation of electric properties of molecules (POL) will also be described.

In addition to this, an
subset of segmented standard basis sets are included, for example, STO-3G,
3-21G 4-31G, 6-31G, 6-31G*, 6-31G**, cc-pVXZ (X=D,T,Q), and aug-cc-pVXZ
(X=D,T). In addition, the library also contains different variants of the
Turbomole RI basis sets. For additional all electron basis set we recommend a
visit to the EMSL Basis Set Exchange
(https://bse.pnl.gov/bse/portal).
All basis sets are stored in the
directory :file:`basis_library`. The different types of available basis sets can be
found in the file :file:`basistype.tbl` in this directory. Aliases for the names are
listed in the file :file:`basis.tbl`. However, the best way to find out which basis sets
are available is to issue the command :command:`molcas help basis X` where :command:`X` is the atom.
Note that a short hand notation can be used for most basis sets: for example
ANO-L-VTZP will give a basis set of valence triple zeta accuracy with
polarization functions.

.. include:: basis_library/ANO-S.inc

.. include:: basis_library/ANO-L.inc

.. include:: basis_library/ANO-RCC.inc

.. 6.1 .. include:: basis_library/ANO-DK3.inc

.. include:: basis_library/POL.inc

.. include:: basis_library/AE_library_structure.inc

.. _UG\:sec\:the_ecp_libraries:

The ECP Library
---------------

|molcas| is able to perform *effective core potential* (ECP) calculations
and *embedded cluster* calculations.
In ECP calculations, only the *valence* electrons of a molecule are
explicitly handled in a quantum mechanical calculation,
at a time that the *core* electrons are kept frozen and are represented by
ECP's.
(An example of this is a calculation on :math:`\ce{HAt}` in which only the 5d, 6s and 6p
electrons of Astatine and the one of Hydrogen are explicitly considered.)
Similarly, in *embedded cluster* calculations,
only the electrons assigned to a piece of the whole system (the *cluster*)
are explicitly handled in the quantum mechanical calculation,
under the assumption that they are the only ones relevant for
some local properties under study;
the rest of the whole system (the *environment*)
is kept frozen and represented by embedding potentials
which act onto the *cluster*.
(As an example, calculations on a :math:`\ce{TlF12^{11-}}` cluster embedded in
a frozen lattice of :math:`\ce{KMgF3}` can be sufficient to calculate spectroscopical
properties of :math:`\ce{Tl^+}`\-doped :math:`\ce{KMgF3}` which are due to the :math:`\ce{Tl+}`
impurity.)

In order to be able to perform ECP calculations in molecules, as well as
*embedded cluster* calculations in ionic solids, with the
Ab Initio Model Potential method
(AIMP) :cite:`Huzinaga:87,Barandiaran:88,Barandiaran:90,Wittborn:95,Rakowitz:99a,Rakowitz:99b`
|molcas| is provided with the library
:file:`ECP` which includes nonrelativistic and relativistic
*core* ab initio model potentials and
*embedding* ab initio model potentials
representing both complete-cations and complete-anions in ionic
lattices :cite:`Barandiaran:88,Barandiaran:92`.

Before we continue we should comment a little bit on the terminology used here.
Strictly speaking, ECP methods are all that use the frozen-core approximation.
Among them, we can distinguish two families: the "pseudopotential" methods
and the "model potential" methods.
The pseudopotential methods are ultimately based on the
Phillips--Kleinman equation :cite:`Phillips:59`
and handle valence nodeless pseudo orbitals.
The model potential methods are based on the Huzinaga
equation :cite:`Huzinaga:71,Huzinaga:73`
and handle node-showing valence orbitals;
the AIMP method belongs to this family.
Here, when we use the general term ECP we will be referring to the more
particular of AIMP.
According to its characteristics,
the AIMP method can be also applied to represent
frozen-ions in ionic lattices in embedded cluster calculations;
in this case,
we will not be very strict in the nomenclature and we will also call ECP's
to the frozen-ion (embedding) *ab initio* model potentials.

The effective potentials in the libraries include the effects of the atomic
core wave functions (embedding ion wave functions) through the
following operators:

* a local representation of the core (ion) Coulomb operator,
* a non-local spectral representation of the core (ion) exchange operator,
* a core (ion) projection operator,
* a spectral representation of the relativistic mass-velocity and
  Darwin operators corresponding to the valence orbitals,
  if the Cowan--Griffin-based scalar relativistic
  CG-AIMP method :cite:`Barandiaran:90` is used.
* a spectral representation of the relativistic no-pair Douglas--Kroll
  operators, if the scalar relativistic no-pair Douglas--Kroll NP-AIMP method
  :cite:`Wittborn:95,Rakowitz:99a,Rakowitz:99b` is used.

Given the quality and non-parametric nature of the operators
listed above, the flexibility of the basis sets to be used
with the AIMP's is crucial, as in any *ab initio* method.

The valence basis sets included in the libraries
have been obtained by energy minimization in atomic valence-electron
calculations,
following standard optimization procedures.
All the experience gathered in the design of
molecular basis sets starting from all-electron atomic basis sets,
and in particular from segmented minimal ones,
is directly applicable to the AIMP valence basis sets included in
the libraries.
They are, for non-relativistic and relativistic Cowan--Griffin AIMPs, minimal
basis sets with added functions,
such as polarization and diffuse functions;
in consequence,
the minimal sets should be split in molecular calculations
in order to get reasonable sets (a splitting pattern
is recommended in the library for every set);
the splitting can be done by means of "the basis set label".
For the relativistic no-pair Douglas--Kroll AIMPs contracted valence basis sets
are given directly in a form which is recommended in molecular calculations,
i.e. they are of triple zeta quality in the outer shells and contain
polarization functions.
In both cases these *valence* basis sets contain very
*inner* primitive GTF's: They are necessary since,
typical to a model potential method,
the valence orbitals will show correct nodal structure.
Finally, it must be noted that
the core AIMP's can be safely mixed together with all-electron basis sets.

In AIMP *embedded cluster calculations*,
the cluster basis set,
which must be decided upon by the user,
should be designed following high quality standard procedures.
Very rigid cluster basis sets should not be used.
In particular,
the presence of the necessary embedding projection operators,
which prevent the cluster densities from collapsing onto the crystal lattice,
demands flexible cluster bases, including, eventually,
components outside the cluster volume :cite:`Pascual:93`.
The use of flexible cluster basis sets is then a
necessary requirement to avoid artificial frontier effects,
not ascribable to the AIMP embedding potentials.
This requirement is unavoidable, anyway, if good correlated wave
functions are to be calculated for the cluster.
Finally, one must remember that
the AIMP method does exclude any correlation between the
cluster electronic group
and the embedding crystal components; in other words, only
intra-cluster correlation effects can be accounted for in
AIMP embedded cluster calculations.
Therefore the cluster-environment partition
and the choice of the cluster wave function
must be done accordingly. In particular, the use of
one-atom clusters is not recommended.

Core- and embedding-AIMP's can be combined in a natural way
in valence-electron, embedded cluster calculations.
They can be used with any of the different types of wave
functions that can be calculated with |molcas|.

.. include:: basis_library/ECP.inc

.. include:: basis_library/ECP.nodeless.inc

.. include:: basis_library/ECP_library_structure.inc
