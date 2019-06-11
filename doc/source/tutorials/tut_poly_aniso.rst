.. index::
   single: Program; Poly_aniso
   single: Poly_aniso

.. _TUT\:sec\:poly_aniso:

:program:`POLY_ANISO` --- Semi-*ab initio* Electronic Structure and Magnetism of Polynuclear Complexes Program
===============================================================================================================

The program :program:`POLY_ANISO` calculates nonperturbatively the temperature- and field- dependent magnetic
properties (Van Vleck susceptibility tensor and function, molar magnetization vector and function) and the
pseudospin Hamiltonians for Zeeman interaction (the :math:`g` tensor and higher rank tensorial components) and the
zero-field splitting (the :math:`D` tensor and higher rank tensorial components) for arbitrary mononuclear complexes
and fragments on the basis of ab initio spin-orbit calculations.
:program:`POLY_ANISO` requires as input file the :file:`RUNFILE` containing all necessary
*ab initio* information: spin orbit eigenstates, angular momentum matrix elements, the states been mixed
by the spin-orbit coupling in :program:`RASSI`, etc. Usually, the :program:`POLY_ANISO`
runs after :program:`RASSI`.

For a proper spin-orbit calculation the relativistic basis sets should be used for the whole calcualtion.
For :program:`SEWARD`, the atomic mean-field (:kword:`AMFI`), Douglas--Kroll (:kword:`DOUG`) must be employed.
To ensure the computation of angular momentum integrals the :kword:`ANGMOM` should be also used, specifying the origin
of angular momentum integrals as the coordinates of the magnetic center of the molecule, i.e. the coordinates of the atom
where the unpaired electrons mainly reside. For program :program:`RASSI` the necessary keywords are: :kword:`SPIN`,
since we need a spin-orbit coupling calculation, and :kword:`MEES`, to ensure the computation of angular momentum
matrix elements in the basis of spin-free states (SFS).

.. index::
   single: Poly_aniso; Input

In the cases where spin-orbit coupling has a minor effect on the low-lying energy spectrum (most of the
isotropic cases: :math:`\ce{Cr^{3+}}`, :math:`\ce{Gd^{3+}}`, etc.) the pseudospin is usually the same as the ground spin. For these cases
the :program:`POLY_ANISO` may run without specifying any keywords in the input file.

::

  &POLY_ANISO

In the cases when spin-orbit coupling play an important role in the low-lying energy spectrum, i.e. in the cases of e.g. octahedral :math:`\ce{Co^{2+}}`,
most of the lanthanide complexes, the pseudospin differs strongly from the spin of the ground state. In these cases,
the dimension of the pseudospin can be found by analysing the spin-orbit energy spectrum obtained at :program:`RASSI`.
The pseudospin is best defined as a group of spin-orbit states close in energy. Once specified, these eigenstates are further used
by the :program:`POLY_ANISO` to build proper pseudospin eigenfunctions. As an example of an input for :program:`POLY_ANISO`
requiring the computation of all magnetic properties (which is the default) and the computation of the :math:`g` tensor for the ground
Kramers doublet (i.e. pseudospin of a Kramers doublet is :math:`\tilde{S}=1/2`).

::

  &POLY_ANISO
   MLTP
   1
   2

.. compound::

  :program:`POLY_ANISO` has implemented pseudospins: :math:`\tilde{S}=1/2`, :math:`\tilde{S}=1`, ..., up to :math:`\tilde{S}=7/2`. The user can also ask for more pseudospins at the same time: ::

    &POLY_ANISO
     MLTP
     3
     2 4 2

  For the above input example, the :program:`POLY_ANISO` will compute the :math:`g` tensor for the ground Kramers doublet
  (spin-orbit states 1 and 2), the :math:`g` tensor, ZFS tensor and coefficients of higher rank ITO for the pseudospin
  :math:`\tilde{S}=3/2` (spin orbit functions 3--6), and the :math:`g` tensor for the third excited Kramers doublet (spin orbit functions 7 and 8).

.. index::
   single: Poly_aniso; Output

:program:`POLY_ANISO` Output
----------------------------

The :program:`POLY_ANISO` section of the |molcas| output is divided in four parts. In the first part, the :math:`g` tensor and higher rank Zeeman tensors are computed. They are followed by :math:`D` tensor and higher rank ZFS tensors. The program also computes the angular moments in the direction of the main magnetic axes.

In the second part, the paramaters of the crystal field acting on the ground atomic multiplet of lanthanides are calculated.

In the third part, the powder magnetic susceptibility is printed, followed by the magnetic susceptibility tensors with and without intermolecular interaction included.

In the fourth part, magnetization vectors (if required) are printed, and then the powder molar magnetization calculated for the :kword:`TMAG`
temperature.

The keywords :kword:`TINT` and :kword:`HINT` control the temperature and field intervals for computation of
magnetic susceptibility and molar magnetization respectively.
Computation of the magnetic properties at the experimental temperature and field points with the estimation of the standard deviation from experiment
is also possible via :kword:`TEXP`, defining the experimental temperature and measured magnetic susceptibility and
:kword:`HEXP`, defining the experimental field and averaged molar magnetization.

::

  &POLY_ANISO
  TITLE
  g tensor and magnetic susceptibility
  TYPE
  4
  MLTP
  2
  3 3
  TINT
  0.0 100 101 0.001

The above input requires computation of the parameters of two pseudospins :math:`\tilde{S}=1`: the ground (spin-orbit functions 1--3)
and first excited (spin-orbit functions 4--6) and the magnetic susceptibility in 101 steps equally distributed in
the temperature domain 0.0--100.0 K.

:program:`POLY_ANISO` --- Basic and Most Common Keywords
--------------------------------------------------------

.. class:: keywordlist

:kword:`MLTP`
  Specifies the number and dimension of the pseudospins Hamiltonians

:kword:`TMAG`
  Sets the temperature for the computation of molar magnetization

:kword:`MVEC`
  Number and radial coordinates of directions for which the magnetization vector will be computed
