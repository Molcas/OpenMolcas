.. _UG\:sec\:structure_of_the_ae_basis_set_libraries:

Structure of the all electron basis set library
...............................................

.. compound::

  The start of a given basis set in the library is given by the line ::

    /label

  where "label" is the basis set label, as defined below in the input
  description to
  :program:`SEWARD`. Then follows two lines with the
  appropriate literature reference for that basis set. These
  cards are read by
  :program:`SEWARD`
  and must thus be included in the library, and may not be blank.
  Next is a set of comment lines, which begin with an asterisk in column 1,
  giving some details of the basis
  sets. A number of lines follow, which specifies the basis set:

#. Charge of the atom and the highest angular momentum.
   For each angular momentum (l) then follows.

   .. _it\:ae-2:

#. Number of primitives and contracted functions for angular momentum l
   (must be identical to those given in the basis set label).

   .. _it\:ae-3:

#. Exponents of the primitive functions.

   .. _it\:ae-4:

#. The contraction matrix (with one CGTO per column). Note that all
   basis sets are given in the generally contracted format, even if they
   happen to be segmented. Note that the number of CGTOs must correspond
   to the data given in the label.

The following is an example of an entry in a basis set library. ::

  * This is the Huzinaga 5s,2p set contracted to 3s,2p        -- Comment
  * according to the Dunning paper.                           -- Comment
  /H.TZ2P.Dunning.5s2p.3s2p.                                  -- Label
  Exponents  : S. Huzinaga, J. Chem. Phys., 42, 1293(1965).   -- First ref line
  Coefficients: T. H. Dunning, J. Chem. Phys., 55, 716(1971). -- Second ref line
   1.0 1                                                      -- Charge, sp
   5 3                                                        -- 5s->3s
   52.56 7.903 1.792 0.502 0.158                              -- s-exponents
   0.025374  0.0  0.0                                       -- contr. matrix
   0.189684  0.0  0.0                                       -- contr. matrix
   0.852933  0.0  0.0                                       -- contr. matrix
   0.0     1.0  0.0                                         -- contr. matrix
   0.0     0.0  1.0                                         -- contr. matrix
   2  2                                                       -- 2p->2p
   1.5  0.5                                                   -- p-exponents
   1.0 0.0                                                  -- contr. matrix
   0.0 1.0                                                  -- contr. matrix
