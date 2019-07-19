.. _sec\:about_this_manual:

The |molcas| Manual
===================

.. _sec\:who_should_read:

Manual in Four Parts
--------------------

This manual is designed for use with the *ab initio* Quantum
chemistry software package |molcas| |molcasversion| developed at the
by the world-wide |molcas| team where its base and origin is the
Department of Theoretical Chemistry, Lund University, Sweden. |molcas| is designed for use
by Theoretical Chemists and requires knowledge of the
Chemistry involved in the calculations in order to produce and
interpret the results correctly. The package can be moderately difficult to use
because of this "knowledge requirement", but the results are often more
meaningful than those produced by :bdit:`blackbox` packages which may not be
sufficiently chemically precise in either input or output.

The |molcas| manual is divided in four parts to facilitate its use.

#. The |molcas| :bdit:`Installation Guide` describes simple and more complex aspects on how to install, tailor, and
   control the |molcas| package.

#. The :bdit:`Short Guide` to |molcas| is a brief introductory guide which addresses the needs of the novice and intermediate users
   and is designed for all those who want to start using |molcas| as soon as possible.
   Only basic environment definitions, simple input examples, and minimal description of output results are included in the short guide.

   Two types of introductory tutorials are given in the short guide: problem-based and program-specific.

   #. Problem-based tutorials are exercises focused on solving a simple Quantum
      Chemical project and contain all the required input files. Examples include
      computing electronic energy of a molecule at different levels of
      theory, optimizing the geometry of a molecule, calculating the transition state in the ground
      state of a chemical system, and computing an excited state.
      The input files for this section can be found in the directory :file:`$MOLCAS/doc/samples/problem_based_tutorials`.
      These examples are also employed in |molcas| workshops that the |molcas| team has organized in recent years.

   #. Another type of tutorial is designed for the first-time user to provide an understanding of program modules
      contained in |molcas| include simple, easy-to-follow examples for many of these modules.

   The systems covered in the short guide are not necessarily calculated with most suitable methods or produce highly significant results,
   but provide both several tips for the beginner and actual input file formats.

   The :bdit:`Short Guide` to |molcas| can be independently printed as a booklet.

#. The |molcas| :bdit:`User's Guide` contains a complete listing of the input
   keywords for each of the program modules and a information regarding
   files used in each calculation. Here the user will find all keywords that can be
   used together with a specific program and thus how to set up the input for a
   |molcas| run.

#. :bdit:`Advanced Examples` and :bdit:`Annexes` include outlines of
   actual research performed using |molcas|.

   The approach to a research project is outlined including input files and shell scripts. More
   importantly, however, the value of the calculations is evaluated and
   advanced features of |molcasversion| are used and explained to improve the
   value of the results.

The complete manual is available on the net in HTML and PDF formats
(|MolcasWWW|).

.. _notation:

Notation
--------

For clarity, some words are printed using special typefaces.

* Keywords, i.e. words used in input files, are typeset in
  the small-caps typeface, for example :kword:`EndOfInput`.

* Programs (or modules) are typeset in the teletype typeface.
  This will eliminate some potential confusion. For example,
  when discussing the RASSCF method, regular uppercase letters
  are used, while the program will look like :program:`RASSCF`.

* Files are typeset in the slanted teletype typeface, like
  :file:`InpOrb`.

* Commands, unix or other, are typeset in a sans serif typeface,
  like :command:`ln -fs`.

* Complete examples, like input files, shell scripts, etc,
  are typeset in the teletype typeface.
