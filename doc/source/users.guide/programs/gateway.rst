.. index::
   single: Program; Gateway
   single: Gateway

.. _UG\:sec\:gateway:

:program:`gateway`
==================

.. only:: html

  .. contents::
     :local:
     :backlinks: none

.. xmldoc:: <MODULE NAME="GATEWAY">
            %%Description:
            <HELP>
            The Gateway module collects information about molecular
            system (geometry, basis sets, symmetry) to be used for future calculations.
            Note that there are two input styles - the old style format and the new XYZ format.
            The input can also contain an embedded reaction field input (starts with the
            "RF-input" keyword and terminates with the "End of RF-input" keyword). Keywords below
            carrying a "(RF)" are associated with that embedded RF-input section.
            Gateway also controls options associated with auxiliary basis sets to be used
            in density fitting procedures.
            </HELP>

The Gateway module collects information about molecular
system (geometry, basis sets, symmetry) to be used for future calculations.

Gateway module is a subset of :program:`seward`. All keywords
for this module can also appear as an input for :program:`SEWARD`, however,
for clarity the information about molecular system can be placed
as an input for this module. Note, that :program:`gateway` does not
compute any integral, and so must be followed by run of :program:`SEWARD`
module.

:program:`GATEWAY` destroys the communication file :file:`RUNFILE`,
if it is used in a combination with geometry optimization it should run
outside the optimization loop.

Input
-----

This sections will describe the various possible input blocks in :program:`Gateway`.
These control

* the molecular structure (coordinates, symmetry and basis sets),
* explicit auxiliary basis sets in terms of CD basis sets (aCD and acCD) or
  external auxiliary basis sets,
* parameters for reaction field calculations, i.e. parameters for the Kirkwood model
  or the PCM model and options for Pauli repulsion integral and external field integrals,
* options for finite nuclear charge distribution models in association with relativistic calculations, and
* the option to use the Saddle method to locate transitions state geometries.

The :program:`Gateway` input section always starts with the program reference: ::

  &GATEWAY

General keywords
................

.. class:: keywordlist

:kword:`TITLE`
  The keyword followed by a title.

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="TITLE" APPEAR="Title" KIND="CUSTOM" LEVEL="BASIC">
              %%Keyword: TITLE <basic>
              <HELP>
              The keyword followed by a title
              </HELP>
              </KEYWORD>

:kword:`TEST`
  :program:`GATEWAY` will only process the input and generate a non-zero return code.

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="TEST" APPEAR="Test" KIND="SINGLE" LEVEL="BASIC">
              %%Keyword: Test <basic>
              <HELP>
              GATEWAY will only process the input and generate a non-zero
              return code.
              </HELP>
              </KEYWORD>

:kword:`EXPErt`
  Activates "expert mode", in which various default settings are
  altered. This will, for example, allow the user to combine
  relativistic and non-relativistic basis sets.

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="EXPE" APPEAR="Expert mode" KIND="SINGLE" LEVEL="ADVANCED">
              %%Keyword: Expert <advanced>
              <HELP>
              Activates "expert mode", in which various default settings are
              altered. This will, for example, allow the user to combine
              relativistic and non-relativistic basis sets.
              </HELP>
              </KEYWORD>

:kword:`BASDIR`
  The keyword allows to set up an extra location for basis set files.
  The value can be either an absolute path (started from /) or relative to
  submit directory, e.g. BASDIR=.
  In order to use a local copy of a basis set file with name FOO --- place
  this file into directory specified in BASDIR

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="BASDIR" APPEAR="BasDir" KIND="STRING" LEVEL="BASIC">
              %%Keyword: BASDIR <basic>
              <HELP>
              The keyword allows to set up an extra location for basis set files.
              The value can be either an absolute path (started from /) or relative to
              submit directory, e.g. BASDIR=.
              In order to use a local copy of a basis set file with name FOO - place
              this file into directory specified in BASDIR
              </HELP>
              </KEYWORD>

:kword:`BASLIB`
  The keyword followed by the absolute path to the basis set library directory. The default
  is the :file:`$MOLCAS/basis_library` directory. Note that this directory must also be host to
  local copies of the .tbl files.

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="BASLIB" APPEAR="BasLib" KIND="STRING" LEVEL="BASIC">
              %%Keyword: BASLIB <basic>
              <HELP>
              The keyword followed by the absolute path to the basis set library directory. The default
              is the $MOLCAS/basis_library directory. Note that this directory must also be host to
              local copies of the .tbl files.
              </HELP>
              </KEYWORD>

:kword:`RTRN`
  Max number of atoms for which bond lengths, angles and dihedral
  angles are listed, and
  the radius defining the maximum length of a bond follows on
  the next line. The latter is used as a threshold when printing out
  angles and dihedral angles. The length can be followed by
  :kword:`Bohr` or
  :kword:`Angstrom` which indicates the unit in which the length
  was specified, the default is
  :kword:`Bohr`.
  The default values are 15 and 3.0 au.

  .. xmldoc:: %%Keyword: RTRN <advanced>
              Max number of atoms for which bond lengths, angles and dihedral
              angles are listed, and
              the radius defining the maximum length of a bond follows on
              the next line. The latter is used as a threshold when printing out
              angles and dihedral angles. The length can be followed by
              "Bohr" or "Angstrom" which indicates the unit in which the length
              was specified, the default is "Bohr".
              The default values are 15 and 3.0 au.

:kword:`ISOTopes`
  Specify isotopic substitutions or atomic masses. By default, the mass of the most
  abundant or stable isotope is used for each atom. With this keyword different
  isotopes or arbitrary masses can be chosen. The keyword is followed by the number
  :math:`n` of isotopic specifications, and then by :math:`n` lines. Each of these
  lines should contain the symmetry-unique index of the atom for which the default
  mass is to be modified and either the mass number of the desired isotope (tabulated
  masses for most known isotopes are available in the code, use ``0`` for the default
  isotope) or the desired mass in dalton, in the latter case the keyword :kword:`Dalton`
  should follow. Note that all atoms belonging to the same "center type" must have the
  same mass. This usually means all atoms of a given element with the same basis set.
  If more fine-grained specifications are wanted, additional center types must be
  created by using several :kword:`BASIs Set` blocks for the same element (native
  input) or by using labels (XYZ input).

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="ISOTOPES" APPEAR="Isotopic specification" KIND="CUSTOM" LEVEL="ADVANCED">
              %%Keyword: Isotopes <advanced>
              <HELP>
              Specifies isotopes or masses. First write the number of atom masses to change,
              then that number of lines, on each: the symmetry-unique index of the atom and
              (a) the mass number of the isotope, or (b) the mass in dalton and the word DALTON.
              </HELP>
              </KEYWORD>

:kword:`ECPShow`
  Force :program:`GATEWAY` to print ECP parameters.

  .. xmldoc:: <GROUP MODULE="GATEWAY" KIND="BOX" NAME="PROPT" APPEAR="Print options" LEVEL="BASIC">

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="ECPS" APPEAR="Print ECP info" KIND="SINGLE" LEVEL="BASIC">
              %%Keyword: ECPSHOW <basic>
              <HELP>
              Force GATEWAY to print ECP parameters.
              </HELP>
              </KEYWORD>

:kword:`AUXShow`
  Force :program:`GATEWAY` to print auxiliary basis set parameters.

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="AUXS" APPEAR="Print auxiliary basis info" KIND="SINGLE" LEVEL="BASIC">
              %%Keyword: AUXSHOW <basic>
              <HELP>
              Force GATEWAY to print auxiliary basis set parameters.
              </HELP>
              </KEYWORD>

:kword:`BSSHow`
  Force :program:`GATEWAY` to print basis set parameters.

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="BSSH" APPEAR="Print basis info" KIND="SINGLE" LEVEL="BASIC">
              %%Keyword: BSSHOW <basic>
              <HELP>
              Force GATEWAY to print basis set parameters.
              </HELP>
              </KEYWORD>

:kword:`VERBose`
  Force :program:`GATEWAY` to print a bit more verbose.

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="VERB" APPEAR="Verbose output" KIND="SINGLE" LEVEL="BASIC">
              %%Keyword: Verbose <basic>
              <HELP>
              Force GATEWAY to print a bit more verbose.
              </HELP>
              </KEYWORD>

  .. xmldoc:: </GROUP>

Molecular structure: coordinates, symmetry and basis sets
.........................................................

There are three different ways to specify the molecular structure, symmetry and
the basis sets in :program:`Gateway`:

* XYZ input,
* the so-called native input (old |molcas| standard).

.. * XYZ input, and
   * Z-matrix input.

Note that only XYZ input for :program:`Gateway` is supported by Graphical User interface.
:program:`Gateway` makes a decision about the type of the input based on keywords.
If :kword:`Coord` is used, it assumes that the input is in XYZ format.

.. , if :kword:`ZMAT` is used,
   it assumes Z-matrix input.

The three different modes will be described below.

Z-matrix and XYZ input
::::::::::::::::::::::

Some times it is more convenient to set up information about coordinates in
a standard form of Z-matrix or Cartesian coordinates. In this case,
the basis set for the atoms should be specified after the :kword:`XBAS`
keyword. After that either :kword:`ZMAT` or :kword:`XYZ` should appear
to specify the coordinates.
Note that coordinates in these formats use ångström as units.

.. class:: keywordlist

:kword:`XBAS`
  A keyword to specify the basis for atoms. The specification is very similar
  to the native format: ``ATOM.BasisSet``. Each new atom is written at a new line.
  The end of the keyword is marked by an :kword:`End of basis` line.

  If all atoms have the same basis, e.g. ANO-S-VDZ, it is possible to use
  this name without element name. In this case there is no need to specify
  :kword:`End of basis`.

  .. compound::

    Example: ::

      XBAS=STO-3G

    or ::

      XBAS
      C.STO-3G
      H.STO-3G
      End of basis

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="XBAS" APPEAR="Basis set (alternate format)" KIND="CUSTOM" LEVEL="ADVANCED">
              %%Keyword: XBAS <basic>
              <HELP>
              A keyword to specify the basis for atoms. The specification is very similar
              to the native format: ATOM.BasisSet. Each new atom is written at a new line.
              The end of the keyword is marked by an 'End of basis' line.
              </HELP>

              If all atoms have the same basis, e.g. ANO-S-VDZ, it is possible to use
              this name without element name. In this case there is no need to specify
              'End of basis'.
              </KEYWORD>

:kword:`ZMAT`
  Alternative format to give coordinates in terms of bond lengths,
  bond angles, and dihedral angles.
  Each line of a Z-matrix gives the
  internal coordinates for one of the atoms within the molecule with the following
  syntax:

  .. container:: list

    **Name  I bond-length  J bond-angle  K dihedral-angle**

    **Name** is the label (atomic symbol + string) for a symmetry distinct center L;

    **I bond-length** distance of L from atom I;

    **J bond-angle** planar angle between atoms L--I--J;

    **K dihedral-angle** dihedral angle between atoms L--I--J--K.

  Note that the first atom only requires the **Name** and defines the origin of
  Cartesian axis.
  The second atom requires **Name  I bond-length** and it will define the Z axis.
  The third atom requires **Name  I bond-length  J bond-angle** and defines the
  XZ plane (and implicitly, the Y axis).

  Only numerical values must be used (no variable names) and ångströms
  and degrees are assumed as units.

  Two types of special atoms are allowed: *dummy* **X** atoms and
  *ghost* **Z** atoms. The former will appear in the calculations,
  they have a nuclear charge of 0 and have not electrons and Basis Set.
  They will also appear in the definition of internal coordinates in :program:`SLAPAF`.
  The latter are used only within the Z-Matrix definition of the geometry but
  they will appear in the final Z-matrix section in :program:`SLAPAF`.
  Both special atoms can be used to define the Cartesian axis and the symmetry elements.

  **End of ZMAT** or a blank line mark the end of the section.

  Here is an example for (S)-1-chloroethanol (:math:`C_1` symmetry): ::

    XBAS
    H.ANO-L...2s1p.
    C.ANO-L...3s2p1d.
    O.ANO-L...3s2p1d.
    Cl.ECP.Huzinaga.7s7p1d.1s2p1d.7e-NR-AIMP.
    End of basis
    ZMAT
    C1
    O2      1   1.40000
    C3      1   1.45000   2   109.471
    H4      1   1.08900   2   109.471     3   120.000
    Cl5     1   1.75000   2   109.471     3  -120.000
    H6      2   0.94700   1   109.471     3   180.000
    H7      3   1.08900   1   109.471     2   180.000
    H8      3   1.08900   1   109.471     7   120.000
    H9      3   1.08900   1   109.471     7   240.000
    End of z-matrix

  In geometry optimizations, :program:`SLAPAF` will regenerate the coordinates as
  Z-matrix in the section with the summary concerning each iteration. This will
  be possible only if *ghost* atoms are used within the first three atoms or
  if they are not used at all.

  Both :kword:`BASIs` and :kword:`ZMAT` keywords can be used at the same time. Here is an example
  for a complex between methanol and water (:math:`C_s` symmetry): ::

    Symmetry
     Y
    XBAS
    H.ANO-L...1s.
    C.ANO-L...2s1p.
    O.ANO-L...2s1p.
    End of basis
    ZMAT
    C1
    O2  1 1.3350
    H3  1 1.0890  2 109.471
    H4  1 1.0890  2 109.471  3 -120.
    H6  2 1.0890  1 109.471  3  180.
    End of z-matrix
    Basis set
    O.ANO-L...2s1p.
     O    -2.828427     0.000000     2.335000  / Angstrom
    End of basis
    Basis set
    H.ANO-L...1s.
     H    -2.748759     0.819593     2.808729  / Angstrom
    End of basis

  In this case :program:`SLAPAF` will not regenerate the Z-matrix.

  .. xmldoc:: %%Keyword: ZMAT <basic>
              Alternative format to give coordinates in the form of Z-matrix.
              Only numerical values must be used (no variable names) and angstroms
              and degrees are assumed as units. Special ghost Z and dummy X atoms
              are allowed. 'End of ZMAT' or a blank line marks the end of the section.

:kword:`XYZ`
  The keyword is followed by XYZ formatted file (a reference to a file),
  or file, inlined into the input.

  .. compound::

    Example: ::

      XBAS=STO-3G
      XYZ=$CurrDir/Water.xyz

    or ::

      XBAS=STO-3G
      XYZ
      1
       note Angstrom units!
      C 0 0 0

  Currently, the :kword:`XYZ` keyword does not operate with symmetry, and
  the calculation is always performed without symmetry.

  .. xmldoc:: %%Keyword: XYZ <basic>
              Alternative format to set up geometry as XYZ formatted file

Advanced XYZ input
::::::::::::::::::

If the geometry is specified in XYZ format, all atoms should be specified.
The default units are ångströms. By default, maximum possible symmetry is used.

"Molcas XYZ" file format is an extension of plain XYZ format.

* First line of this file contains the number of atoms.

* Second line (a comment line) can contain "a.u." or "bohr" to
  use atomic units, instead of default ångströms.
  Also this line can contain keyword TRANS, followed by 3 numbers,
  and/or ROT, followed by 9 numbers (in this case coordinates
  will be Translated by specified vector, and/or Rotated), and SCALE (or
  SCALEX, SCALEY, SCALEZ) followed by a scale factor.

* Remaining lines are used to specify Element and cartesian
  coordinates.

  Element name might be optionally followed by a Number (e.g. ``H7``),
  a Label (separated by ``_`` sign: e.g. ``H_INNER``), or Basis Set (separated by ``.``,
  e.g. ``H.STO-3G``)

.. class:: keywordlist

:kword:`COORD`
  The keyword followed on the next line by the name of an HDF5 (created by any module), or the name of an XYZ file,
  or inline coordinates in XYZ format. If the file is located in the same directory, where
  |molcas| job was submitted there is no need to specify the PATH to this file.
  The keyword may appear several times. In this case all coordinate files
  will be concatenated, and considered as individual fragments.

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="COORD" APPEAR="Coord" KIND="CUSTOM" INPUT="REQUIRED" LEVEL="BASIC">
              %%Keyword: COORD (XYZ format) <basic>
              <HELP>
              The keyword followed on the next line by the name of an HDF5 or XYZ file,
              or inline coordinates in XYZ format.
              The keyword may appear several times. In this case all coordinate files
              will be concatenated, and considered as individual fragments.
              </HELP>
              </KEYWORD>

:kword:`BASIS`
  The keyword can be used to specify global basis set for all atoms, or for a group of atoms.
  The keyword followed by a label of basis set, or by comma separated list of basis sets for
  individual atoms.

  Note! The basis set definition in XYZ mode does not allow to use
  inline basis set.

  Example: ::

    COORD
    4

    C           0.00000 0.00000 0.00000
    H           1.00000 0.00000 0.00000
    H           0.00000 1.00000 0.00000
    H           0.00000 0.00000 1.00000
    BASIS
    STO-3G, H.6-31G*

  In this example, the C atom (in the origin) will have the basis set STO-3G and
  the H atoms 6-31G*.

  If keyword BASIS never appears in the input, the default basis,
  ANO-S-MB, will be used.

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="BASIS (XYZ)" APPEAR="Basis set" KIND="STRING" LEVEL="BASIC">
              %%Keyword: BASIS (XYZ format) <basic>
              <HELP>
              The keyword followed on the next line by the name of global basis set for
              all atoms, or by comma separated list of basis sets for individual atoms.
              Note! The basis set definition in XYZ mode does not allow to use
              inline basis set.
              </HELP>
              </KEYWORD>

:kword:`GROUP`
  The keyword can be used to specify the symmetry of the molecule.

  The keyword must be followed by one of:

  * FULL (default) --- use maximum possible subgroup of :math:`D_{2h}`
  * NOSYM (same as E, or C1)
  * space separated list of generators: e.g. X XY (for more details see SYMMETRY keyword)

  .. Limitations: in the current implementation atom labels, and basis sets are ignored
     during symmetry recognition.

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="GROUP" APPEAR="Group" KIND="STRING" LEVEL="BASIC">
              %%Keyword: GROUP (XYZ format) <basic>
              <HELP>
              The keyword followed on the next line by the list of group generators
              (with the same syntax as SYMMETRY keyword),
              or by FULL (highest possible group), or by NOSYM, if no symmetry operations
              should be used. The keyword can be used only with XYZ format of input,
              after COORD keyword.
              </HELP>
              </KEYWORD>

If XYZ input has been used in :program:`gateway`, a file with native |molcas| input will be
produced and stored in working directory under the name :file:`findsym.std`.

Note that choosing XYZ input you are expecting that the coordinates might be transformed.
It can be shown by the following example: ::

  &gateway
  coord
  3

  O 0 0 0
  H 1.0000001 0 0
  H 0 1 0.0000001
  *nomove
  *group=c1

The geometry of the molecule is slightly distorted, but within a threshold it is :math:`C_{2v}`.
Thus by default (keywords :kword:`nomove` and :kword:`group` are not active), the
coordinates will be transformed to maintain the highest possible symmetry.
If keyword :kword:`nomove` is active, the molecule is not allowed to rotate, and
a mirror plane :math:`xy` is the only symmetry element. Since the third hydrogen atom is
slightly out of this plane, it will be corrected. Only activation of the keyword :kword:`group=C1`
will ensure that the geometry is unchanged.

Native input
::::::::::::

If the geometry is specified in a native |molcas| format, only symmetry
inequivalent atoms should be specified. The default units are atomic units.
By default, symmetry is not used in the calculation.

.. class:: keywordlist

:kword:`SYMMetry`
  Symmetry specification follows on next line. There may be up to
  three different point group generators specified on that line. The
  generators of a point group is the minimal set of symmetry operators
  which is needed to generate all symmetry
  operators of a specific point group. A generator is in the input
  represented as a sequence of up to three of the characters x, y, and
  z. The order within a given sequence is arbitrary and the generators
  can be given in any sequence. Observe that the order of the irreps
  is defined by the order of the generators as
  (:math:`E`, :math:`g_1`, :math:`g_2`, :math:`g_1g_2`, :math:`g_3`, :math:`g_1g_3`, :math:`g_2g_3`,
  :math:`g_1g_2g_3`)! Note that :math:`E` is always assumed and should never
  be specified.

  Below is listed the possible generators.

  * **x** --- Reflection in the :math:`yz`-plane.
  * **y** --- Reflection in the :math:`xz`-plane.
  * **z** --- Reflection in the :math:`xy`-plane.
  * **xy** --- Twofold rotation around the :math:`z`-axis.
  * **xz** --- Twofold rotation around the :math:`y`-axis.
  * **yz** --- Twofold rotation around the :math:`x`-axis.
  * **xyz** --- Inversion through the origin.

  The default is no symmetry.

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="SYMMETRY" APPEAR="Symmetry" KIND="STRING" LEVEL="BASIC" EXCLUSIVE="COORD">
              %%Keyword: Symmetry (non-XYZ format) <basic>
              Symmetry point group is specified by up to three group generators.
              Possible generators are "x", "y", "z", "xy", "xz", "yz", and "xyz".
              The order of the irreps depends on the order of the generators.
              The keyword can be used only in 'native' input format.
              </KEYWORD>

:kword:`BASIs Set`
  This notes the start of a basis set definition.
  The next line always contains a basis set label.
  The basis set definition is alway terminated with the "End of Basis" keyword.
  For the definitions of basis set labels see the subsequent sections.
  Below follows a description of the options associated with the
  basis set definition.

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="BASIS (NATIVE)" APPEAR="Basis set" KIND="CUSTOM" LEVEL="ADVANCED" INPUT="REQUIRED" EXCLUSIVE="COORD">
              %%Keyword: BASIS (non-XYZ format) <basic>
              This notes the start of a basis set definition.
              The next line always contains a basis set label.
              The basis set definition is alway terminated with the "End of Basis" keyword.
              For details consult the manual.
              </KEYWORD>

  * **Label [/ option]** ---
    The label is a specification of a specific basis set, e.g.
    C.ANO...4s3p2d., which is an ANO basis set.
    If no option is specified
    :program:`GATEWAY` will look for the basis
    set in the default basis directory. If an option is specified it
    could either be the name of an alternative basis directory or
    the wording "Inline" which defines
    that the basis set will follow in the current input
    file. For the format of the
    **Inline** option see the section
    "Basis set format". Observe that the label is arbitrary for this
    option and will not be decoded.
    The **Label** card is mandatory.

  * **Name x, y, z (Angstrom or Bohr)** ---
    This card specifies an arbitrary (see next sentence!) name
    for a symmetry distinct center and its Cartesian coordinates.
    Observe, that the
    name "DBAS" is restricted to assign the center of the
    diffuse basis functions required to model the continuum
    orbitals in R-matrix calculations.
    The label is truncated to four characters. Observe that this
    label must be unique to each center. The coordinate unit can
    be specified as an option. The default unit is Bohr.
    There should at least be one card of this type in a basis set
    definition.

  * **Charge** ---
    The real entry on the subsequent line defines
    the charge associated with
    this basis set. This will override the default which is defined in
    the basis set library. The option can be used to put in ghost
    orbitals as well as to augment the basis sets of the library.
    The **Charge** card is optional.

    .. xmldoc:: %%Keyword: Charge (non-XYZ format) <advanced>
                The real entry on the subsequent line defines
                the charge associated with
                this basis set. This will override the default which is defined in
                the basis set library. The option can be used to put in ghost
                orbitals as well as to augment the basis sets of the library.
                The "Charger" card is optional.

  * **Spherical** [option] ---
    Specifying which shells will be in real spherical Gaussians. Valid options
    are "all" or a list of the shell characters separated by a blank. The
    shell characters are s, p, d, f, etc. All shells after p are by
    default in real spherical Gaussians, except for the d-functions in the
    6-31G family of basis sets which are in Cartesian.
    The **Spherical** card is optional. The s and p shells and the d-functions of
    the 6-31G family of basis sets are by default in Cartesian Gaussians.

    .. xmldoc:: %%Keyword: Spherical (non-XYZ format) <advanced>
                Specifying which shells will be in real spherical Gaussians. Valid options
                are "all" or a list of the shell characters separated by a blank. The
                shell characters are s, p, d, f, etc. All shells after p are by
                default in real spherical Gaussians, except for the d-functions in the
                6-31G family of basis sets which are in Cartesian.
                The "Spherical" card is optional. The s and p shells and the d-functions of
                the 6-31G family of basis sets are by default in Cartesian Gaussians.

  * **Cartesian** [option] ---
    Specifying which shells will be in a Cartesian Gaussian representation. For syntax
    consult the corresponding **Spherical** keyword.

    .. xmldoc:: %%Keyword: Cartesian (non-XYZ format) <advanced>
                Specifying which shells will be in a Cartesian Gaussian representation. For syntax
                consult the corresponding Spherical keyword.

  * **Contaminant** [option] ---
    Specifying for which shells the contaminant will be kept.
    The contaminants are functions of lower rank which are generated
    when a Cartesian shell is transformed to a spherical representation
    (e.g. :math:`r^2=x^2+y^2+z^2` for d-shells, p contaminants for f-shells,
    s and d contaminants for g-shells, etc.).
    Valid options are the same as for the **Spherical** keyword.
    The default is no contaminant in any shell. The **Contaminant** card is optional.

    .. xmldoc:: %%Keyword: Contaminant (non-XYZ format) <advanced>
                Specifying for which shells the contaminant will be kept.
                The contaminants are functions of lower rank which are generated
                when a Cartesian shell is transformed to a spherical representation
                (e.g. r^2=x^2+y^2+z^2 for d-shells, p contaminants for f-shells,
                s and d contaminants for g-shells, etc.).
                Valid options are the same as for the Spherical keyword.
                The default is no contaminant in any shell. The "Contaminant" card is optional.

  * **Muon** ---
    Specifying that the basis set is muonic.

    .. xmldoc:: %%Keyword: Muon (non-XYZ format) <advanced>
                Specifying that the basis set is muonic.

  * **End of Basis set** ---
    Marks the end of the basis set specification.
    This card is mandatory.

    .. xmldoc:: %%Keyword: End of Basis set (non-XYZ format) <advanced>
                Marks the end of the basis set specification.
                This card is mandatory.

Example of an input in native |molcas| format: ::

  &GATEWAY
  Title
  formaldehyde
  SYMMETRY
  X Y
  Basis set
  H.STO-3G....
  H1           0.000000    0.924258   -1.100293 /Angstrom
  End of basis

  Basis set
  C.STO-3G....
  C3           0.000000    0.000000   -0.519589 /Angstrom
  End of basis

  Basis set
  O.STO-3G....
  O            0.000000    0.000000    0.664765 /Angstrom
  End of basis

  End of input

Advanced keywords:

.. class:: keywordlist

:kword:`SYMThreshold`
  followed by a real number --- threshold for symmetry recognition (default is 0.01 Å)

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="SYMT" APPEAR="Symmetry Thr" KIND="REAL" LEVEL="ADVANCED" DEFAULT_VALUE="0.01" REALTIME_UPDATE="YES">
              %%Keyword: SYMThreshold (XYZ format) <advanced>
              <HELP>
              The keyword followed on the next line by the threshold for symmetry recognition code (default is 0.01)
              </HELP>
              </KEYWORD>

:kword:`MOVE`
  allow to translate and rotate molecule in order to find highest possible symmetry.
  (this is a default for all groups, except of :math:`C_1`)

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="MOVE" APPEAR="MOVE" KIND="SINGLE" LEVEL="ADVANCED" EXCLUSIVE="NOMOVE">
              %%Keyword: MOVE (XYZ format) <advanced>
              <HELP>
              Allow to translate and rotate molecule in order to find highest possible symmetry.
              (this is a default for all groups, except of C1)
              </HELP>
              </KEYWORD>

:kword:`NOMOVE`
  do not allow to transform coordinates while searching for highest group (default for :math:`C_1` group)

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="NOMOVE" APPEAR="NoMOVE" KIND="SINGLE" LEVEL="ADVANCED" EXCLUSIVE="MOVE">
              %%Keyword: NOMOVE (XYZ format) <advanced>
              <HELP>
              Do not allow to transform coordinates while searching for highest group (default for C1 group)
              </HELP>
              </KEYWORD>

:kword:`BSSE`
  followed by an integer. Indicates which XYZ-file that should be
  treated like ghost atoms.

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="BSSE" APPEAR="BSSE" KIND="INT" LEVEL="ADVANCED">
              %%Keyword: BSSE (XYZ format) <advanced>
              <HELP>
              Followed by an integer. Indicates which xyz-file that should be treated like ghost atoms.
              </HELP>
              </KEYWORD>

:kword:`VART`
  Specifies that the energy should not be considered invariant to translations.
  Translational variance is detected automatically, but sometimes it may be useful to enforce it.

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="VART" APPEAR="Var Trans" KIND="SINGLE" LEVEL="ADVANCED">
              %%Keyword: VarT <advanced>
              <HELP>
              Specifies that the energy should not be considered invariant to translations.
              Translational variance is detected automatically, but sometimes it may be useful to enforce it.
              </HELP>
              </KEYWORD>

:kword:`VARR`
  Specifies that the energy should not be considered invariant to rotations.
  Rotational variance is detected automatically, but sometimes it may be useful to enforce it.

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="VARR" APPEAR="Var Rot" KIND="SINGLE" LEVEL="ADVANCED">
              %%Keyword: VarR <advanced>
              <HELP>
              Specifies that the energy should not be considered invariant to rotations.
              Rotational variance is detected automatically, but sometimes it may be useful to enforce it.
              </HELP>
              </KEYWORD>

:kword:`NUMErical`
  Forces the calculation of numerical gradients even when analytic gradients are available.

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="NUMERICAL" APPEAR="Numerical gradients" KIND="SINGLE" LEVEL="ADVANCED">
              %%Keyword: Numerical <advanced>
              <HELP>
              Forces the calculation of numerical gradients even when analytic gradients are available.
              </HELP>
              </KEYWORD>

:kword:`SHAKe`
  Randomly modifies the initial coordinates of the atoms, maintaining the input (or computed)
  symmetry. This can be useful to avoid a geometry optimization converging to a higher-symmetry
  saddle point. The maximum displacement in the axes :math:`x`, :math:`y` and :math:`z` is read from the following
  real number. This number can be followed by :kword:`Bohr` or :kword:`Angstrom`, which indicates
  the unit in which the displacement is specified, the default is :kword:`Bohr`.

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="SHAKE" APPEAR="Shake" KIND="REAL" LEVEL="ADVANCED">
              %%Keyword: Shake <advanced>
              <HELP>
              Randomly modifies the initial coordinates of the atoms, maintaining the input (or computed)
              symmetry. This can be useful to avoid a geometry optimization converging to a higher-symmetry
              saddle point. The maximum displacement in the axes x, y and z is read from the following
              real number. This number can be followed by Bohr or Angstrom, which indicates
              the unit in which the displacement is specified, the default is Bohr.
              </HELP>
              </KEYWORD>

.. compound::

  Example: ::

    &GATEWAY
    COORD
    water.xyz
    BASIS
    STO-3G

  or, in short EMIL notation: ::

    &GATEWAY
    COORD=water.xyz; BASIS=STO-3G

Coordinate file may contain variables, as demonstrated in an example: ::

  >>FILE H2.input
  2
  scale $DD
  H 0.35 0 0
  H -0.35 0 0
  >>EOF

  >> FOREACH DD IN ( 0.9 1.0 1.1 )
  &GATEWAY
  COORD=$WorkDir/H2.input
  BASIS=STO-3G
  &SEWARD
  &SCF
  >>ENDDO

The atom name in XYZ file can contain an orbitrary label (to simplify assigning of different
basis sets). To indicate the label, use ``_``: e.g. ``C_SMALL``. The same label should be
defined in the basis section: ``BASIS=C_SMALL.ANO-S-MB``. The basis set label can be also
added into the name of an element: ::

  COORD
  1

  O.ANO-S-VDZP 0 0 0

XYZ file can also contain information about point charges. There are three possibilities to
setup atomic charges in XYZ file. The main option is to use ``Q`` as an element name, and in this
case the forth number, the charge, should be specified. Another possibility is to use element
names ended with minus sign. In this case, a formal charge for the element will be used.
E.g. ``H-``, ``Li-``, ``Na-``, ``K-`` defines :math:`+1` charge located in the corresponding location.
``Mg-``, ``Ca-`` --- defines charge :math:`+2`, ``Al-`` --- :math:`+3`, ``C-``, ``Si-`` :math:`+4`, for anions, ``F-``, ``Cl-``, ``Br-``, ``I-`` defines :math:`-1`,
``O-``, ``S-`` --- :math:`-2`. Finally, a label at the comment line of XYZ file --- CLUSTER followed by
an integer number can specify how many atoms are "real", so the rest will be treated as
charges with default values for this element.

Constraints
...........

In case of optimizations with constraints these are defined in the :program:`GATEWAY` input.
For a complete description of this keyword see :numref:`UG:sec:definition_of_internal_coordinates`.

.. class:: keywordlist

:kword:`CONStraints`
  This marks the start of the definition of the constraints
  which the optimization is subject to.
  This section is always ended by the keyword
  :kword:`End of Constraints`.
  This option can be used in conjunction with any definition of the
  internal coordinates.

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="CONSTRAINTS" APPEAR="Constraints" KIND="CUSTOM" LEVEL="BASIC">
              %%Keyword: Constraints <basic>
              <HELP>
              This marks the start of the definition of the constraints
              which the optimization is subject to.
              This section is always ended by the keyword "End of Constraints".
              Consult the manual for the details.
              </HELP>
              </KEYWORD>

:kword:`NGEXclude`
  This marks the start of the definition of additional restrictions for numerical differentiation.
  This section is always ended by the keyword :kword:`End of NGExclude`.
  The syntax of this section is like that of normal constraints, and the degrees of
  freedom specified here will be excluded from numerical differentiation (like phantom constraints).
  If a line containing only "Invert" is included inside the section,
  the definition is reversed and only these degrees of freedom are differentiated.
  :kword:`NGEXclude` is intended for use with the :kword:`KEEPOldGradient` keyword in :program:`ALASKA`,
  and can be combined with :kword:`CONStraints`, which will further reduce
  the numerical differentiation subspace :cite:`Stenrup2015`.
  Note that the value assigned to the constraints in this section is unused, but a ``Value`` block
  must still be included.

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="NGEXCLUDE" APPEAR="Exclude from numerical differentiation" KIND="CUSTOM" LEVEL="BASIC">
              %%Keyword: NGExclude <basic>
              <HELP>
              This marks the start of the definition of additional restrictions for numerical differentiation.
              This section is always ended by the keyword "End of NGExclude".
              </HELP>
              </KEYWORD>

Explicit auxiliary basis sets
.............................

The so-called Resolution of Identity (RI) technique (also called Density
Fitting, DF) is implemented in the |molcas| package. This option involves the use
of an auxiliary basis set in the effective computation of the 2-electron
integrals. |molcas| incorporates both the use of conventionally computed,
externally provided, auxiliary basis sets (RIJ, RIJK, and RIC types), and
on-the-fly generated auxiliary basis sets. The latter are atomic CD (aCD) or the
atomic compact CD (acCD) basis
sets, based on the Cholesky decomposition method. The externally provided
auxiliary basis sets are very compact, since they are tailored for special
wave function methods. However, they are not provided for all available valence
basis sets. The aCD or acCD RI auxiliary basis sets are a more general option and
provides auxiliary basis sets for any wave function model and valence basis set.

.. xmldoc:: <GROUP MODULE="GATEWAY" KIND="BOX" NAME="AUX" APPEAR="RI/DF options (optional)" LEVEL="BASIC">
            <HELP>
            Options of RI/DF definition of auxiliary basis sets.
            Set various thresholds and parameters for atomic CD auxiliary basis sets.
            </HELP>

.. class:: keywordlist

:kword:`RIJ`
  Use the RI-J basis in the density fitting (DF) approach to treat the two-electron integrals. Note that the valence
  basis set must have a supporting auxiliary basis set for this to work.

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="RIJ" APPEAR="RI-J option" KIND="SINGLE" EXCLUSIVE="RIJK,RIC,RICD,LOW,MEDI,HIGH" LEVEL="BASIC">
              %%Keyword: RIJ <basic>
              <HELP>
              Use the RI-J auxiliary basis in the density fitting (DF) approach to treat the two-electron integrals.
              Note that the valence basis set must have a supporting auxiliary basis set for this to work.
              </HELP>
              </KEYWORD>

:kword:`RIJK`
  Use the RI-JK auxiliary basis in the density fitting (DF) approach to treat the two-electron integrals. Note that the valence
  basis set must have a supporting auxiliary basis set for this to work.

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="RIJK" APPEAR="RI-JK option" KIND="SINGLE" EXCLUSIVE="RIJ,RIC,RICD,LOW,MEDI,HIGH" LEVEL="BASIC">
              %%Keyword: RIJK <basic>
              <HELP>
              Use the RI-JK auxiliary basis in the density fitting (DF) approach to treat the two-electron integrals.
              Note that the valence basis set must have a supporting auxiliary basis set for this to work.
              </HELP>
              </KEYWORD>

:kword:`RIC`
  Use the RI-C auxiliary basis in the density fitting (DF) approach to treat the two-electron integrals. Note that the valence
  basis set must have a supporting auxiliary basis set for this to work.

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="RIC" APPEAR="RI-C option" KIND="SINGLE" EXCLUSIVE="RIJ,RIJK,RICD,LOW,MEDI,HIGH" LEVEL="BASIC">
              %%Keyword: RIC <basic>
              <HELP>
              Use the RI-C auxiliary basis in the density fitting (DF) approach to treat the two-electron integrals.
              Note that the valence basis set must have a supporting auxiliary basis set for this to work.
              </HELP>
              </KEYWORD>

:kword:`RICD`
  Use the aCD or acCD approach :cite:`Aquilante:07b` to treat the two-electron integrals.
  This procedure will use an on-the-fly generated auxiliary basis set.

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="RICD" APPEAR="RI-aCD option" KIND="SINGLE" EXCLUSIVE="RIJ,RIJK,RIC,LOW,MEDI,HIGH" LEVEL="BASIC">
              %%Keyword: RICD <basic>
              <HELP>
              Use the aCD or acCD approach to treat the two-electron integrals.
              This procedure will use an on-the-fly generated auxiliary basis set.
              </HELP>
              </KEYWORD>

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="XRICD" KIND="SINGLE" LEVEL="UNDOCUMENTED" />

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="NOCD" KIND="SINGLE" LEVEL="UNDOCUMENTED" />

:kword:`CDTHreshold`
  Threshold for on-the-fly generation of aCD or acCD auxiliary basis sets for RI calculations
  (default value 1.0d-4).

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="CDTH" APPEAR="aCD threshold" KIND="REAL" DEFAULT_VALUE="1.0D-4" REQUIRE="RICD" LEVEL="ADVANCED">
              %%Keyword: CDThreshold <advanced>
              <HELP>
              Threshold for on-the-fly generation of aCD or acCD auxiliary basis sets for RI calculations
              (default value 1.0d-4).
              </HELP>
              </KEYWORD>

:kword:`SHAC`
  Skip high angular combinations à la Turbomole when creating on-the-fly basis sets
  (default off).

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="SHAC" APPEAR="Skip high angular combinations" KIND="SINGLE" REQUIRE="RICD" EXCLUSIVE="KHAC" LEVEL="ADVANCED">
              %%Keyword: SHAC <advanced>
              <HELP>
              Skip high angular combinations a la Turbomole when creating on-the-fly basis sets
              (default off).
              </HELP>
              </KEYWORD>

:kword:`KHAC`
  Keep high angular combinations when creating on-the-fly basis sets
  (default on).

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="KHAC" APPEAR="Keep high angular combinations" KIND="SINGLE" REQUIRE="RICD" EXCLUSIVE="SHAC" LEVEL="ADVANCED">
              %%Keyword: KHAC <basic>
              <HELP>
              Keep high angular combinations when creating on-the-fly basis sets
              (default on).
              </HELP>
              </KEYWORD>

:kword:`aCD basis`
  Generate an atomic CD (aCD) auxiliary basis sets (default off).

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="ACD" APPEAR="aCD auxiliary basis" KIND="SINGLE" REQUIRE="RICD" EXCLUSIVE="ACCD" LEVEL="ADVANCED">
              %%Keyword: aCD basis <basic>
              <HELP>
              Generate an atomic CD (aCD) auxiliary basis sets (default off).
              </HELP>
              </KEYWORD>

:kword:`acCD basis`
  Generate an atomic compact CD (acCD) auxiliary basis sets (default on).

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="ACCD" APPEAR="acCD auxiliary basis" KIND="SINGLE" REQUIRE="RICD" EXCLUSIVE="ACD" LEVEL="ADVANCED">
              %%Keyword: acCD basis <basic>
              <HELP>
              Generate an atomic compact CD (acCD) auxiliary basis sets (default on).
              </HELP>
              </KEYWORD>

.. xmldoc:: </GROUP>

.. index::
   single: Reaction field
   single: Cavity
   single: Solvent

.. _UG\:sec\:rfield:

Reaction field calculations
...........................

The effect of the solvent on the quantum chemical calculations has been
introduced in |molcas| through the reaction field created
by the surrounding environment, represented by a polarizable dielectric
continuum outside the boundaries of a cavity containing the solute molecule.
|molcas| supports Self Consistent Reaction Field (SCRF) and Multi
Configurational Self Consistent Reaction Field (MCSCRF) calculations within
the framework of the :program:`SCF` and the :program:`RASSCF` programs.
The reaction field, computed in a self-consistent fashion, can be
later added as a constant perturbation for the remaining programs, as
for example :program:`CASPT2`.

The purpose of this facility is to incorporate
the effect of the environment (a solvent or a solid matrix) on the studied molecule.
The utility itself it is not a program, but requires
an additional input which has to be provided to the
:program:`GATEWAY` program.
Two methods are available for SCRF calculations: one is based on the Kirkwood
model, the other is the so called Polarizable Continuum Model (PCM).
The reaction field is computed as the response of a dielectric medium polarized
by the solute molecule: the solute is placed in a "cavity" surrounded by the
dielectric. In Kirkwood model the cavity is always spherical, whereas in PCM the
cavity is modeled on the actual solute shape.

The possible set of parameters controlled by input are:

* the Kirkwood model,
* the PCM model, and
* one-electron integrals representing
  Pauli repulsion and external fields.

First a brief presentation of the Kirkwood and the PCM models.

.. index::
   single: Kirkwood model

The Kirkwood Model
::::::::::::::::::

The Kirkwood model is an expansion of the so-called Onsager model where
the surrounding will be characterized by its dielectric
permitivity and a radius describing a spherical cavity,
indicating where the dielectric medium starts.
(Note that all atoms in the studied molecule must be inside the spherical cavity.) The Pauli repulsion
due to the medium can be introduced by use of the spherical well
integrals which are generated by :program:`SEWARD`.
The charge distribution of the molecule will introduce
an electric field acting on the dielectric medium. This reaction field will interact with the
charge distribution of the molecule. This interaction will manifest itself as
a perturbation to the one-electron Hamiltonian. The perturbation will be
automatically computed in a direct fashion (no multipole integrals are stored on
disk) and added to the one-electron Hamiltonian. Due to the direct way in which
this contribution is computed rather high terms in the multipole expansion of the
charge can be afforded.

.. index::
   single: PCM

The Polarizable Continuum Model, PCM
::::::::::::::::::::::::::::::::::::

The PCM has been developed in order to describe the solvent reaction field in a
more realistic way, basically through the use of cavities of general shape, modeled on the
solute. The cavity is built as the envelope of spheres centered on solute atoms or atomic
groups (usually, hydrogen atoms are included in the same sphere as the heavy atoms they are
bonded to). The reaction field is described by means of apparent charges (solvation
charges) spread on the cavity surface, designed to reproduced the electrostatic potential
due to the polarized dielectric inside the cavity.
Such charges are used both to compute solute-solvent interactions (modifying the total energy
of the solute), and to perturb the molecular Hamiltonian through a suitable operator
(thus distorcing the solute wave-function, and affecting all the electronic properties).
The PCM operator contains both one- and two-electron terms: it is computed using
atomic integrals already present in the program, through a "geometry matrix"
connecting different points lying on the cavity surface. It can be shown that
with this approach the SCF and RASSCF variational procedures lead to the
free energy of the given molecule in solution: this is the thermodynamic meaning
of the SCF or CI energy provided by the program. More precisely, this is the
solute-solvent electrostatic contribution to the free energy
(of course, other terms depending on solute atomic motions, like vibrational and
rotational free energies, should be included separately);
it can be used to get a good
approximation of the solvation free energy, by subtracting the SCF or CI energy
computed in vacuo, and also to compute directly energy surfaces and reaction paths
in solution. On the other hand, the solute wave-function perturbed by the
reaction field can be used to compute any electronic property in solution.

Also other quantities can be computed, namely the cavitation free energy (due
the the work spent to create the cavity in the dielectric) and the
dispersion-repulsion free energy: these terms affect only the total free energy of the molecule,
and not its electronic distribution. They are collectively referred to as
non-electrostatic contributions.

Note that two other keywords are defined for the :program:`RASSCF`
program:
they refer to the CI root selected for the calculation of the reaction field (RFROOT), and
to the possibility to perform a non-equilibrium calculation (NONEQ) when vertical electronic
transitions are studied in solution. These keywords are referenced in the
:program:`RASSCF` section. To include the reaction field perturbation in a :program:`SCF`, :program:`RASSCF`, :program:`CASPT2` or :program:`RASSI`
calculation, another keyword must be specified (RFPERT), as explained in the
respective program sections.

Complete and detailed examples of how to add a reaction field,
through the Kirkwood or the PCM model, into quantum chemical
calculations in |molcas| is presented in :numref:`TUT:sec:cavity` of the
examples manual. The user is encouraged to read that section for further details.

Input for the Kirkwood and PCM models
:::::::::::::::::::::::::::::::::::::

.. _UG\:sec\:rfield_files:

Files
"""""

.. *********************************** this part should be revised

The reaction field calculations will store the information in the following files, which
will be used by the following programs

.. class:: filelist

:file:`ONEINT`
  One-electron integral file used to store the Pauli repulsion integrals

:file:`RUNFILE`
  Communications file. The last computed self-consistent reaction field (SCF or RASSCF)
  will be stored here to be used by following programs

:file:`GV.off`
  Input file for the external program "geomview" (see Tutorial section
  "Solvent models"), for the visualization of PCM cavities

Input
"""""
Below follows a description of the input to the reaction field utility in the
:program:`GATEWAY` program. The :program:`RASSCF` program has
its own keywords to compute reaction fields for excited states.

Compulsory keywords

:kword:`RF-Input`
  Activate reaction field options.

.. xmldoc:: <GROUP MODULE="GATEWAY" NAME="RF-INPUT" APPEAR="Reaction Field Options" KIND="BLOCK" LEVEL="ADVANCED">
            <HELP>
            Reaction field options.
            </HELP>

.. class:: keywordlist

:kword:`END Of RF-Input`
  This marks the end of the input to the reaction field utility.

  .. xmldoc:: %%Keyword: End of RF-input <compulsory>
              This marks the end of the input to the reaction field utility.

Optional keywords for the Kirkwood Model

.. class:: keywordlist

:kword:`REACtion Field`
  This command is exclusive to the Kirkwood model.
  It indicates the beginning of the specification of the
  reaction field parameters. The subsequent line will contain
  the dielectric constant of the medium, the radius of the
  cavity in Bohrs (the cavity is always centered around the
  origin), and the angular quantum number of the highest multipole
  moment used in the expansion of the change distribution of
  the molecule (only charge is specified as 0, charge and dipole
  moments as 1, etc.).
  The input specified below specifies that
  a dielectric permitivity of 80.0 is used, that the cavity radius is 14.00 a.u.,
  and that the expansion of the charge distribution is truncated after :math:`l=4`, i.e. hexadecapole
  moments are the last moments included in the expansion.
  Optionally a fourth argument can be added giving the value of the dielectric constant of the
  fast component of the solvent (default value 1.0).

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="REACTION" APPEAR="Onsager-Kirkwoord Model" KIND="CUSTOM" LEVEL="ADVANCED" EXCLUSIVE="PCM-MODEL">
              %%Keyword: Reaction field (RF) <basic>
              <HELP>
              This command is exclusive to the Kirkwood model.
              This indicated the beginning of the specification of the
              reaction field parameters. The subsequent line will contain
              the dielectric constant of the medium, the radius of the
              cavity in Bohrs (the cavity is always centered around the
              origin), and the angular quantum number of the highest multipole
              moment used in the expansion of the change distribution of
              the molecule (only charge is specified as 0, charge and dipole
              moments as 1, etc.).
              The input specified below specifies that
              a dielectric permitivity of 80.0 is used, that the cavity radius is 14.00 a.u.,
              and that the expansion of the charge distribution is truncated after l=4, i.e. hexadecapole
              moments are the last moments included in the expansion.
              Optionally a fourth argument can be added giving the value of the dielectric constant of the
              fast component of the solvent (default value 1.0).
              </HELP>
              </KEYWORD>

Sample input for the reaction field part (Kirkwood model) ::

  RF-Input
  Reaction field
  80.0 14.00 4
  End Of RF-Input

Sample input for a complete reaction field calculation using the Kirkwood model.
The :program:`SCF` computes the reaction field in a
self consistent manner while the :program:`MRCI`
program adds the effect as a constant perturbation.

.. extractfile:: ug/RF.input

  &GATEWAY
  Title = HF molecule
  Symmetry
  X Y
  Basis set
  F.ANO-S...3S2P.
  F      0.00000   0.00000   1.73300
  End of basis
  Basis set
  H.ANO-S...2S.
  H      0.00000   0.00000   0.00000
  End of basis
  Well integrals
   4
   1.0 5.0  6.75
   1.0 3.5  7.75
   1.0 2.0  9.75
   1.0 1.4 11.75

  RF-Input
  Reaction field
   80.0 4.75 4
  End of RF-Input

  &SEWARD

  &SCF
  Occupied =  3 1 1 0

  &MOTRA
  LumOrb
  Frozen   =  1 0 0 0
  RFPert

  &GUGA
  Electrons =     8
  Spin      =     1
  Inactive  =     2    1    1    0
  Active    =     0    0    0    0
  CiAll     =     1

  &MRCI
  SDCI

Optional keywords for the PCM Model

.. xmldoc:: <GROUP MODULE="GATEWAY" NAME="PCM" APPEAR="PCM Options" KIND="BOX" LEVEL="ADVANCED">

.. class:: keywordlist

:kword:`PCM-model`
  If no other keywords are specified, the program will execute a standard PCM calculation
  with water as solvent. The solvent reaction field will be included in all the
  programs (:program:`SCF`, :program:`RASSCF`, :program:`CASPT2`, etc.)
  invoked after :program:`SEWARD`: note that in some cases additional keywords are required
  in the corresponding program sections. Some PCM parameters can be changed through the following
  keywords.

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="PCM-MODEL" APPEAR="PCM Model" KIND="SINGLE" LEVEL="ADVANCED" EXCLUSIVE="REACTION">
              %%Keyword: PCM-model (RF) <basic>
              <HELP>
              If no other keywords are specified, the program will execute a standard PCM calculation
              with water as solvent. The solvent reaction field will be included in all the
              programs (SCF, RASSCF, CASPT2, etc.)
              invoked after SEWARD: note that in some cases additional keywords are required
              in the corresponding program sections. Many PCM parameters can be changed through the following
              keywords.
              </HELP>
              </KEYWORD>

:kword:`SOLVent`
  Used to indicate which solvent is to be simulated. The name of the requested solvent
  must be written in the line below this keyword. Find implemented solvents in the PCM model below this section.

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="SOLVENT" APPEAR="Solvent" KIND="CHOICE" LEVEL="ADVANCED" REQUIRE="PCM-MODEL" LIST="WATER,ACETONITRILE,METHANOL,ETHANOL,ISOQUINOLINE,QUINOLINE,CHLOROFORM,ETHYLETHER,METHYLENECHLORIDE,DICHLOROETHANE,CARBONTETRACHLORIDE,BENZENE,TOLUENE,CHLOROBENZENE,NITROMETHANE,HEPTANE,CYCLOHEXANE,ANILINE,ACETONE,TETRAHYDROFURAN,DIMETHYLSULFOXIDE,ARGON,KRYPTON,XENON">
              %%Keyword: Solvent (RF) <basic>
              <HELP>
              Used to indicate which solvent is to be simulated.
              </HELP>
              The name of the requested solvent
              must be written in the line below this keyword. Allowed solvents are:
              WATER, ACETONITRILE, METHANOL, ETHANOL, ISOQUINOLINE,
              QUINOLINE, CHLOROFORM, ETHYLETHER, METHYLENECHLORIDE,
              DICHLOROETHANE, CARBONTETRACHLORIDE, BENZENE, TOLUENE,
              CHLOROBENZENE, NITROMETHANE, HEPTANE, CYCLOHEXANE, ANILINE,
              ACETONE, TETRAHYDROFURAN, DIMETHYLSULFOXIDE, ARGON, KRYPTON,
              XENON
              </KEYWORD>

:kword:`DIELectric constant`
  Defines a different dielectric constant for the selected solvent; useful to describe
  the system at temperatures other that 298 K, or to mimic solvent mixtures.
  The value is read in the line below the keyword.
  An optional second value might be added on the same line which
  defines a different value for the infinite frequency dielectric constant for
  the selected solvent (this is used in non-equilibrium calculations; by
  default it is defined for each solvent at 298 K).

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="DIELECTRIC" APPEAR="Dielectric constant" KIND="REALS" SIZE="2" LEVEL="ADVANCED" REQUIRE="PCM-MODEL">
              %%Keyword: Dielectric constant (RF) <basic>
              <HELP>
              Defines a different dielectric constant for the selected solvent; useful to describe
              the system at temperatures other that 298 K, or to mimic solvent mixtures.
              The value is read in the line below the keyword.
              An optional second value might be added on the same line which
              defines a different value for the infinite frequency dielectric constant for
              the selected solvent (this is used in non-equilibrium calculations; by
              default it is defined for each solvent at 298 K).
              </HELP>
              </KEYWORD>

  .. :kword:`INFInite frequency dielectric constant`
       Defines a different value for the infinite frequency dielectric constant for
       the selected solvent (this is used in non-equilibrium calculations; by
       default it is defined for each solvent at 298 K).
       The value is read in the line below the keyword.

       .. .. xmldoc:: %%Keyword: Infinite frequency dielectric constant (RF) <advanced>
                      Defines a different value for the infinite frequency dielectric constant for
                      the selected solvent (this is used in non-equilibrium calculations; by
                      default it is defined for each solvent at 298 K).
                      The value is read in the line below the keyword.

:kword:`CONDuctor version`
  It requires a PCM calculation where the solvent is represented as a polarized conductor:
  this is an approximation to the dielectric model which works very well for
  polar solvents (i.e. dielectric constant greater than about 5), and it has some
  computational advantages being based on simpler equations. It can be useful in cases
  when the dielectric model shows some convergence problems.

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="CONDUCTOR" APPEAR="Conductor Model" KIND="SINGLE" LEVEL="ADVANCED" REQUIRE="PCM-MODEL">
              %%Keyword: Conductor version (RF) <advanced>
              <HELP>
              It requires a PCM calculation where the solvent is represented as a polarized conductor:
              this is an approximation to the dielectric model which works very well for
              polar solvents (i.e. dielectric constant greater than about 5), and it has some
              computational advantages being based on simpler equations. It can be useful in cases
              when the dielectric model shows some convergence problems.
              </HELP>
              </KEYWORD>

:kword:`AAREa`
  It is used to define the average area (in Å\ |2|)
  of the small elements on the cavity surface
  where solvation charges are placed; when larger elements are chosen, less charges
  are defined, what speeds up the calculation but risks to worsen the results. The
  default value is 0.4 Å\ |2| (i. e. 60 charges on a sphere of radius 2 Å).
  The value is read in the line below the keyword.

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="AARE" APPEAR="Tessera Average Area" KIND="REAL" LEVEL="ADVANCED" DEFAULT_VALUE="0.4" REQUIRE="PCM-MODEL">
              %%Keyword: AAREa (RF) <advanced>
              <HELP>
              It is used to define the average area (in A^2) of the small elements on the cavity surface
              where solvation charges are placed; when larger elements are chosen, less charges
              are defined, what speeds up the calculation but risks to worsen the results. The
              default value is 0.4 A^2 (i.e. 60 charges on a sphere of radius 2 A).
              The value is read in the line below the keyword.
              </HELP>
              </KEYWORD>

:kword:`R-MIn`
  It sets the minimum radius (in Å) of the spheres that the program adds to the atomic
  spheres in order to smooth the cavity surface (default 0.2 Å).
  For large solute, if the programs
  complains that too many sphere are being created, or if computational times
  become too high, it can be useful to enlarge this value (for example to 1 or 1.5
  Å), thus reducing the number of added spheres.
  The value is read in the line below the keyword.

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="R-MIN" APPEAR="Minimum sphere radius" KIND="REAL" LEVEL="ADVANCED" DEFAULT_VALUE="2.0" REQUIRE="PCM-MODEL">
              %%Keyword: R-min (RF) <advanced>
              <HELP>
              It sets the minimum radius (in A) of the spheres that the program adds to the atomic
              spheres in order to smooth the cavity surface (default 0.2 A).
              For large solute, if the programs
              complains that too many sphere are being created, or if computational times
              become too high, it can be useful to enlarge this value (for example to 1 or 1.5
              A), thus reducing the number of added spheres.
              The value is read in the line below the keyword.
              </HELP>
              </KEYWORD>

:kword:`PAULing`
  It invokes the use of Pauling's radii to build the solute cavity: in
  this case, hydrogens get their own sphere (radius 1.2 Å).

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="PAULING" APPEAR="Pauling radii" KIND="SINGLE" LEVEL="ADVANCED" DEFAULT_VALUE="2.0" REQUIRE="PCM-MODEL">
              %%Keyword: Pauling (RF) <advanced>
              <HELP>
              It invokes the use of Pauling's radii to build the solute cavity: in
              this case, hydrogens get their own sphere (radius 1.2 A).
              </HELP>
              </KEYWORD>

:kword:`SPHEre radius`
  It is used to provide sphere radii from input: for each sphere given
  explicitly by the user, the keyword "Sphere radius" is required,
  followed by a line containing two numbers: an integer indicating the
  atom where the sphere has to be centered, and a real indicating its
  radius (in Å). For example, "Sphere radius" followed by "3 1.5"
  indicates that a sphere of radius 1.5 Å is placed around atom #3;
  "Sphere radius" followed by "4 2.0" indicates that another sphere of
  radius 2 Å is placed around atom #4 and so on.

  .. xmldoc:: %%Keyword: Sphere radius (RF) <advanced>
              It is used to provide sphere radii from input: for each sphere given
              explicitly by the user, the keyword 'Sphere radius' is required,
              followed by a line containing two numbers: an integer indicating the
              atom where the sphere has to be centered, and a real indicating its
              radius (in A). For example, 'Sphere radius' followed by '3 1.5'
              indicates that a sphere of radius 1.5 A is placed around atom 3;
              'Sphere radius' followed by '4 2.0' indicates that another sphere of
              radius 2 A is placed around atom 4 and so on.

.. xmldoc:: </GROUP>

.. xmldoc:: </GROUP>

Solvents implemented in the PCM model are

.. %---- Table of allowed solvents ------

.. _tab\:pcm_solvents:

=================== ==========
Name                Dielectric
                    constant
=================== ==========
water                    78.39
dimethylsulfoxide        46.70
nitromethane             38.20
acetonitrile             36.64
methanol                 32.63
ethanol                  24.55
acetone                  20.70
isoquinoline             10.43
dichloroethane           10.36
quinoline                 9.03
methylenchloride          8.93
tetrahydrofuran           7.58
aniline                   6.89
chlorobenzene             5.62
chloroform                4.90
ethylether                4.34
toluene                   2.38
benzene                   2.25
carbontetrachloride       2.23
cyclohexane               2.02
heptane                   1.92
xenon                     1.71
krypton                   1.52
argon                     1.43
=================== ==========

Sample input for the reaction field part (PCM model): the solvent is
water, a surface element average area of 0.2 Å\ |2| is requested. ::

  RF-input
  PCM-model
  Solvent
  water
  AAre
  0.2
  End of RF-input

Sample input for a standard PCM calculation in water.
The :program:`SCF` and :program:`RASSCF` programs compute the reaction field
self consistently and add its contribution to the Hamiltonian. The :program:`RASSCF` is
repeated twice: first the ground state is determined, then a non-equilibrium
calculation on the first excited state is performed.

.. extractfile:: ug/RF.formaldehyde.input

  &GATEWAY
  Coord
  4
  formaldehyde
  O 0.000000 0.000000 -1.241209
  C 0.000000 0.000000 0.000000
  H 0.000000 0.949585 0.584974
  H 0.000000 -0.949585 0.584974

  Basis = STO-3G
  Group = C1
  RF-input
  PCM-model
  solvent = water
  End of RF-input

  &SEWARD ; &SCF

  &RASSCF
  nActEl   = 4 0 0
  Symmetry = 1
  Inactive = 6
  Ras2     = 3
  CiRoot
  1 1
  1
  LumOrb

  &RASSCF
  nActEl   = 4 0 0
  Symmetry = 1
  Inactive = 6
  Ras2     = 3
  CiRoot
  2 2
  1 2
  0 1
  JOBIPH
  NonEq
  RFRoot  = 2

Again the user is recommended to read :numref:`TUT:sec:cavity` of the
examples manual for further details.

Keywords associated to one-electron integrals
.............................................

.. xmldoc:: <GROUP MODULE="GATEWAY" NAME="ONE-ELECTRON" APPEAR="1-electron integral options" KIND="BOX" LEVEL="ADVANCED">

.. class:: keywordlist

:kword:`FNMC`
  Request that the so-called Finite Nuclear Mass Correction, excluded by the Born--Oppenheimer approximation,
  be added to the one-electron Hamiltonian.

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="FNMC" APPEAR="Finite nuclear mass correction" KIND="SINGLE" LEVEL="ADVANCED">
              %%Keyword: FNMC <advanced>
              <HELP>
              Request that the so-called Finite Nuclear Mass Correction, excluded by the Born-Oppenheimer approximation,
              be added to the one-electron Hamiltonian.
              </HELP>
              </KEYWORD>

:kword:`WELL integrals`
  Request computation of Pauli repulsion integrals for dielectric
  cavity reaction field calculations.
  The first line specifies the total number of primitive well integrals in the
  repulsion integral. Then follows a number of lines, one for each
  well integral, specifying the coefficient of the well integral in the
  linear combination of the well integrals which defines the repulsion integral,
  the exponent of the well integral, and the distance of the center of the
  Gaussian from the origin. In total three entries on each line.
  All entries in atomic units.
  If zero or a negative number is specified for the number of well integrals
  a standard set of 3 integrals with their position adjusted for the radius of
  the cavity will be used.
  If the distance of the center of the Gaussian from the origin is
  negative displacements relative to the cavity radius is assumed.

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="WELL" APPEAR="Well integrals" KIND="REALS_COMPUTED" SIZE="3" LEVEL="BASIC">
              <ALTERNATE KIND="INT" />
              %%Keyword: Well integrals <basic>
              <HELP>
              Request computation of Pauli repulsion integrals for dielectric
              cavity reaction field calculations.
              </HELP>
              The first line specifies the total number of primitive well integrals in the
              repulsion integral. Then follows a number of lines, one for each
              well integral, specifying the coefficient of the well integral in the
              linear combination of the well integrals which defines the repulsion integral,
              the exponent of the well integral, and the distance of the center of the
              Gaussian from the origin. In total three entries on each line.
              All entries in atomic units.
              If zero or a negative number is specified for the number of well integrals
              a standard set of 3 integrals with their position adjusted for the radius of
              the cavity will be used.
              If the distance of the center of the Gaussian from the origin is
              negative displacements relative to the cavity radius is assumed.
              </KEYWORD>

:kword:`XFIEld integrals`
  Request the presence of an external electric field represented by a
  number of partial charges and dipoles. Optionally, polarisabilities may be specified whose
  induced dipoles are determined self-consistently during the SCF iteration.
  The first line may contain, apart from the first integer [nXF] (number of centers), up to
  four additional integers. The second integer [nOrd] specifies the maximum multipole order,
  or -1 signifying no permanent multipoles. Default is 1 (charges and dipoles). The third
  integer [p] specifies the type of external polarisabilities: 0 (default) no polarisabilities,
  1 (isotropic), or 2 (anisotropic). The fourth integer [nFrag] specifies the number of fragments one
  multipole may contribute to (relevant only if polarisabilities are present). The default is 0,
  meaning that each permanent multipole is only excluded in the calculation of the field at its own
  polarisability, 1 means that one gives a fragment number to each multipole and that the static
  multipoles do not contribute to the polarising field within the same fragment, whereas 2 can be
  used in more complex situations, e.g. polymers, allowing you to specify a second fragment number
  so that junction atoms does not contribute to either of the neighbouring fragments.
  Finally, the fifth and last integer [nRead] (relevant only if Langevin dipoles are used) may
  be 0 or 1 (where 0 is default), specifying whether an element number (e.g. 8 for oxygen) should be
  read for each multipole. In that case the default radius for that element is used to determine which
  Langevin grid points should be annihilated. A negative element number signifies that a particular
  radius should be used for that multipole, in thousands of a Bohr (-1400 meaning 1.4 Bohr).
  Then follows nXF lines, one for each center. On each line is first nFrag+nRead (which may equal 0)
  integers, specifying the fragments that the multipole should not contribute to (the first fragment is
  taken as the fragment that the polarisability belongs to) and the element number. Then follows
  the three coordinates of the center, followed by the multipoles and polarisabilities. The number of
  multipole entries is 0 for nOrd=-1, 1 for nOrd=0, 4 for nOrd=1, and 10 for nOrd=2. The number of
  polarisability entries are 0 for p=0, 1 for p=1, and 6 for p=2. The order of quadrupole moment and
  anisotropic polarisability entries is xx, xy, xz, yy, yz, zz. If default is used, i.e. only specifying
  the number of centers on the first line, each of these lines will contain 7 entries (coordinates,
  charge, and dipole vector). All entries are in atomic units, if not otherwise requested by the :kword:`Angstrom`
  keyword that must be placed between nXF and nOrd. All these data can be stored in a separate file whose
  name must be passed as an argument of the :kword:`XField` keyword.

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="XFIELD" APPEAR="External field" KIND="CUSTOM" LEVEL="BASIC">
              <ALTERNATE KIND="STRING" />
              %%Keyword: Xfield integrals <basic>
              <HELP>
              Request the presence of an external electric field represented by a
              number of partial charges and dipoles. Optionally, polarisabilities may be specified whose
              induced dipoles are determined self-consistently during the SCF iteration.
              </HELP>
              The first line may contain, apart from the first integer (nXF) (number of centers), up to
              four additional integers. The second integer (nOrd) specifies the maximum multipole order,
              or -1 signifying no permanent multipoles. Default is 1 (charges and dipoles). The third
              integer (p) specifies the type of external polarisabilities: 0 (default) no polarisabilities,
              1 (isotropic), or 2 (anisotropic). The fourth integer (nFrag) specifies the number of fragments one
              multipole may contribute to (relevant only if polarisabilities are present). The default is 0,
              meaning that each permanent multipole is only excluded in the calculation of the field at its own
              polarisability, 1 means that one gives a fragment number to each multipole and that the static
              multipoles do not contribute to the polarising field within the same fragment, whereas 2 can be
              used in more complex situations, e.g. polymers, allowing you to specify a second fragment number
              so that junction atoms does not contribute to either of the neighbouring fragments.
              Finally, the fifth and last integer (nRead) (relevant only if Langevin dipoles are used) may
              be 0 or 1 (where 0 is default), specifying whether an element number (e.g. 8 for oxygen) should be
              read for each multipole. In that case the default radius for that element is used to determine which
              Langevin grid points should be annihilated. A negative element number signifies that a particular
              radius should be used for that multipole, in thousands of a Bohr (-1400 meaning 1.4 Bohr).
              Then follows nXF lines, one for each center. On each line is first nFrag+nRead (which may equal 0)
              integers, specifying the fragments that the multipole should not contribute to (the first fragment is
              taken as the fragment that the polarisability belongs to) and the element number. Then follows
              the three coordinates of the center, followed by the multipoles and polarisabilities. The number of
              multipole entries is 0 for nOrd=-1, 1 for nOrd=0, 4 for nOrd=1, and 10 for nOrd=2. The number of
              polarisability entries are 0 for p=0, 1 for p=1, and 6 for p=2. The order of quadrupole moment and
              anisotropic polarisability entries is xx, xy, xz, yy, yz, zz. If default is used, i.e. only specifying
              the number of centers on the first line, each of these lines will contain 7 entries (coordinates,
              charge, and dipole vector). All entries are in atomic units, if not otherwise requested by the Angstrom
              keyword that must be placed between nXF and nOrd. All these data can be stored in a separate file whose
              name must be passed as an argument of the XField keyword.
              </KEYWORD>

:kword:`SDIPole`
  Requests computation of velocity integrals. This is usually enabled by default.

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="SDIPOLE" APPEAR="Velocity integrals" KIND="SINGLE" LEVEL="ADVANCED">
              %%Keyword: Sdipole <basic>
              <HELP>
              Requests computation of velocity integrals.
              </HELP>
              </KEYWORD>

:kword:`ANGM`
  Supplement
  :file:`ONEINT` for transition angular momentum calculations.
  Entry which specifies the angular momentum origin (in au).
  By default this is enabled with the origin at the center of mass.

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="ANGM" APPEAR="Angular momentum" KIND="REALS" SIZE="3" LEVEL="ADVANCED">
              %%Keyword: Angm <basic>
              <HELP>
              Supplement the file for transition angular momentum calculations.
              Enter the angular momentum operator origin (in au).
              </HELP>
              The keyword is followed by a card which specifies the angular momentum
              origin (in au).
              </KEYWORD>

:kword:`OMQI`
  Supplement
  :file:`ONEINT` for transition orbital magnetic quadrupole calculations.
  Entry which specifies the orbital magnetic quadrupole origin (in au).

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="OMQI" APPEAR="Orbital magnetic quadrupole" KIND="REALS" SIZE="3" LEVEL="ADVANCED">
              %%Keyword: OMQI <basic>
              <HELP>
              Supplement the file for transition orbital magnetic quadrupole calculations.
              Enter the orbital magnetic quadrupole operator origin (in au).
              </HELP>
              The keyword is followed by a card which specifies the orbital magnetic quadrupole
              origin (in au).
              </KEYWORD>

:kword:`AMPR`
  Request the computation of angular momentum product integrals.
  The keyword is followed by values which specifies the angular momentum
  origin (in au).

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="AMPR" APPEAR="Angular momentum product" KIND="REALS" SIZE="3" LEVEL="ADVANCED">
              <HELP>
              Request the computation of angular momentum product integrals and specify the
              angular momentum origin (in au).
              </HELP>
              %%Keyword: Ampr <basic>
              Request the computation of angular momentum product integrals.
              The keyword is followed by a card which specifies the angular momentum
              origin (in au).
              </KEYWORD>

:kword:`DSHD`
  Requests the computation of diamagnetic shielding integrals. The first
  entry specifies the gauge origin. Then follows an integer
  specifying the number of points at which the diamagnetic
  shielding will be computed. If this entry is zero, the diamagnetic
  shielding will be computed at each nucleus. If nonzero, then the
  coordinates (in au) for each origin has to be supplied, one entry for each
  origin.

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="DSHD" APPEAR="Diamagnetic shielding" KIND="STRINGS" SIZE="10" LEVEL="ADVANCED">
              %%Keyword: DSHD <basic>
              <HELP>
              Activate the computation of diamagnetic shielding integrals. The first entry
              specifies the gauge origin. On the subsequent entries an
              integer specifying the number of points at which the diamagnetic
              shielding will be computed. If this entry is zero, the diamagnetic
              shielding will be computed at each nucleus. If nonzero, then the
              coordinates (in au) for each origin has to be supplied, one entry for each
              origin.
              </HELP>
              </KEYWORD>

  .. :kword:`DOUGlas-kroll`
     Explicit request that the one-electron Hamiltonian include the scalar relativistic
     effects according to the so-called Douglas--Kroll transformation.

  ..   .. xmldoc:: %%Keyword: Douglas-Kroll <basic>
                 Explicit request that the one-electron Hamiltonian include the scalar relativistic
                 effects according to the so-called Douglas-Kroll transformation.
                 This option is automatically invoked for the ANO-RCC and ANO-DK3 basis sets.

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="DOUGLAS-KROLL" KIND="SINGLE" LEVEL="UNDOCUMENTED" />

:kword:`RX2C`
  Request the scalar relativistic X2C (eXact-two-Component) corrections to the
  one-electron Hamiltonian as well as the property integrals.

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="RX2C" APPEAR="Relativistic X2C integrals" KIND="SINGLE" EXCLUSIVE="RBSS" LEVEL="BASIC">
              %%Keyword: RX2C <basic>
              <HELP>
              Request the scalar relativistic X2C (eXact-two-Component) corrections to the
              one-electron Hamiltonian as well as the property integrals.
              </HELP>
              </KEYWORD>

:kword:`RBSS`
  Request the scalar relativistic BSS (Barysz--Sadlej--Snijders) corrections to the
  one-electron Hamiltonian as well as the property integrals. The non-iterative
  scheme is employed for the construction of BSS transformation.

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="RBSS" APPEAR="Relativistic BSS integrals" KIND="SINGLE" EXCLUSIVE="RX2C" LEVEL="BASIC">
              %%Keyword: RBSS <basic>
              <HELP>
              Request the scalar relativistic BSS (Barysz-Sadlej-Snijders) corrections to the
              one-electron Hamiltonian as well as the property integrals. The non-iterative
              scheme is employed for the construction of BSS transformation.
              </HELP>
              </KEYWORD>

:kword:`NOAMfi`
  Explicit request for no computation of atomic mean-field integrals.

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="NOAM" APPEAR="No AMFI integrals" KIND="SINGLE" EXCLUSIVE="AMFI" LEVEL="BASIC">
              %%Keyword: NOAMFI <basic>
              <HELP>
              Explicit request for no computation of atomic mean-field integrals.
              </HELP>
              </KEYWORD>

:kword:`AMFI`
  Explicit request for the computation of atomic mean-field integrals (used in
  subsequent spin--orbit calculations). These integrals are computed by default for the
  ANO-RCC and ANO-DK3 basis sets.

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="AMFI" APPEAR="AMFI integrals option" KIND="SINGLE" EXCLUSIVE="NOAM" LEVEL="BASIC">
              %%Keyword: AMFI <basic>
              <HELP>
              Explicit request for the computation of atomic mean-field integrals (used in
              subsequent spin-orbit calculations). These integrals are computed by default for
              relativistic basis sets like the ANO-RCC and ANO-DK3 basis sets.
              </HELP>
              </KEYWORD>

:kword:`EPOT`
  An integer follows which represents the
  number of points for which the electric potential will be computed. If
  this number is zero, the electric potential acting on each nucleus will be
  computed. If nonzero, then the coordinates (in au) for each point have to be
  supplied, one entry for each point.
  This keyword is mutually exclusive with :kword:`EFLD` and :kword:`FLDG`.

  .. xmldoc:: <SELECT MODULE="GATEWAY" NAME="EF" APPEAR="Electric potential, field and field gradient options" LEVEL="BASIC" CONTAINS="EPOT,EFLD,FLDG">

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="EPOT" APPEAR="Electric potential" KIND="CUSTOM" LEVEL="ADVANCED">
              <HELP>
              Activate the computation of the electric potential at some points.
              The first entry is the number of points at which this should be computed.
              The coordinates (in au) for each point have to be
              supplied on the subsequent entries.
              If the number of points is zero, the electric potential on each nucleus will be computed.
              </HELP>
              %%Keyword: EPOT <basic>
              An integer follows which represents the
              number of points for which the electric potential will be computed. If
              this number is zero, the electric potential acting on each nucleus will be
              computed. If nonzero, then the coordinates (in au) for each point have to be
              supplied, one entry for each point.
              This keyword is mutually exclusive with EFLD and FLDG.
              </KEYWORD>

:kword:`EFLD`
  An integer follows which represents the
  number of points for which the electric potential and electric field will be computed. If
  this number is zero, the electric field acting on each nucleus will be
  computed. If nonzero, then the coordinates (in au) for each point have to be
  supplied, one entry for each point.
  This keyword is mutually exclusive with :kword:`EPOT` and :kword:`FLDG`.

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="EFLD" APPEAR="Electric field" KIND="CUSTOM" LEVEL="ADVANCED">
              <HELP>
              Activate the computation of the electric potential and field at some points.
              The first entry is the number of points at which this should be computed.
              The coordinates (in au) for each point have to be
              supplied on the subsequent entries.
              If the number of points is zero, the electric field on each nucleus will be computed.
              </HELP>
              %%Keyword: EFLD <basic>
              Followed by a card with an integer entry which represents the
              number of points for which the electric potential and electric field will be computed. If
              this number is zero, the electric field acting on each nucleus will be
              computed. If nonzero, then the coordinates (in au) for each point have to be
              supplied, one entry for each point.
              This keyword is mutually exclusive with EPOT and FLDG.
              </KEYWORD>

:kword:`FLDG`
  An integer required which represents the
  number of points for which the electric potential, electric field and electric field gradient will be
  computed. If this number is zero, the electric field gradient acting
  on each nucleus will be computed. If nonzero, then either the coordinates (in au) for
  each point or labels for each atom center have to be supplied, one entry for each point.
  In case a label is supplied it must match one of those given previous in the input during specification
  of the coordinates of the atom centers. Using a label instead of a coordinate can e.g. be useful
  in something like a geometry optimization where the coordinate isn't known when the input is written.
  This keyword is mutually exclusive with :kword:`EPOT` and :kword:`EFLD`.

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="FLDG" APPEAR="Electric field gradient" KIND="CUSTOM" LEVEL="ADVANCED">
              <HELP>
              Activate the computation of the electric potential, field and field gradient at some points.
              The first entry is the number of points at which this should be computed.
              The coordinates (in au) for each point have to be
              supplied on the subsequent entries.
              If the number of points is zero, the electric field gradient on each nucleus will be computed.
              </HELP>
              %%Keyword: FLDG <basic>
              An integer required which represents the
              number of points for which the electric potential, electric field and electric field gradient will be
              computed. If this number is zero, the electric field gradient acting
              on each nucleus will be computed. If nonzero, then either the coordinates (in au) for
              each point or labels for each atom center have to be supplied, one entry for each point.
              In case a label is supplied it must match one of those given previous in the input during specification
              of the coordinates of the atom centers. Using a label instead of a coordinate can e.g. be useful
              in something like a geometry optimization where the coordinate isn't known when the input is written.
              This keyword is mutually exclusive with EPOT and EFLD.
              </KEYWORD>

  .. xmldoc:: </SELECT>

:kword:`EMPC`
  Use point charges specified by the keyword :kword:`XField` when calculating the Orbital-Free Embedding potential.

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="EMPC" APPEAR="Embedded Point Charges" KIND="SINGLE" LEVEL="BASIC">
              %%Keyword: EMPC <basic>
              <HELP>
              Use point charges specified by the keyword XFIELD when calculating the Orbital-Free Embedding potential.
              </HELP>
              </KEYWORD>

:kword:`RF-Input`
  Specification of reaction field parameters, consult the reaction field section of this manual.

  .. xmldoc:: %%Keyword: RF-input <basic>
              Specification of reaction field parameters, consult the reaction field section of this
              manual.

.. xmldoc:: </GROUP>

Keywords associated with nuclear charge distribution models
...........................................................

Input parameters associated with different models of the nuclear charge distribution. The
default is to use a point charge representation.

.. xmldoc:: <GROUP MODULE="GATEWAY" NAME="NUCLEAR" APPEAR="Nuclear Models" KIND="BOX" LEVEL="ADVANCED">

.. class:: keywordlist

:kword:`FINIte`
  Request a finite center representation of the nuclei by a single exponent s-type Gaussian.

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="FINITE" APPEAR="Activate Gaussian Nuclear Charge Distribution" KIND="SINGLE" EXCLUSIVE="MGAUSSIAN" LEVEL="ADVANCED">
              %%Keyword: Finite <basic>
              <HELP>
              Request a finite center representation of the nuclei by a single exponent s-type Gaussian.
              </HELP>
              </KEYWORD>

:kword:`MGAUSsian`
  Request a finite center representation of the nuclei by a modified Gaussian.

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="MGAUSSIAN" APPEAR="Activate Modified Gaussian Charge Distribution" KIND="SINGLE" EXCLUSIVE="FINITE" LEVEL="ADVANCED">
              %%Keyword: MGauss <basic>
              <HELP>
              Request a finite center representation of the nuclei by a modified Gaussian.
              </HELP>
              </KEYWORD>

.. xmldoc:: </GROUP>

The Saddle method for transition state optimization
...................................................

The Saddle method :cite:`Saddle_method` is a method to locate transition states (TS). The method, in practice, can be viewed as a
series of constrained optimization along the reaction path, which connects two starting structure (could be
the reactants and products of a reaction), to locate the region of the TS and a subsequent unconstrained optimization
to locate the TS. The only data needed for the procedure are the energies and coordinates of the two structures.
**Note** that this option will overwrite the
coordinates which have already been specified with the normal input of the molecular geometry. However, this does
not make that input section redundant and should always be included.

.. xmldoc:: <GROUP MODULE="GATEWAY" NAME="SADDLEMETHOD" APPEAR="Saddle Method" KIND="BOX" LEVEL="ADVANCED">

.. class:: keywordlist

:kword:`RP-Coordinates`
  This activates the Saddle method for TS geometry optimization.
  The line is followed by an integer specifying the number of symmetry unique coordinates to be specified. This
  is followed by two sets of input --- one line with the energy and then the Cartesian coordinates in bohr --- for
  each of the two starting structures of the Saddle method. Note that the order of the coordinates must always
  match the order specified with the conventional input of the coordinates of the molecular system.
  Alternatively, two lines with the filenames containing the coordinates of reactants and products, respectively,
  (in XYZ format) can be given.

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="RP-COORD" APPEAR="Reactants and Products coordinates" KIND="STRINGS" SIZE="2" LEVEL="ADVANCED">
              %%Keyword: RP-Coordinates <advanced>
              <HELP>
              This activates the Saddle method for TS geometry optimization.
              The line is followed by an integer specifying the number of symmetry unique coordinates to be specified. This
              is followed by two sets of input - one line with the energy and then the Cartesian coordinates in bohr - for
              each of the two starting structures of the Saddle method. Note that the order of the coordinates must always
              match the order specified with the conventional input of the coordinates of the molecular system.
              Alternatively, two lines with the filenames containing the coordinates of reactants and products, respectively,
              (in XYZ format) can be given.
              </HELP>
              </KEYWORD>

:kword:`NOALign`
  By default, the two starting structures are aligned to minimize the root mean square distance (RMSD) between them,
  in particular, the first structure is moved and the second structure remains fixed.
  If this keyword is given, the starting structures are used as given.

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="NOALIGN" APPEAR="No align" KIND="SINGLE" LEVEL="ADVANCED" EXCLUSIVE="ALIGNONLY">
              %%Keyword: NoAlign <advanced>
              <HELP>
              By default, the two starting structures are aligned to minimize the root mean square distance (RMSD) between them,
              in particular, the first structure is moved and the second structure remains fixed.
              If this keyword is given, the starting structures are used as given.
              </HELP>
              </KEYWORD>

:kword:`ALIGn only`
  The two starting structures are aligned, but nothing more is done.
  An input block for :program:`seward` is still needed, but no integrals are computed.

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="ALIGNONLY" APPEAR="Align only" KIND="SINGLE" LEVEL="ADVANCED" EXCLUSIVE="NOALIGN">
              %%Keyword: AlignOnly <advanced>
              <HELP>
              The two starting structures are aligned, but nothing more is done.
              An input block for SEWARD is still needed, but no integrals are computed.
              </HELP>
              </KEYWORD>

:kword:`WEIGhts`
  Relative weights of each atom to use for the alignment and for the calculations of the
  "distance" between structures. The possibilities are:

  .. container:: list

    **MASS**. This is the default. Each atom is given a weight proportional to its mass. Equivalent
    to mass-weighted coordinates.

    **EQUAL**. All atoms have an equal weight.

    **HEAVY**. Only heavy atoms are considered, with equal weights. Hydrogens are given zero weight.

  .. compound::

    A list of :math:`N` numbers can also be provided, and they will be used as weights for the :math:`N`
    symmetry-unique atoms. For example: ::

      WEIGhts
      0 0 0 0 0 0 1 1 1 1 1 1 0 0 0 0

    will align only atoms 7--12 out of 16.

  Note that, in any case, weights of 0 are likely to cause problems with constraints, and they will
  be increased automatically.

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="WEIGHTS" APPEAR="Weights" KIND="STRING" DEFAULT_VALUE="Mass" LEVEL="ADVANCED">
              %%Keyword: Weights <advanced>
              <HELP>
              Relative weights of each atom to use for the alignment and for the calculation of the
              "distance" between structures. The possibilities are:
              MASS: This is the default. Each atom is given a weight proportional to its mass. Equivalent to mass-weighted coordinates.
              EQUAL: All atoms have an equal weight.
              HEAVY: Only heavy atoms are considered, with equal weights. Hydrogens are given zero weight.
              A list of N numbers can also be provided, and they will be used as weights for the N symmetry-unique atoms.
              </HELP>
              </KEYWORD>

:kword:`SADDle`
  Step size reduction for each macro iteration of the saddle method.
  The value is given in weighted coordinates, divided by the square root of the total weight
  (see the :kword:`WEIGHTS` keyword).
  Default value is 0.1 au.

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="SADDLE" APPEAR="Saddle Step" KIND="REAL" DEFAULT_VALUE="0.1" LEVEL="ADVANCED">
              %%Keyword: SaddleStep <advanced>
              <HELP>
              Step size reduction for each macro iteration of the saddle method.
              The value is given in weighted coordinates, divided by the square root of the total weight
              (see the WEIGHTS keyword).
              Default value is 0.1 au.
              </HELP>
              </KEYWORD>

.. xmldoc:: </GROUP>

Geometry optimization using constrained internal coordinates
............................................................

.. compound::

  These keyword are used together with the :program:`geo` to optimize the relative position of two or
  more rigid fragments. The starting geometry can either be defined by supplying an xyz-file for each
  fragment using the keyword :kword:`coord` or by placing a file named :file:`$Project.zmt` in a directory
  named :file:`$Project.GEO`. The z-matrix should be in the following format: ::

    O     0.982011 0                                 1
    H     0.982013 0   104.959565 0                  2   1
    H     1.933697 1   107.655494 1   114.496053 1   2   3   1
    O     0.988177 0   173.057942 1   -56.200750 1   4   2   3
    H     0.979890 0   104.714572 0   179.879745 1   5   4   2

  where the three columns of real numbers are internal coordinates, and the last
  three columns of integers indicate which other atoms that are used to define
  the coordinate. The type of coordinates from left to right are bond distances,
  bond angles and dihedral angels, both for the coordinates and the link. The
  column of integers just to the right of each coordinate indicate if this
  coordinate should be optimized or not (1 = optimize, 0 = do not optimize).

There are also two utility-keywords used to create a z-matrix or to write out
a constraint-definition for :program:`slapaf` and keywords to rotate and translate
fragments. (See documentation for :program:`GEO` for more details)

.. class:: keywordlist

:kword:`HYPER`
  This keyword is used to specify that a geometry optimization with constrained
  internal coordinates shall be performed later, a z-matrix and a set of
  displaced geometries are therefore constructed. The keyword should be followed by three
  real numbers defining the maximum displacement for each coordinate type.
  The order from left to right is bond distances, bond angles and dihedral angles.
  To use default values for the parameters the mutually exclusive keyword
  :kword:`geo` should be entered instead.

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="HYPER" APPEAR="hyper" KIND="REALS" SIZE="3" LEVEL="ADVANCED" EXCLUSIVE="GEO">
              <HELP>
              Perform a geometry optimization in constrained internal coordinates using
              user-defined parameters for hypersurface gridpoints.
              </HELP>
              %%Keyword: hyper <advanced>
              Followed by three real numbers to define hypersurface gridpoint
              parameters for bond distance, bond angles and dihedral angles. Allows for a
              geometry optimization in constrained internal coordinates.
              </KEYWORD>

:kword:`GEO`
  This keyword is used to specify that a geometry optimization with constrained
  internal coordinates shall be performed later, a z-matrix and a set of displaced
  geometries are therefore constructed. Default values of 0.15 Å, 2.5 degrees,
  and 2.5 degrees are used for the maximum displacement of bond distances, bond
  angles and dihedral angles respectively. To enter other values for the parameters
  the mutually exclusive keyword :kword:`hyper` should be used.

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="GEO" APPEAR="geo" KIND="SINGLE" LEVEL="ADVANCED" EXCLUSIVE="HYPER">
              %%Keyword: geo <advanced>
              <HELP>
              Perform a geometry optimization in constrained internal coordinates using
              default parameters for hypersurface gridpoints (bond=0.15, bond angle=2.5, and
              dihedral angle=2.5)
              </HELP>
              </KEYWORD>

:kword:`OPTH`
  This keyword is used to define the specific details of the optimization algorithm used
  for the geometry optimization in constrained internal coordinates.
  This keyword should be followed by two to three lines of parameter. The first line should
  contain an integer indicating optimization type (1 = steepest descent, 2 = a mix of
  steepest descent and Newton's method, and 3 = Newton's method). The second line
  should contain a real number defining a step factor.
  This number is multiplied with the gradient to obtain the step length.
  For optimization type 2 a third line containing a real number that defines a gradient limit
  should be entered. This limit determines how large the gradient must be for the steepest
  descent algorithm to be used. When the gradient is smaller than this limit Newton's method
  is used instead.

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="OPTH" APPEAR="OptH" KIND="STRINGS" SIZE="3" LEVEL="ADVANCED">
              %%Keyword: opth <advanced>
              <HELP>
              Followed by one line with an integer specifying the optimization type (1 = steepest
              descent, 2 = mixed, 3 = Newton's method), a second line with a real number specifying
              a step factor and if using type "mixed" a third line with a real number specifying
              the maximum gradient size for which steepest descent is used.
              </HELP>
              </KEYWORD>

:kword:`OLDZ`
  This keyword is used both to start a new calculation from a user-defined z-matrix and
  to restart calculations. When using the keyword for a new calculation a directory
  :file:`$Project.GEO` must exist and contain a file called :file:`$Project.zmt` with a z-matrix in
  the format defined above. The directory must not contain any files with the suffix :file:`.info`
  when performing a fresh calculation since these files contain restart information.

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="OLDZ" APPEAR="Old Z-matrix" KIND="SINGLE" LEVEL="ADVANCED">
              %%Keyword: OldZ <advanced>
              <HELP>
              Start new calculation based on $Project.GEO/$Project.zmt
              </HELP>
              </KEYWORD>

:kword:`ZONLY`
  This keyword is used to construct a z-matrix from a set of xyz-files (fragments)
  and store it in the directory :file:`$Project.GEO`. The optimization parameters
  of the resulting z-matrix are set so that only coordinates linking fragments are
  set to 1 (= optimize coordinate).

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="ZONLY" APPEAR="z-constraints" KIND="SINGLE" LEVEL="ADVANCED">
              <HELP>
              Prints a z-matrix ($Project.zmt) in the directory $Project.GEO.
              </HELP>
              %%Keyword: zonly <advanced>
              A z-matrix ($Project.zmt) is printed in the $Project.GEO-directory. The optimization parameters
              are set so that each fragment is kept rigid and only coordinates linking fragments are
              optimized.
              </KEYWORD>

:kword:`ZCONS`
  This keyword is used to define constraints from a set of xyz-files (fragments)
  on a form that could be supplied to the
  :program:`slapaf` in order to keep the fragments rigid. The resulting constraints-file
  is named :file:`$Project.cns` and stored in the directory :file:`$Project.GEO`. The
  atom-numbers in this constraint-file will not match those of your original xyz-file and
  should not be used together with this. Instead a new xyz-file named :file:`cons.xyz` is created
  and placed into the directory :file:`$Project.GEO`, this has the proper numbering to use together with the constraints.

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="ZCONS" APPEAR="z-constraints" KIND="SINGLE" LEVEL="ADVANCED">
              <HELP>
              Prints a constraints-file ($Project.cns) and an xyz-file (cons.xyz)
              with matching atom numbering in the directory $Project.GEO .
              </HELP>
              %%Keyword: zcons <advanced>
              Prints a file with a constraints-definition for rigid fragments formatted for
              use in slapaf ($Project.cns) and an xyz-file (cons.xyz) with the same
              atom number. Both files are printed in the $Project.GEO-directory.
              </KEYWORD>

:kword:`ORIGIN`
  This keyword is used to translate and rotate a set of xyz-files. The keyword must be entered
  before the xyz-files is entered with :kword:`coord`.
  The keyword should be followed by two lines for each fragment in the input.
  The first row should contain 3 real numbers defining a translation (x, y, z),
  the second row should contain 9 numbers defining a rotation (row1, row2, row3 of
  rotation matrix). The keyword :kword:`origin` is mutually exclusive with the keyword :kword:`frgm`
  which is an alternative way to express the same rotations and translations.

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="ORIG" APPEAR="origin" KIND="UNKNOWN" LEVEL="ADVANCED" EXCLUSIVE="FRGM">
              %%Keyword: origin <advanced>
              <HELP>
              Followed by two lines for each fragment.
              The first line should have 3 real numbers defining a translation and the
              second 9 real numbers defining a rotation.
              <!--
              (See ROT and TRANS.)
              Must occur before the xyz-files are entered with coord.
              -->
              </HELP>
              </KEYWORD>

:kword:`FRGM`
  This keyword is used together with the keywords :kword:`rot` and :kword:`trans` to define
  rotation and translation of a specific fragment. :kword:`Frgm` defines an active fragment (each xyz-file is considered a fragment, the files are numbered based on
  order of appearance in the input from top to bottom). The keyword must be entered before the xyz-file it is supposed to modify is
  entered with :kword:`coord`. Each occurence of
  :kword:`frgm` should be followed by either one of or both of the keywords :kword:`rot` and :kword:`trans`
  to define rotation and translation. This keyword is mutually exclusive with the keyword :kword:`origin`

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="FRGM" APPEAR="fragment" KIND="INT" LEVEL="ADVANCED" EXCLUSIVE="ORIGIN">
              %%Keyword: frgm <advanced>
              <HELP>
              Followed by a fragment number and either or both of ROT and TRANS to define
              rotation and translation of this fragment.
              </HELP>
              Each xyz-file is considered a
              fragment, numbering is from top to bottom of input. Must occur before the modified xyz-file
              is entered with coord.
              </KEYWORD>

:kword:`ROT`
  This keyword should be followed by nine real numbers defining the rotation for the fragment defined by
  the preceeding :kword:`frgm`. The numbers should be the nine elements of a rotation matrix
  listed with one full row at the time.

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="ROT" APPEAR="rotation" KIND="REALS" SIZE="9" LEVEL="ADVANCED" REQUIRE="FRGM">
              <HELP>
              The nine numbers define a rotation matrix.
              </HELP>
              %%Keyword: rot <advanced>
              The keyword should be followed by nine real numbers defining a rotation matrix.
              Should only be used together with the FRGM keyword.
              </KEYWORD>

:kword:`TRANS`
  This keyword should be followed by three real numbers defining the translation for the fragment defined
  by the preceeding :kword:`frgm`. The numbers should be the x, y and z coordinates of the translation
  in that order.

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="TRANS" APPEAR="translation" KIND="REALS" SIZE="3" LEVEL="ADVANCED" REQUIRE="FRGM">
              <HELP>
              The three numbers define a translation. (x y z)
              </HELP>
              %%Keyword: trans <advanced>
              The keyword should be followed by three real numbers defining a translation (x y z).
              Should only be used together with the FRGM keyword.
              </KEYWORD>

Example of an input: ::

  &GATEWAY
  Title
  Water Dimer
  frgm=2
  trans=3.0 0.0 0.0
  Coord=water_monomer.xyz
  Coord=water_monomer.xyz
  Group=c1
  basis=cc-pVTZ
  hyper
  0.2 3.0 3.0
  opth
  3
  15.0d0

In this example a water dimer is constructed from a single monomer by translating
it 3.0 Å with the keyword trans. An optimization in constrained internal
coordinates using newtons method with a step-factor of 15.0d0 are prepared for. For
more details on these optimization see the manual entry for the module
:program:`geo`.

.. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="TINKER" KIND="SINGLE" LEVEL="UNDOCUMENTED" />

QM/MM calculations with |molcas|/:program:`Gromacs`
...................................................

The following keywords apply to QM/MM calculations performed with the |molcas|/:program:`GROMACS` interface (see :numref:`UG:sec:espf` for more details).

.. class:: keywordlist

:kword:`GROMacs`
  Requests that the definition of the full QM+MM system should be imported from :program:`GROMACS`. The keyword should be followed by one of the options :kword:`SIMPLE` or :kword:`CASTMM` on the next line. In the case of :kword:`SIMPLE`, all MM atoms defined in the :program:`GROMACS` input will be treated as *outer* MM atoms in |molcas|. This means, for example, that in a geometry optimization, their positions will be updated using microiterations rather than the conventional optimization scheme. Conversely, :kword:`CASTMM` requests that certain MM atoms should be treated as *inner* MM atoms in |molcas|. Their positions will be updated with the same scheme as used for the QM atoms. The :kword:`CASTMM` option should be followed by two additional input lines, the first one containing the number of MM atoms to convert from outer to inner type, and the second containing a list of those atoms (using their corresponding :program:`GROMACS` indices).

  .. xmldoc:: <GROUP MODULE="GATEWAY" NAME="GROMACS" APPEAR="Gromacs" KIND="RADIO" LEVEL="ADVANCED">
              %%Keyword: Gromacs <basic>
              Requests that the definition of the full QM+MM system should be imported from GROMACS.
              The keyword should be followed by one of the options SIMPLE or CASTMM on the next line.
              In the case of SIMPLE, all MM atoms defined in the GROMACS input will be treated as outer MM atoms in MOLCAS.
              This means, for example, that in a geometry optimization, their positions will be updated using microiterations rather than the conventional optimization scheme.
              Conversely, CASTMM requests that certain MM atoms should be treated as inner MM atoms in MOLCAS.
              Their positions will be updated with the same scheme as used for the QM atoms.
              The CASTMM option should be followed by two additional input lines,
              the first one containing the number of MM atoms to convert from outer to inner type,
              and the second containing a list of those atoms (using their corresponding GROMACS indices).

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="SIMPLE" APPEAR="Simple" KIND="SINGLE" LEVEL="UNDOCUMENTED" />

  .. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="CASTMM" APPEAR="CastMM" KIND="INTS" SIZE="2" LEVEL="UNDOCUMENTED" />

  .. xmldoc:: </GROUP>

:kword:`LINKatoms`
  Defines link atoms for use with the Morokuma updating scheme. The desired number of link atoms should be given as an integer on the next line. This should be followed by additional input lines, one for each link atom to be defined. Each definition should be of the form ILA, IQM, IMM, SCALE, where ILA, IQM and IMM are the :program:`GROMACS` indices of the link atom and the corresponding QM and MM frontier atoms, respectively. SCALE is the scaling factor to be used in the Morokuma scheme. Note that each link atom must be defined as a QM atom in the :program:`GROMACS` input. In addition, the frontier MM atom must be an inner MM atom specified as discussed above.

  .. xmldoc:: %%Keyword: LinkAtoms <advanced>
              Defines link atoms for use with the Morokuma updating scheme (MOLCAS/GROMACS calculations only).
              The desired number of link atoms should be given as an integer on the next line.
              This should be followed by additional input lines, one for each link atom to be defined.
              Each definition should be of the form ILA, IQM, IMM, SCALE, where ILA, IQM and IMM
              are the GROMACS indices of the link atom and the corresponding QM and MM frontier atoms, respectively.
              SCALE is the scaling factor to be used in the Morokuma scheme.
              Note that each link atom must be defined as a QM atom in the GROMACS input.
              In addition, the frontier MM atom must be an inner MM atom specified with the GROMACS keyword in GATEWAY.

.. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="FPCO" KIND="SINGLE" LEVEL="UNDOCUMENTED" />

.. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="FPPR" KIND="SINGLE" LEVEL="UNDOCUMENTED" />

.. xmldoc:: <KEYWORD MODULE="GATEWAY" NAME="PRINT" KIND="INTS_COMPUTED" SIZE="2" LEVEL="UNDOCUMENTED" />

.. xmldoc:: </MODULE>
