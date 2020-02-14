.. index::
   single: Program; WFA
   single: WFA

.. _UG\:sec\:wfa:

:program:`wfa`
================

.. warning::

   This program requires a submodule

.. only:: html

  .. contents::
     :local:
     :backlinks: none

.. xmldoc:: <MODULE NAME="WFA">
            %%Description:
            <HELP>
            The WFA program of the molcas program system provides various
            visual and quantitative wavefunction analysis methods.
            </HELP>

The :program:`WFA` program of the |molcas| program system provides various
visual and quantitative wavefunction analysis methods.
It is based on the libwfa :cite:`libwfa` wavefunction analysis library.
The interface to |molcas| is described in Ref. :cite:`Molcas_libwfa`.

.. Quantitative analysis methods are printed to the standard output, orbital coefficients are
   written to the HDF5 file that is used for input and output, and input files
   for the external TheoDORE program are created.

The program computes natural transition
orbitals (NTO) :cite:`Martin2003,Plasser2014`, which provide a compact description of
one-electron excited states. Natural difference orbitals (NDO) :cite:`Plasser2014` can be
computed to visualize many-body effects and orbital relaxation effects :cite:`Plasser2014b`. A module for the
statistical analysis of exciton wavefunctions is included :cite:`Bappler2014,Plasser2015`,
which provides various quantitative descriptors to describe the excited states.
Output is printed for the 1-electron transition density matrix (1TDM) and for the 1-electron difference density matrix (1DDM).
A decomposition into local and charge transfer contributions on different chromophores
is possible through the charge transfer number analysis :cite:`Plasser2012`,
which has been integrated into |molcas| recently.
Postprocessing is possible through the external `TheoDORE <http://theodore-qc.sourceforge.net/>`_ :cite:`TheoDORE` program.

Installation
------------

The :program:`WFA` module is currently not installed by default.
Its installation occurs via CMake.
It requires a working HDF5 installation and access to the include files of the Armadillo C++ linear algebra library.
In the current settings, external BLAS/LAPACK libraries have to be used.
Use, e.g., the following commands for installation: ::

  FC=ifort cmake -D LINALG=MKL -D WFA=ON -D ARMADILLO_INC=../armadillo-7.300.0/include ..

.. _UG\:sec\:wfa_dependencies:

Dependencies
------------

The :program:`WFA` program requires HDF5 files, which are written by either
:program:`SCF`, :program:`RASSCF`, or :program:`RASSI`. In the case of :program:`RASSI`,
the :kword:`TRD1` keyword has to be activated.

.. _UG\:sec\:wfa_files:

Files
-----

Input files
...........

.. class:: filelist

:file:`WFAH5`
  All information that the :program:`WFA` program needs is contained in this HDF5 file.
  The name can be adjusted with the :kword:`H5FIle` option.

Output files
............

.. class:: filelist

:file:`WFAH5`
  The orbital coefficients of NOs, NTOs, and NDOs are written to the same HDF5 file that
  is also used for input.

:file:`*.om`
  These are input files for the external TheoDORE program.

:file:`OmFrag.txt`
  Input file for TheoDORE.

For a seamless interface to TheoDORE, you can also create the :file:`tden_summ.txt` file via ::

  grep '^|' molcas.log > tden_summ.txt

Extraction of the NOs, NTOs, and NDOs from the HDF5 file occurs with the external `Molpy program <https://github.com/steabert/molpy>`_. Call, e.g.: ::

  penny molcas.rassi.h5 --wfaorbs molden

.. _UG\:sec\:wfa_input:

Input
-----

The input for the :program:`WFA` module is preceded by: ::

  &WFA

Keywords
........

Basic Keywords:

.. class:: keywordlist

:kword:`H5FIle`
  Specifies the name of the HDF5 file used for reading and writing
  (e.g. :file:`$Project.scf.h5`, :file:`$Project.rasscf.h5`, :file:`$Project.rassi.h5`).
  You either have to use this option or rename the file of
  interest to :file:`WFAH5`.

  .. xmldoc:: <KEYWORD MODULE="WFA" NAME="H5FILE" APPEAR="HDF5 file" KIND="STRING" LEVEL="BASIC">
              %%Keyword:H5FIle <basic>
              <HELP>
              Specifies the name of the HDF5 file used for reading and writing
              (e.g. $Project.scf.h5, $Project.rasscf.h5, $Project.rassi.h5).
              You either have to use this option or rename the file of
              interest to WFAH5.
              </HELP>
              </KEYWORD>

:kword:`WFALevel`
  Select how much output is produced (0-4, default: 3).

  .. xmldoc:: <KEYWORD MODULE="WFA" NAME="WFALEVEL" APPEAR="Print level" KIND="CHOICE" LIST="0,1,2,3,4" LEVEL="BASIC" DEFAULT_VALUE="3">
              %%Keyword:WFALevel <basic>
              <HELP>
              Select how much output is produced (0-4, default: 3).
              </HELP>
              </KEYWORD>

:kword:`CTNUmmode`
  Specifies what properties are computed in a `TheoDORE <http://theodore-qc.sourceforge.net/>`_-style fragment-based analysis (0-3, default: 1).
  This requires defining fragments via :kword:`ATLIsts`.

  0 --- none

  1 --- Basic: POS, PR, DEL, CT, CTnt

  2 --- Extended:  POS, POSi, POSf, PR, PRi, PRf, DEL, COH, CT, CTnt

  3 --- For transition metal complexes: POSi, POSf, PR, CT, MC, LC, MLCT, LMCT, LLCT

  The definition of the descriptors is provided
  `here <https://sourceforge.net/p/theodore-qc/wiki/Transition%20density%20matrix%20analysis/attachment/Om_desc.pdf>`_.
  For a more fine-grained input use :kword:`PROPlist`.

  .. xmldoc:: <KEYWORD MODULE="WFA" NAME="CTNUMMODE" APPEAR="Computed properties" KIND="CHOICE" LIST="0: None,1: Basic,2: Extended,3: Metal complexes" LEVEL="BASIC" DEFAULT_VALUE="1" REQUIRE="ATLISTS">
              %%Keyword:CTNUmmode <basic>
              <HELP>
              Define what properties are computed in a TheoDORE-style analysis. (0-3, default: 1).
              </HELP>
              </KEYWORD>

:kword:`ATLIsts`
  Define the fragments in a `TheoDORE <http://theodore-qc.sourceforge.net/>`_-style analysis.

  The first entry is the number of fragments.
  Then enter the atomic indices of the fragment followed by a \*.
  Example: ::

    ATLISTS
    2
    1 2 4 *
    3 *

  .. xmldoc:: <KEYWORD MODULE="WFA" NAME="ATLISTS" APPEAR="Fragment definition" KIND="CUSTOM" LEVEL="BASIC">
              %%Keyword:ATLIsts <basic>
              <HELP>
              Define the fragments in a TheoDORE-style analysis.
              </HELP>
              </KEYWORD>

:kword:`REFState`
  Index of the reference state for 1TDM and 1DDM analysis (default: 1).

  .. xmldoc:: <KEYWORD MODULE="WFA" NAME="REFSTATE" APPEAR="Reference state" KIND="INT" LEVEL="BASIC" DEFAULT_VALUE="1">
            : %%Keyword:REFState <basic>
              <HELP>
              Index of the reference state for 1TDM and 1DDM analysis.
              </HELP>
              </KEYWORD>

Advanced keywords for fine grain output options and debug information:

.. class:: keywordlist

:kword:`MULLiken`
  Activate Mulliken population analysis (also for CT numbers).

  .. xmldoc:: <KEYWORD MODULE="WFA" NAME="MULLIKEN" APPEAR="Mulliken population analysis" KIND="SINGLE" LEVEL="ADVANCED">
              %%Keyword:MULLiken <advanced>
              <HELP>
              Activate Mulliken population analysis.
              </HELP>
              </KEYWORD>

:kword:`LOWDin`
  Activate Löwdin population analysis (also for CT numbers).

  .. xmldoc:: <KEYWORD MODULE="WFA" NAME="LOWDIN" APPEAR="Lowdin population analysis" KIND="SINGLE" LEVEL="ADVANCED">
              %%Keyword:LOWDin <advanced>
              <HELP>
              Activate Lowdin population analysis.
              </HELP>
              </KEYWORD>

:kword:`NXO`
  Activate NO, NTO, and NDO analysis.

  .. xmldoc:: <KEYWORD MODULE="WFA" NAME="NXO" APPEAR="NXO analysis" KIND="SINGLE" LEVEL="ADVANCED">
              %%Keyword:NXO <advanced>
              <HELP>
              Activate NO, NTO, and NDO analysis.
              </HELP>
              </KEYWORD>

:kword:`EXCIton`
  Activate exciton and multipole analysis.

  .. xmldoc:: <KEYWORD MODULE="WFA" NAME="EXCITON" APPEAR="Exciton analysis" KIND="SINGLE" LEVEL="ADVANCED">
              %%Keyword:EXCIton <advanced>
              <HELP>
              Activate exciton and multipole analysis.
              </HELP>
              </KEYWORD>

:kword:`DOCTnumbers`
  Activate charge transfer number analysis and creation of :file:`*.om` files.

  .. xmldoc:: <KEYWORD MODULE="WFA" NAME="DOCTNUMBERS" APPEAR="Charge transfer numbers" KIND="SINGLE" LEVEL="ADVANCED">
              %%Keyword:DOCTnumbers <advanced>
              <HELP>
              Activate charge transfer number analysis and creation of \*.om files.
              </HELP>
              </KEYWORD>

:kword:`H5ORbitals`
  Print the NOs, NTOs, and/or NDOs to the HDF file.

  .. xmldoc:: <KEYWORD MODULE="WFA" NAME="H5ORBITALS" APPEAR="Save orbitals in HDF5" KIND="SINGLE" LEVEL="ADVANCED">
              %%Keyword:H5ORbitals <advanced>
              <HELP>
              Print the NOs, NTOs, and/or NDOs to the HDF file.
              </HELP>
              </KEYWORD>

:kword:`PROPlist`
  Manual input of properties to be printed out in a `TheoDORE <http://theodore-qc.sourceforge.net/>`_-style fragment based analysis.
  Use only if :kword:`CTNUMMODE` does not provide what you want.

  .. xmldoc:: <KEYWORD MODULE="WFA" NAME="PROPLIST" APPEAR="Property list" KIND="CUSTOM" LEVEL="ADVANCED">
              %%Keyword:PROPlist <advanced>
              <HELP>
              Manual input of properties to be printed out in a TheoDORE-style analysis.
              </HELP>
              </KEYWORD>

  Enter as a list followed by a \*, e.g. ::

    PROPLIST
    Om POS PR CT COH CTnt *

  The full list of descriptors is provided
  `here <https://sourceforge.net/p/theodore-qc/wiki/Transition%20density%20matrix%20analysis/attachment/Om_desc.pdf>`_.

:kword:`DEBUg`
  Print debug information.

  .. xmldoc:: <KEYWORD MODULE="WFA" NAME="DEBUG" APPEAR="Print debug information" KIND="SINGLE" LEVEL="ADVANCED">
              %%Keyword:DEBUg <advanced>
              <HELP>
              Print debug information.
              </HELP>
              </KEYWORD>

:kword:`ADDInfo`
  Add info for verification runs with :command:`molcas verify`.

  .. xmldoc:: <KEYWORD MODULE="WFA" NAME="ADDINFO" APPEAR="Add info" KIND="SINGLE" LEVEL="ADVANCED">
              %%Keyword:ADDInfo <advanced>
              <HELP>
              Add info for verifications runs with molcas verify.
              </HELP>
              </KEYWORD>

Input example
.............

::

  * Analysis of SCF job
  &SCF

  &WFA
  H5file = $Project.scf.h5

::

  * Analysis of RASSCF job
  * Reduced output
  &RASSCF

  &WFA
  H5file = $Project.rasscf.h5
  wfalevel = 1

::

  * Analysis of RASSI job, use the TRD1 keyword
  &RASSI
  TRD1

  &WFA
  H5file = $Project.rassi.h5
  ATLISTS
  2
  1 2 4 *
  3 *

.. _UG\:sec\:wfa_output:

Output
------

State/difference density matrix analysis (:program:`SCF`/:program:`RASSCF`/:program:`RASSI`)
............................................................................................

.. compound::

  ::

    RASSCF analysis for state 2 (3) A

  or ::

    RASSI analysis for state R_2

.. _tab\:wfa_dm:

======================= ===================================================================================================
Descriptor              Explanation
======================= ===================================================================================================
``n_u``                 Number of unpaired electrons :math:`n_u=\sum_i\min(n_i, 2-n_i)` :cite:`Head-Gordon2003,Plasser2014`
``n_u,nl``              Number of unpaired electrons :math:`n_{u,nl}=\sum_i n_i^2(2-n_i)^2`
``PR_NO``               NO participation ratio :math:`\text{PR}_{\text{NO}}`
``p_D`` and ``p_A``     Promotion number :math:`p_D` and :math:`p_A`
``PR_D`` and ``PR_A``   D/A participation ratio :math:`\text{PR}_D` and :math:`\text{PR}_A`
``<r_h> [Ang]``         Mean position of detachment density :math:`\vec{d}_D` :cite:`Plasser2015`
``<r_e> [Ang]``         Mean position of attachment density :math:`\vec{d}_A`
``|<r_e - r_h>| [Ang]`` Linear D/A distance :math:`\vec{d}_{D\rightarrow A} = \vec{d}_A - \vec{d}_D`
``Hole size [Ang]``     RMS size of detachment density :math:`\sigma_D`
``Electron size [Ang]`` RMS size of attachment density :math:`\sigma_A`
======================= ===================================================================================================

Transition density matrix analysis (:program:`RASSI`)
.....................................................

::

  RASSI analysis for transiton from state 1 to 2 (Tr_1-2)

.. _tab\:wfa_tdm:

====================================== =============================================================================================================================
Output listing                         Explanation
====================================== =============================================================================================================================
``Leading SVs``                        Largest NTO occupation numbers
``Sum of SVs (Omega)``                 :math:`\Omega`, Sum of NTO occupation numbers
``PR_NTO``                             NTO participation ratio :math:`\text{PR}_{\text{NTO}}` :cite:`Plasser2012`
``Entanglement entropy (S_HE)``        :math:`S_{H|E}=-\sum_i\lambda_i\log_2\lambda_i` :cite:`Plasser2016`
``Nr of entangled states (Z_HE)``      :math:`Z_{HE}=2^{S_{H|E}}`
``Renormalized S_HE/Z_HE``             Replace :math:`\lambda_i\rightarrow \lambda_i/\Omega`
``omega``                              Norm of the 1TDM :math:`\Omega`, single-exc. character
``<Phe>``                              Exp. value of the particle-hole permutation operator, measuring de-excitations :cite:`Kimber2020`
``<r_h> [Ang]``                        Mean position of hole :math:`\langle\vec{x}_h\rangle_{\text{exc}}` :cite:`Plasser2015`
``<r_e> [Ang]``                        Mean position of electron :math:`\langle\vec{x}_e\rangle_{\text{exc}}`
``|<r_e - r_h>| [Ang]``                Linear e/h distance :math:`\vec{d}_{h\rightarrow e} = \langle\vec{x}_e - \vec{x}_h\rangle_{\text{exc}}`
``Hole size [Ang]``                    RMS hole size: :math:`\sigma_h = (\langle\vec{x}_h^2\rangle_{\text{exc}} - \langle\vec{x}_h\rangle_{\text{exc}}^2)^{1/2}`
``Electron size [Ang]``                RMS electron size: :math:`\sigma_e = (\langle\vec{x}_e^2\rangle_{\text{exc}} - \langle\vec{x}_e\rangle_{\text{exc}}^2)^{1/2}`
``RMS electron-hole separation [Ang]`` :math:`d_{\text{exc}} = (\langle \left|\vec{x}_e - \vec{x}_h\right|^2\rangle_{\text{exc}})^{1/2}` :cite:`Bappler2014`
``Covariance(r_h, r_e) [Ang^2]``       :math:`\text{COV}\left(\vec{x}_h,\vec{x}_e\right) = \langle\vec{x}_h\cdot\vec{x}_e\rangle_{\text{exc}} -
                                       \langle\vec{x}_h\rangle_{\text{exc}}\cdot\langle\vec{x}_e\rangle_{\text{exc}}`
``Correlation coefficient``            :math:`R_{eh} = \text{COV}\left(\vec{x}_h,\vec{x}_e\right)/\sigma_h\cdot\sigma_e` :cite:`Plasser2015`
``Center-of-mass size``                :math:`(\langle \left|\vec{x}_e + \vec{x}_h\right|^2\rangle_{\text{exc}}-\langle \vec{x}_e + \vec{x}_h\rangle_{\text{exc}}^2)^{1/2}`
====================================== =============================================================================================================================

.. xmldoc:: </MODULE>
