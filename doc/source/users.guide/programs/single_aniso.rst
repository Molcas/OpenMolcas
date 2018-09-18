.. index::
   single: Program; Single_aniso
   single: Single_aniso

.. _UG\:sec\:single_aniso:

:program:`single_aniso` |extramark|
===================================

.. warning::

   This program is not available in OpenMolcas

.. only:: html

.. contents::
    :local:
    :backlinks: none

.. xmldoc:: <MODULE NAME="SINGLE_ANISO" APPEAR="Single_Aniso">
            %%Description:
            <HELP>
            The SINGLE_ANISO program allows the non-perturbative calculation of
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
            ||   c) Higher order zero-field splitting tensors (D^2, D^4, D^6, ... etc.)
            ||   d) Higher order Zeeman splitting tensors (G^1, G^3, G^5, ... , etc.)
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

The :program:`SINGLE_ANISO` program is a routine which allows the non-perturbative calculation of effective spin (pseudospin) Hamiltonians and static magnetic properties of mononuclear complexes and fragments completely *ab initio*, including the spin-orbit interaction. As a starting point it uses the results of :program:`RASSI` calculation for the ground and several excited spin-orbital multiplets. A short description of methodology and applications can be found in :cite:`Chibotaru:1`, :cite:`Chibotaru:2`. The second version of the :program:`SINGLE_ANISO` program is able to calculate the following quantities:

* Parameters of pseudospin magnetic Hamiltonians (the methodology is described in :cite:`Chibotaru:3`):

  #. First rank (linear after pseudospin) Zeeman splitting tensor :math:`g_{\alpha\beta}`, its main values, including the sign of the product :math:`g_{X} \cdot g_{Y} \cdot g_{Z}`, and the main magnetic axes.
  #. Second rank (bilinear after pseudospin) zero-field splitting tensor :math:`D_{\alpha\beta}`, its main values and the anisotropy axes. The anisotropy axes are given in two coordinate systems: a) in the initial Cartesian coordinate system (:math:`x, y, z`) and b) in the coordinate system of the main magnetic axes (:math:`X_{\text{m}}, Y_{\text{m}}, Z_{\text{m}}`).
  #. Higher rank ZFS tensors (:math:`D^4`, :math:`D^6`, etc.) and Zeeman splitting tensors (:math:`G^3`, :math:`G^5`, etc.) for complexes with moderate and strong spin-orbit coupling.
  #. Angular moments along the main magnetic axes.

* All (27) parameters of the *ab initio* Crystal field acting on the ground atomic multiplet of lanthanides, and the decomposition of the CASSCF/RASSI wave functions into functions with definite projections of the total angular moment on the quantization axis.

* Static magnetic properties:

  #. Van Vleck susceptibility tensor :math:`\chi_{\alpha\beta}(T)`.
  #. Powder magnetic susceptibility function :math:`\chi(T)`.
  #. Magnetization vector :math:`\vec M (\vec H)` for specified directions of the applied magnetic field :math:`\vec H`.
  #. Powder magnetization :math:`M_{\text{mol}}(H)`.

The magnetic Hamiltonians are defined for a desired group of :math:`N` electronic states obtained in :program:`RASSI` calculation to which a pseudospin :math:`\tilde{S}` (it reduces to a true spin :math:`S` in the absence of spin-orbit coupling) is subscribed according to the relation :math:`N=2\tilde{S}+1`. For instance, the two wave functions of a Kramers doublet correspond to :math:`\tilde{S}=1/2`. The implementation is done for :math:`\tilde{S}=1/2, 1, 3/2, \ldots ,15/2`.

.. The second version of the :program:`SINGLE_ANISO` program allows the calculation of all 27 parameters of the exact Crystal-Field acting on the ground atomic multiplet for lanthanides. Moreover, the *ab initio* wave functions corresponding to the lowest atomic multiplet :math:`\ket{J,M_J}` are decomposed in a linear combination of functions with definite projection of the total moment on the quantization axis.

The calculation of magnetic properties takes into account the contribution of excited states (the ligand-field and charge transfer states of the complex or mononuclear fragment included in the RASSI calculation) via their thermal population and Zeeman admixture. The intermolecular exchange interaction between magnetic molecules in a crystal can be taken into account during the simulation of magnetic properties by a phenomenological parameter :math:`z_J` specified by the user (see keyword MLTP).

.. index::
   pair: Dependencies; Single_aniso

.. _UG\:sec\:single_aniso_dependencies:

Dependencies
------------

The :program:`SINGLE_ANISO` program takes all needed *ab initio* information from the :file:`RUNFILE`: i.e. matrix elements of angular momentum, spin-orbit energy spectrum and mixing coefficients, number of mixed states and their multiplicity, etc. In order to find the necessary information in the :file:`RUNFILE`, the keywords MEES and SPIN are mandatory for :program:`RASSI`. The :program:`SEWARD` keyword ANGM is also compulsory.

.. index::
   pair: Files; Single_aniso

.. _UG\:sec\:single_aniso_files:

Files
-----

Input files
...........

.. class:: filelist

:file:`RUNFILE`
  The file of communication between different modules in |molcas|. Its presence is mandatory when the calculation is not a restart one from a data file.

Restart files & options
.......................

.. class:: filelist

:file:`RUNFILE`
  The file of communication between different modules in |molcas|. Normally it is already present in :file:`i$WorkDir`.
  The :program:`SINGLE_ANISO` may be restarted as many times as necessary in the same working directory where the previous :program:`RASSI` was succesfully executed. The :file:`RUNFILE` contains then all necessary data.

:file:`ANISOINPUT`
  The program may be restarted from the ASCII text file :file:`ANISOINPUT` generated by a previous succesful run of the :program:`SINGLE_ANISO` (the name of this file may be specified during execution, see :kword:`REST` keyword below). This file contains all necessary data for :program:`SINGLE_ANISO` as well as for the :program:`POLY_ANISO`. In this case the initial :file:`$WorkDir` may be empty (:file:`RUNFILE` is not necessary).

:file:`$Project.aniso`
  The :program:`SINGLE_ANISO` may be restarted from the binary file :file:`$Project.aniso` produced in a previous run. The initial :file:`$WorkDir` may be empty (:file:`RUNFILE` is not necessary).

Output files
............

.. class:: filelist

:file:`$Project.aniso`
  This binary file may be used for restart. It is produced by any successful run of the code.

:file:`ANISOINPUT`
  This file is intended to be as input for the :program:`POLY_ANISO` module in |molcas|. It is an ASCII formated file. It is produced by any successful run of the code.

:file:`zeeman_energy_xxx.txt`
  Zeeman eignestates for the applied field in the direction # *xxx* are placed in the corresponding text file. It may be used directly with external plotting programs like gnuplot to visualize the data.

:file:`XT_compare.txt`
  In case :kword:`TEXP` is employed (experimental XT(T) data points), the :program:`SINGLE_ANISO` produces a data file used to directly plot the comparison between experimental and calculated magnetic susceptibility.

:file:`MH_compare_xxx.txt`
  In case :kword:`HEXP` is employed (experimental M(H,T) data points), the :program:`SINGLE_ANISO` produces one or several data file(s) used to directly plot the comparison(s) between experimental and calculated molar magnetization at each temperature.

.. index::
   pair: Input; Single_aniso

.. _UG\:sec\:single_aniso_input:

Input
-----

Normally :program:`SINGLE_ANISO` runs without specifying any of the following keywords. The only unknown variable for :program:`SINGLE_ANISO` is the dimension (multiplicity) of the pseudospin. By default one multiplet is selected, which has the dimension equal to the multiplicity of the ground term. For example, in cases where spin-orbit coupling is weak, the multiplicity of the effective spin Hamiltonian is usually the same as the multiplicity of the lowest term, while in the cases with strong anisotropy (lanthanide or actinide complexes, :math:`\ce{Co^{2+}}` complexes, etc...) the lowest energy levels of the complexes form a group of states which can differ quite strong from the spin multiplicity of the lowest term. In these cases the user should specify the multiplicity corresponding to a chosen value of pseudospin :math:`(2\tilde{S}+1)`. For instance, in :math:`\ce{Dy^{3+}}` the spin of the ground state term is :math:`S=5/2`, but in many situations only the ground Kramers doublet is considered; then the user should set the multiplicity of the pseudospin equal to 2 (see MLTP keyword).
The calculation of the parameters of the crystal field corresponding to the ground atomic multiplet for lanthanides should be requested by the CRYS keyword. ::

  &SINGLE_ANISO

Argument(s) to a keyword are always supplied on the next line of the
input file.

Optional general keywords to control the input
..............................................

.. class:: keywordlist

:kword:`TITLe`
  One line following this one is regarded as title.

  .. xmldoc:: <KEYWORD MODULE="SINGLE_ANISO" NAME="TITLE" KIND="STRING" LEVEL="BASIC">
              %%Keyword: TITLE <basic>
              <HELP>
              One line following this one is regarded as title.
              </HELP>
              </KEYWORD>

:kword:`TYPE`
  This keyword is obsolete

  .. xmldoc:: <KEYWORD MODULE="SINGLE_ANISO" NAME="TYPE" KIND="INT" LEVEL="BASIC">
              %%Keyword: TYPE <basic>
              <HELP>
              This keyword is obsolete
              </HELP>
              </KEYWORD>

:kword:`MLTP`
  The number of molecular multiplets (i.e. groups of spin-orbital eigenstates) for
  which :math:`g`, :math:`D` and higher magnetic tensors will be calculated (default :kword:`MLTP`\=1).
  The program reads two lines: the first is the number of multiplets (:math:`n_{\text{mult}}`) and
  the second the array of :math:`n_{\text{mult}}` numbers specifying the dimension of each multiplet.
  By default, the code will first analyze the energy spectra by itself and will
  compute the :math:`g` and :math:`D` tensors for ten low-lying groups of states. By using this
  keyword the user overwrites the default.

  Example: ::

    MLTP
    4
    4 4 2 2

  :program:`SINGLE_ANISO` will compute the :math:`g` tensor for four groups of states:
  the first two groups having the effective spin :math:`\tilde{S}=\ket{3/2}` each, while
  the other two groups of states being Kramers doublets.

  .. xmldoc:: <KEYWORD MODULE="SINGLE_ANISO" NAME="MLTP" KIND="INT" LEVEL="BASIC" DEFAULT_VALUE="1">
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

  :program:`SINGLE_ANISO` will compute temperature dependence of the magnetic susceptibility in 331 points evenly distributed in temperature interval: 0.0 K -- 330.0 K.

  .. xmldoc:: <KEYWORD MODULE="SINGLE_ANISO" NAME="TINT" KIND="REAL" LEVEL="BASIC">
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
  Specifies the field points for the evaluation of the magnetization in a certain direction. The program will read four numbers: :math:`H_{\text{min}}`, :math:`H_{\text{max}}`, :math:`n_H`.

  .. container:: list

    :math:`H_{\text{min}}` --- the minimal field (Default 0.0 T)

    :math:`H_{\text{max}}` --- the maximal filed (Default 10.0 T)

    :math:`n_H` --- number of field points (Default 101)

  Example: ::

    HINT
    0.0  20.0  201

  :program:`SINGLE_ANISO` will compute the molar magnetization in 201 points evenly distributed in field interval: 0.0 T -- 20.0 T.

  .. xmldoc:: <KEYWORD MODULE="SINGLE_ANISO" NAME="HINT" KIND="REAL" LEVEL="BASIC">
              %%Keyword: HINT <basic>
              <HELP>
              Specifies the field points for the evaluation of the molar magnetization.
              The program will read four numbers: Hmin, Hmax, nH, and dltH0. Units of magnetic field = Tesla (T).
              ||Hmin  -- the minimal field (Default 0.0 T)
              ||Hmax  -- the maximal field (Default 300.0 T)
              ||nH    -- number of field points (Default 101)
              </HELP>
              </KEYWORD>

:kword:`TMAG`
  Specifies the temperature(s) at which the field-dependent magnetization is calculated. The program will read the number of temperature points (:math:`N_{\text{temp}}`) and then an array of real numbers specifying the temperatures (in kelvin) at which magnetization is to be computed.
  Default is to compute magnetization at one temperature point (2.0 K).
  Example: ::

    TMAG
    5   1.8  2.0  3.4  4.0  5.0

  :program:`SINGLE_ANISO` will compute the molar magnetization at 5 temperature points (1.8 K, 2.0 K, 3.4 K, 4.0 K, and 5.0 K).

  .. xmldoc:: <KEYWORD MODULE="SINGLE_ANISO" NAME="TMAG" KIND="REAL" LEVEL="BASIC">
              %%Keyword: TMAG <basic>
              <HELP>
              Specifies the temperature(s) at which the field-dependent magnetization is calculated. The program will read the number of temperature points (NTemp) and then an array of
              </HELP>
              </KEYWORD>

:kword:`ENCU`
  This flag is used to define the cut-off energy for the lowest states for which
  Zeeman interaction is taken into account exactly. The contribution to the magnetization
  arising from states that are higher in energy than :math:`E` (see below) is done by
  second-order perturbation theory. The program will read two integer
  numbers: :math:`N_K` and :math:`M_G`. Default values are: :math:`N_K`\=100, :math:`M_G`\=100.

  .. math:: E=N_K \cdot k_{\text{Boltz}} \cdot \text{TMAG} + M_G \cdot \mu_{\text{Bohr}} \cdot H_{\text{max}}

  The field-dependent magnetization is calculated at the temperature value TMAG.
  Example: ::

    ENCU
    250  150

  If :math:`H_{\text{max}}` = 10 T and :kword:`TMAG` = 1.8 K, then the cut-off energy is:

  .. math:: E=100 \cdot 250 \cdot k_{\text{Boltz}} \cdot 1.8\,\text{K} + 150 \cdot \mu_{\text{Bohr}} \cdot 10\,\text{T} = 1013.06258\,\text{cm}^{-1}

  This means that the magnetization coming from all spin-orbit states with energy lower
  than :math:`E=1013.06258\,\text{cm}^{-1}` will be computed exactly. The contribution from the
  spin-orbit states with higher energy is accounted by second-order perturbation.

  .. xmldoc:: <KEYWORD MODULE="SINGLE_ANISO" NAME="ENCU" KIND="INT" LEVEL="BASIC">
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

:kword:`NCUT`
  This flag is used to define the cut-off energy for the lowest states for which
  Zeeman interaction is taken into account exactly. The contribution to the magnetization
  arising from states that are higher in energy than lowest :math:`N_{\text{CUT}}` states, is done by
  second-order perturbation theory. The program will read one integer number. In case the number
  is larger than the total number of spin-orbit states(:math:`N_{\text{SS}}`, then the :math:`N_{\text{CUT}}` is set to :math:`N_{\text{SS}}`
  (which means that the molar magnetization will be computed exactly, using full Zeeman
  diagonalization for all field points). The field-dependent magnetization is calculated at
  the temperature value(s) defined by :kword:`TMAG`.

  Example: ::

    NCUT
    32

  .. xmldoc:: <KEYWORD MODULE="SINGLE_ANISO" NAME="NCUT" KIND="INT" LEVEL="BASIC">
              %%Keyword: NCUT <basic>
              <HELP>
              This keyword is used to define the cut-off energy for the lowest states for which
              Zeeman interaction is taken into account exactly. The contribution to the
              magnetization coming from states that are higher in energy than E (see below)
              is done by second order perturbation theory. The program will read two integer
              numbers: NK and MG. The field-dependent magnetization
              is calculated at the temperature value TMAG.
              </HELP>
              </KEYWORD>

:kword:`MVEC`
  Defines the number of directions for which the magnetization vector will be computed.
  On the first line below the keyword, the number of directions should be mentioned (NDIR. Default 0).
  The program will read NDIR lines for cartesian coordinates specifying the direction :math:`i` of the
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

  .. xmldoc:: <KEYWORD MODULE="SINGLE_ANISO" NAME="MVEC" KIND="REAL" LEVEL="BASIC">
              %%Keyword: MVEC <basic>
              <HELP>
              Defines the number of directions for which the magnetization vector will be computed.
              On the first line below the keyword, the number of directions should be mentioned (NDIR. Default 0).
              The program will read NDIR lines for spherical coordinates specifying the direction
              "i" of the magnetic field (theta_i and phi_i). These values should be in radians.
              </HELP>
              </KEYWORD>

:kword:`MAVE`
  This keyword specifies the grid density used for the computation of powder molar
  magnetization. The program uses Lebedev-Laikov distribution of points on the unit sphere.
  The program reads two integer numbers: :math:`n_{\text{sym}}` and :math:`n_{\text{grid}}`. The :math:`n_{\text{sym}}` defines which
  part of the sphere is used for averaging. It takes one of the three values: 1 (half-sphere),
  2 (a quater of a sphere) or 3 (an octant of the sphere). :math:`n_{\text{grid}}` takes values from 1
  (the smallest grid) till 32 (the largest grid, i.e. the densiest). The default is to
  consider integration over a half-sphere (since :math:`M(H)=-M(-H)`): :math:`n_{\text{sym}}=1` and :math:`n_{\text{sym}}=15`
  (i.e 185 points distributed over half-sphere). In case of symmetric compounds, powder
  magnetization may be averaged over a smaller part of the sphere, reducing thus the number
  of points for the integration. The user is responsible to choose the appropriate integration scheme.
  Note that the program's default is rather conservative.

  .. container:: list

    :math:`N_\theta` --- number of :math:`\theta` points in the interval :math:`(0, \pi/2)`. (Default 12)

    :math:`N_\phi` --- number of :math:`\phi` points in the interval :math:`(0, 2\pi)`. (Default 24)

  The number of directions over which the actual averaging will take place is roughly the product of :math:`N_\theta` and :math:`N_\phi`.

  .. xmldoc:: <KEYWORD MODULE="SINGLE_ANISO" NAME="MAVE" KIND="INT" LEVEL="BASIC">
              %%Keyword: MAVE <basic>
              <HELP>
              This keyword specifies the grid density used for the computation of powder molar
              magnetization. The program uses Lebedev-Laikov distribution of points on the unit sphere.
              The program reads two integer numbers: NSYM and NGRID. The NSYM defines which
              part of the sphere is used for averaging. It takes one of the three values: 1 (half-sphere),
              2 (a quater of a sphere) or 3 (an octant of the sphere). NGRID takes values from 1
              (the smallest grid) till 32 (the largest grid, i.e. the densiest). The default is to
              consider integration over a half-sphere (since M(H)=-M(-H)): NSYM=1 and NGRID=15
              (i.e 185 points distributed over half-sphere). In case of symmetric compounds, powder
              magnetization may be averaged over a smaller part of the sphere, reducing thus the number
              of points for the integration. The user is responsible to choose the appropriate integration scheme.
              Note that the program's default is rather conservative.
              </HELP>
              </KEYWORD>

:kword:`TEXP`
  This keyword allows computation of the magnetic susceptibility :math:`\chi T(T)` at experimental points.
  On the line below the keyword, the number of experimental points :math:`N_T` is defined, and on
  the next :math:`N_T` lines the program reads the experimental temperature (in kelvin) and the
  experimental magnetic susceptibility (in :math:`\text{cm}^3\,\text{K}\,\text{mol}^{-1}`).
  :kword:`TEXP` and :kword:`TINT` keywords are mutually exclusive. The magnetic susceptibility
  routine will also print the total average standard deviation from the experiment.

  .. xmldoc:: <KEYWORD MODULE="SINGLE_ANISO" NAME="TEXP" KIND="REAL" LEVEL="BASIC">
              %%Keyword: TEXP <basic>
              <HELP>
              This keyword allows computation of the magnetic susceptibility at experimental
              temperature points. On the line below the keyword, the number of experimental
              points NT is defined, and on the next NT lines the program reads the experimental
              temperature (in K) and the experimental magnetic susceptibility (in cm^3Kmol^{-1} ).
              TEXP and TINT keywords are mutually exclusive. The SINGLE_ANISO will also print the
              standard deviation from the experiment.
              </HELP>
              </KEYWORD>

:kword:`HEXP`
  This keyword allows computation of the molar magnetization :math:`M_{\text{mol}} (H)` at experimental points.
  On the line below the keyword, the number of experimental points :math:`N_H` is defined, and on the next :math:`N_H` lines
  the program reads the experimental field strength (in tesla) and the experimental magnetization (in :math:`\mu_{\text{Bohr}}`).
  :kword:`HEXP` and :kword:`HINT` are mutually exclusive. The magnetization routine will print the standard deviation from the experiment.

  .. xmldoc:: <KEYWORD MODULE="SINGLE_ANISO" NAME="HEXP" KIND="REAL" LEVEL="BASIC">
              %%Keyword: HEXP <basic>
              <HELP>
              This keyword allows computation of the molar magnetization at experimental field points.
              On the line below the keyword,the number of experimental points NH is defined, and on
              the next NH lines the program reads the experimental field strength (Tesla) and the
              experimental magnetization (in Bohr magnetons). HEXP and HINT are mutually exclusive.
              The SINGLE_ANISO will print the standard deviation from the experiment.
              </HELP>
              </KEYWORD>

:kword:`ZJPR`
  This keyword specifies the value (in :math:`\text{cm}^{-1}`) of a phenomenological parameter of a mean
  molecular field acting on the spin of the complex (the average intermolecular exchange
  constant). It is used in the calculation of all magnetic properties (not for pseudo-spin
  Hamiltonians) (Default is 0.0)

  .. xmldoc:: <KEYWORD MODULE="SINGLE_ANISO" NAME="ZJPR" KIND="REAL" LEVEL="BASIC">
              %%Keyword: ZJPR <basic>
              <HELP>
              This keyword specifies the value (in cm^-1) of a phenomenological parameter of a
              mean molecular field acting on the spin of the complex (the average intermolecular
              exchange constant). It is used in the calculation of all magnetic properties (not for
              spin Hamiltonians) (Default is 0.0)
              </HELP>
              </KEYWORD>

:kword:`PRLV`
  This keyword controls the print level.

  .. container:: list

    2 --- normal. (Default)

    3 or larger (debug)

  .. xmldoc:: <KEYWORD MODULE="SINGLE_ANISO" NAME="PRLV" KIND="INT" LEVEL="BASIC">
              %%Keyword: PRLV <basic>
              <HELP>
              This keyword controls the print level.
              ||2 -- normal. (Default)
              ||3 or larger (debug)
              </HELP>
              </KEYWORD>

:kword:`POLY`
  The keyword is obsolete. The :program:`SINGLE_ANISO` creates by default one ASCII formated text file named :file:`ANISOINPUT`
  and also a binary file named :file:`$Project.Aniso`. Both may be used to restart (or re-run again) the :program:`SINGLE_ANISO` calculation.

  .. xmldoc:: <KEYWORD MODULE="SINGLE_ANISO" NAME="POLY" KIND="STRING" LEVEL="ADVANCED">
              %%Keyword: POLY <basic>
              <HELP>
              SINGLE_ANISO will prepare an input file (binary) for the future POLY_ANISO program. The default is not to create it.
              </HELP>
              </KEYWORD>

:kword:`CRYS`
  This keyword will enables the computation of the parameters of the crystal-field acting on the ground atomic multiplet of a
  lanthanide from the *ab initio* calculation performed. The implemented methodology is described :cite:`Ungur2017` and :cite:`Chibotaru:3`.
  Two types of crystal field parametererization are implemented:

  1. Parameterisation of the ground :math:`\ket{J,M_J}` group of spin-orbit states (e.g. parameterisation of the ground :math:`J=15/2` of a :math:`\ce{Dy^{3+}}` complex).
  2. Parameterisation of the ground :math:`\ket{L,M_L}` group of spin-free states (e.g. parameterisation of the ground :math:`^6H` multiplet of a :math:`\ce{Dy^{3+}}`).

  For each of the above cases, the parameters of the crystal field are given in terms of irreducible tensor
  operators defined in :cite:`Chibotaru:3`, in terms of Extended Stevens Operators defined in :cite:`Rudowicz1985` and also
  employed in the EasySpin function of MATLAB.
  On the next line the program will read the chemical symbol of the metal ion.
  The code understands the labels of: lanthanides, actinides and first-row transition metal ions. For transition metal ions, the oxidation state
  should be indicated as well.
  By default the program will not compute the parameters of the crystal-field.

  .. xmldoc:: <KEYWORD MODULE="SINGLE_ANISO" NAME="CRYS" KIND="STRING" LEVEL="BASIC">
              %%Keyword: CRYS <basic>
              <HELP>
              This keyword will enable computation of all 27 Crystal-Field parameters acting on the ground atomic multiplet of a lanthanide. On the next line the program wil read the chemical symbol of the lanthanide. By default the program will not compute the parameters of the Crystal Field.
              </HELP>
              </KEYWORD>

:kword:`QUAX`
  This keyword controls the quantization axis for the computation of the Crystal-Field parameters acting on the ground atomic multiplet of a lanthanide. On the next line, the program will read one of the three values: 1, 2 or 3.

  .. container:: list

    1 --- quantization axis is the main magnetic axis :math:`Z_{\text{m}}` of the ground pseudospin multiplet, whose size is specified within the :kword:`MLTP` keyword. (Default)

    2 --- quantization axis is the main magnetic axis :math:`Z_{\text{m}}` of the entire atomic multiplet :math:`\ket{J,M_J}`.

    3 --- the direction of the quantization axis is given by the user: on the next line the program will read three real numbers: the projections (:math:`p_x`, :math:`p_y`, :math:`p_z`) of the specified direction on the initial Cartesian axes. Note that :math:`p_x^2 + p_y^2 + p_z^2 = 1`.

  .. xmldoc:: <KEYWORD MODULE="SINGLE_ANISO" NAME="QUAX" KIND="STRING" LEVEL="BASIC">
              %%Keyword: QUAX <basic>
              <HELP>
              This keyword controls the quantization axis for the computation of the Crystal-Field parameters acting on the ground atomic multiplet of a lanthanide. On the next line, the program will read one of the three values:
              ||1 -- Zm of the ground pseudospin multiplet
              ||2 -- Zm of the ground atomic multiplet
              ||3 -- defined by the user on the following line
              </HELP>
              </KEYWORD>

:kword:`UBAR`
  This keyword allows estimation of the structuere of the blocking barier of a single-molecule magnet. The default is not to compute it.
  The method prints transition matix elements of the magnetic moment according to the :numref:`fig:ubar`.

  .. figure:: ubar.*
     :name: fig:ubar
     :width: 75%
     :align: center

     Pictorial representation of the low-lying energy structure of a single-molecule magnet. A qualitative performance picture of the investigated single-molecular magnet is estimated by the strengths of the transition matrix elements of the magnetic moment connecting states with opposite magnetizations (:math:`n{+} \to n{-}`). The height of the barrier is qualitatively estimated by the energy at which the matrix element (:math:`n{+} \to n{-}`) is large enough to induce significant tunnelling splitting at usual magnetic fields (internal) present in the magnetic crystals (0.01 -- 0.1 tesla). For the above example, the blocking barrier closes at the state (:math:`8{+} \to 8{-}`).

  All transition matrix elements of the magnetic moment are given as
  (:math:`(\vert\mu_{X}\vert+\vert\mu_{Y}\vert+\vert\mu_{Z}\vert)/3`).
  The data is given in Bohr magnetons (:math:`\mu_{\text{Bohr}}`).
  The keyword is used with no arguments.

  .. xmldoc:: <KEYWORD MODULE="SINGLE_ANISO" NAME="UBAR" KIND="STRING" LEVEL="BASIC">
              %%Keyword: UBAR <basic>
              <HELP>
              This keyword allows estimation of the structuere of the blocking barier of a single-molecule magnet. The default is not to compute it.
              The method prints transition matix elements of the magnetic moment connecting states with opposite magnetisations.
              The keyword is used with no arguments.
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

  .. xmldoc:: <KEYWORD MODULE="SINGLE_ANISO" NAME="ABCC" KIND="STRING" LEVEL="BASIC">
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

An input example
................

::

  &SINGLE_ANISO
  MLTP
  3
  4 4 2
  ZJPR
  -0.2
  ENCU
  250 400
  HINT
  0.0  20.0  100
  TINT
  0.0  330.0  331
  MAVE
  1  12
  MVEC
  3
  0.0000  0.0000   0.1000
  1.5707  0.0000   0.5000
  1.5707  1.5707   1.0000

.. xmldoc:: </MODULE>
