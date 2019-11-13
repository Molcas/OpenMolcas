.. index::
   single: Program; NEMO
   single: NEMO

.. _UG\:sec\:nemo:

:program:`nemo` |extramark|
===========================

.. warning::

   This program is not available in |openmolcas|

.. only:: html

  .. contents::
     :local:
     :backlinks: none

.. _UG\:sec\:nemo_description:

Description
-----------

.. xmldoc:: %%Description:
            The Nemo program of the molcas program system generates:
            fitting of potential surfaces, energy optimizations, potential curves and simulation parameters.

The :program:`NEMO` is a potential analyzis package that calculates interaction energies between molecules. The package uses input files from :program:`MPPROP` and :program:`MKNEMOPOT`.
The package was originally a set of programs that has been totaly rewritten and put together into one program. The package
are capable of doing fitting of potential surfaces, energy optimization between molecules, calculate some specific potential curves and generate simulation parameters for rigid molecules.

The theoretical background stands in perturbation theory. The interaction
energy between two molecules can be described by three quantum chemical
calculations. One quantum chemical calculation for each of the monomers
and one calculation for the two molecules together i.e. the dimer. The
energy for the two monomers are then subtracted from the dimer
calculation. That is done for each configuration, i.e. coordinate set
for a dimer calculation, given by the input to :program:`MKNEMOPOT`.
The calculations are set up by the :program:`MKNEMOPOT` package and to
performe those calculations it is recommended to read the manual for
the :program:`MKNEMOPOT` package. The interaction energy can also be
described in the classical energy terms electrostatic, induction dispersion,
repulsion and chargetransfer. Where, a good description the first three
energy terms can be given by distributed multipole expansions and distributed
polarizabilities. The last two energy terms are harder to predict and
are of quantum chemical origin. The reason for calculating the interaction
energy quantum chemically is that this reference energy will be used for
the description of the repulsion(/chargetransfer) parameters. The repulsive
reference energy term is achieved by subtracting the energy for the electrostatic,
induction and dispersion from the reference energy. Note here that the
dispersion energy is only added if the reference energy is performed
with a method that includes true dynamic correlation, i.e. when the energy
includes the London dispersion term. The reference energy will also include
a charge transfer term if it is defined by the user. An estimation of the
repulsive energy term can now be fitted to the reference repulsive energy
term by using the FITPar subprogram in the :program:`NEMO` program. The
fitted parameters are classified in elements and type. Where a hydrogen
atom is element 1 and can be classified in different types depending on
their chemical environment. This information is supplied with the MpProp
file together with coordinates, multipole and polarizabilities of a
molecule. The MpProp file is an output from the :program:`MPPROP` or
:program:`LOPROP` program. A MpProp file does not always contain all the
information needed to run the :program:`NEMO` program. Thus, it is important
have look directly in the file and do your prefered changes before using
it. It can for example be to change the type of a hydrogen atom. If we
take the ethanol molecule as an example. It is composed of two carbons,
five hydrogens and one oxygen atoms. Here we can define three different
type of hydrogens that are bonded to C1, C2 and O1 respectively. The two
carbons in the molecule can of course also be defined to be of different
type. For each defined type there excists two corresponding parameters for
the repulsion energy. These are the ones that are varied to in the fitting
procedure.

The input file comming from the :program:`MKNEMOPOT` program can contain
a cluster definition. A cluster is defined as a supermolecule containing
one/several different/equal molecules. The interaction energy is thus
defined as the interaction between different clusters.

The POTSurf subprogram produces potential energy curves between two
clusters. This is normally used to compare the fitted potential with the
result from a quantum chemical calculation. Whats happening is that one
of clusters are translated and rotated to a certain position. The moved
cluster is then translated along a displacement vector.

In the DIMEr subprogram an optimizition/minimizition of the energy between
two/several molecules is performed. The routine is not good and practical
for many molecules. Because, it was originally written to do the job for
two molecules which works pretty good.

The SIMPar program can produce input files for the :program:`MOLSIM` package.

.. _UG\:sec\:nemo_dependencies:

Dependencies
------------

The :program:`NEMO` program requires a nemo library.
The library is just a concatenation of several different :file:`nemo` files.
In order to run the FitPar subprogram in :program:`NEMO` a :file:`NEMO` file is required.
The :file:`NEMO` file is either autogenerated through the :program:`MKNEMOPOT` or it might be
generated by hand from some other potential.

.. index::
   pair: Files; NEMO

.. _UG\:sec\:nemo_files:

Files
-----

Below is a list of the files that are used/created by the program
:program:`NEMO`.

Input files
...........

.. class:: filelist

:file:`NEMO`
  This file will be opened in the $WorkDir/ directory and is composed of several :file:`Nemo` files
  generated by :program:`MKNEMOPOT`.

:file:`ATOMPAR`
  This file will be opened in the $WorkDir/ directory and it holds the atomic parameters for repulsion, scaling constants for the
  dispersion, valence of the atoms. It will originaly be stored in the $MOLCAS/nemo_libary directory. It's definition is:
  two dummy lines, nElements=103 of lines and all this taken nType=4 times. The signifacant nElements of lines will hold 12 columns.
  Where the first column is the element number, the second column is the element label, the third column

  Columns in the ATOMPAR file:

  * **Column=1**
    element number

  * **Column=2**
    element label

  * **Column=3**
    Alpha

  * **Column=4**
    Kappa

  * **Column=5**
    Charge Transfer Alpha

  * **Column=6**
    Charge Transfer Kappa

  * **Column=7**
    Valence of the atom

  * **Column=8**
    RepExp an integer for the :math:`r^{-n}` type potential.

  * **Column=9**
    RepFac

  * **Column=10**
    DispFac

  * **Column=11**
    K1/Sigma

  * **Column=12**
    K2/Epsilon

Output files
............

.. class:: filelist

:file:`POTSURF`
  This file holds the potential curve. The columns of the PotSurf file will be:

  * **Column=1**
    Coordinate 1

  * **Column=2**
    Electrostatic+Induction+Repulsion

  * **Column=3**
    Electrostatic+Induction+Repulsion+Dispersion

  * **Column=4**
    Electrostatic

  * **Column=5**
    Induction

  * **Column=6**
    Dispersion

  * **Column=7**
    Repulsion

  * **Column=8**
    Charge Transfer

:file:`MOLSIM`
  The input file in molsim format for the particle part.

:file:`MOLSIMLIB`
  The library file in molsim format for the repulsive and dispersive part.

:file:`ATOMFIT`
  This is the same file as ATOMPAR, but it is written to the $WorkDir directory

.. index::
   pair: Input; NEMO

.. _UG\:sec\:nemo_input:

Input
-----

.. compound::

  Below follows a description of the input to :program:`NEMO`. The keywords
  are always significant to four characters, but in order to make the
  input more transparent, it is recommended to use the full keywords.
  The :program:`NEMO` program section of the |molcas| input is bracketed by
  a preceding dummy namelist reference ::

    &NEMO &END

  and an "end of input" statement ::

    End of Input

Argument(s) to a keyword are always supplied on the next line of the
input file, except explicitly stated otherwise.

Optional general keywords
.........................

.. class:: keywordlist

:kword:`ALPHa`
  Use this Keyword to define the alpha parameter for a specific atom and atomtype.
  The keyword should be followed by a line/lines composed of the element number,
  the atomtype and the value for alpha.This Keyword should be ended by a END statement
  in the last line. The example below means that uran type 1 will have the value 0.1 .
  The alpha parameter will be used in the exponent for the repulsion. ::

    ALPHa
    92 1 0.1
    END

  .. xmldoc:: %%Keyword: ALPHa <basic>
              Use this Keyword to define the alpha parameter for a specific atom and atomtype.
              The should be followed by a line/lines composed of the element number,
              the atomtype and the value for alpha.This Keyword should be ended by a END statement
              in the last line. The example below means that uran type 1 will have the value 0.1 .

:kword:`KAPPa`
  Use this Keyword to define the kappa parameter for a specific atom and atomtype.
  The keyword should be followed by a line/lines composed of the element number,
  the atomtype and the value for kappa.This Keyword should be ended by a END statement
  in the last line. The example below means that uran type 1 will have the value 10.0 .
  The kappa parameter will be used as a prefactor to the exponent expression for the repulsion. ::

    KAPPa
    92 1 10.0
    END

  .. xmldoc:: %%Keyword: KAPPa <basic>
              Use this Keyword to define the kappa parameter for a specific atom and atomtype.
              The keyword should be followed by a line/lines composed of the element number,
              the atomtype and the value for kappa.This Keyword should be ended by a END statement
              in the last line. The example below means that uran type 1 will have the value 10.0 .
              The kappa parameter will be used as a prefactor to the exponent expression for the repulsion.

:kword:`ALCT`
  This keyword is for the charge transfer term that can be used if one specifies that in the NEMO keyword.
  The energy term is exactly the same expression as the repulsion, but with a minus sign instead.
  Use this Keyword to define the charge transfer alpha parameter for a specific atom and atomtype.
  The keyword should be followed by a line/lines composed of the element number,
  the atomtype and the value for charge transfer alpha.This Keyword should be ended by a END statement
  in the last line. The example below means that uran type 1 will have the value 0.1 .
  The charge transfer alpha parameter will be used in the exponent for the repulsion. ::

    ALCT
    92 1 0.1
    END

  .. xmldoc:: %%Keyword: ALCT <basic>
              This keyword is for the charge transfer term that can be used if one specifies that in the NEMO keyword.
              The energy term is exactly the same expression as the repulsion, but with a minus instead.
              Use this Keyword to define the charge transfer alpha parameter for a specific atom and atomtype.
              The keyword should be followed by a line/lines composed of the element number,
              the atomtype and the value for charge transfer alpha.This Keyword should be ended by a END statement
              in the last line. The example below means that uran type 1 will have the value 0.1 .
              The charge transfer alpha parameter will be used in the exponent for the repulsion.

:kword:`KACT`
  This keyword is for the charge transfer term that can be used if one specifies that in the NEMO keyword.
  The energy term is exactly the same expression as the repulsion, but with a minus sign instead.
  Use this Keyword to define the charge transfer kappa parameter for a specific atom and atomtype.
  The keyword should be followed by a line/lines composed of the element number,
  the atomtype and the value for charge transfer kappa. This Keyword should be ended by a END statement
  in the last line. The example below means that uran type 1 will have the value 10.0 .
  The charge transfer kappa parameter will be used as a prefactor to the exponent expression for the repulsion. ::

    KACT
    92 1 10.0
    END

  .. xmldoc:: %%Keyword: KACT <basic>
              This keyword is for the charge transfer term that can be used if one specifies that in the NEMO keyword.
              The energy term is exactly the same expression as the repulsion, but with a minus sign instead.
              Use this Keyword to define the charge transfer kappa parameter for a specific atom and atomtype.
              The keyword should be followed by a line/lines composed of the element number,
              the atomtype and the value for charge transfer kappa. This Keyword should be ended by a END statement
              in the last line. The example below means that uran type 1 will have the value 10.0 .
              The charge transfer kappa parameter will be used as a prefactor to the exponent expression for the repulsion.

:kword:`REPFactor`
  If a repulsion of type :math:`\sqrt{F_1 F_2}r^{-n}` is to be used.
  Check the NEMO keyword for information. This keyword is specified in the same way as kappa.

:kword:`DISPfactor`
  Two factors are multiplied with the dispersion energy. They work in the same way as the REPFactor does and
  are specified in the same way.

:kword:`VALEnce`
  Set the number of valence electrons. The keyword should be followed by a line/lines composed of the element number,
  the atomtype and the value for kappa.This Keyword should be ended by a END statement
  in the last line. The example below means that oxygen type 2 will have 6 valence electrons. ::

    VALEnce
    8 2 6.0
    END

  .. xmldoc:: %%Keyword: VALEnce <basic>
              Set the number of valence electrons. The keyword should be followed by a line/lines composed of the element number,
              the atomtype and the value for kappa.This Keyword should be ended by a END statement
              in the last line. The example below means that oxygen type 2 will have 6 valence electrons.

:kword:`NOISotropicPolarizabilities`
  The default is to use isotropic polarizabilities for the induction energy.
  This is due to the fact that we use Thole damping as default, which require isotropic
  polarizabilities.

  .. xmldoc:: %%Keyword: NOISotropicPolarizabilities <basic>
              The default is to use isotropic polarizabilities for the induction energy.
              This is due to the fact that we use Thole damping as default, which require isotropic
              polarizabilities.

:kword:`NOMOve`
  The default interactions sites are not placed in the atoms. If this keyword is used
  the interactions sites are not moved to a new location.

  .. xmldoc:: %%Keyword: NOMove <basic>
              Do not move the interactions sites which is the default.

:kword:`NOQUadrupoleDelete`
  The default is to replace the quadrupoles with local dipoles to get the correct total quadrupole.
  If this keyword is used, the quadrupoles will be truncated at the dipole level.

  .. xmldoc:: %%Keyword: NOQUadrupoleDelete <basic>
              The default is to replace the quadrupoles with local dipoles to get the correct total quadrupole.
              If this keyword is used, the quadrupoles will be truncated at the dipole level.

:kword:`NODAmping`
  As default the Thole damping is used, but using this heyword that is overruled.

  .. xmldoc:: %%Keyword: NODAmping <basic>
              As default the Thole damping is used, but using this heyword that is overruled.

:kword:`REPLace`
  Use this keyword to specify that some atomic quadrupoles should be replaced by charges.

:kword:`MOLD`
  The new local atomic dipole will be used when calculating the new interaction center.
  The default is to use the original local atomic dipole.

:kword:`NOLM`
  The new local atomic quadrupole will be used when estamating the size of the atom.
  This is used when calculating the repulsion and dispersive energy.
  The default is to use the original local atomic quadrupole which is the correct way.

:kword:`RETY`
  REpTYpe: The keyword should be followed by a line, specifying the expression to use for the repulsion type.

  Optional RETY parameters:

  * **m=0**
    (Default) Here the exponent is described by :math:`-r_{12}/(\sqrt{\Tr(Q_1)/3/qv_1+\Tr(Q_2)/3/qv_2}(\alpha_1\alpha_2))`.

  * **m=1**
    Here the exponent is described by :math:`-r_{12}/(\alpha_1\sqrt{\Tr(Q_1)/3/qv_1}+\alpha_2\sqrt{\Tr(Q_2)/3/qv_2})`.

:kword:`NEMO`
  The keyword should be followed by a line, what kind of energy expression to use.
  The parameters for the energies are read from the :file:`nemo` and :file:`ATOMPAR`

  Optional NEMO parameters:

  * **m=0**
    (Default) Electrostatic, inductive, dispersive and a exponetial repulsion energy term is used.

  * **m=1**
    Here a :math:`\sqrt{F_1 F_2}r^{-n}` type repulsion is added to the default energy.

  * **m=2**
    Here dispersion factors are used to scale the energy.

  * **m=3**
    This number means that default energy is used, plus the repulsive term of type 1 and the dispersive scaling of type 2.

  * **m=4**
    A charge transfer term is added to the default energy, which has the same expression
    as the repulsion term only differing in the sign.

  .. :kword:`AMBEr`
       Not implemented. For future use.

  ..   .. xmldoc:: %%Keyword: AMBEr <basic>
                   Not implemented. For future use.

  .. :kword:`SIGMa`
       Not implemented. For future use.

  ..   .. xmldoc:: %%Keyword: SIGMa <basic>
                   Not implemented. For future use.

  .. :kword:`EPSIlon`
       Not implemented. For future use.

  ..   .. xmldoc:: %%Keyword: EPSIlon <basic>
                   Not implemented. For future use.

:kword:`SEED`
  The seed to the random generator.

  .. xmldoc:: %%Keyword: SEED <basic>
              The seed to the random generator.

:kword:`FITPar`
  This is the start keyword for the subprogram :program:`FITPAR`. It should consist of the Keyword plus a END statement.
  Inbetween there should be :program:`FITPAR` specific keywords.
  The subprogram to do the fitting of parameters.

:kword:`DIMEr`
  This is the start keyword for the subprogram :program:`DIMER`. It should consist of the Keyword plus a END statement.
  Inbetween there should be :program:`DIMEr` specific keywords.
  The subprogram do an energy minimisation for two monomers.

:kword:`POTSurf`
  This is the start keyword for the subprogram :program:`POTSURF`. It should consist of the Keyword plus a END statement.
  Inbetween there should be :program:`POTSURF` specific keywords.
  The subprogram generates potential curves.

:kword:`SIMPar`
  This is the start keyword for the subprogram :program:`SIMPAR`. It should consist of the Keyword plus a END statement.
  Inbetween there should be :program:`SIMPAR` specific keywords.

Optional FITPar specific keywords
.................................

These keywords should begin by a FITPar keyword and end with a END statement.

.. class:: keywordlist

:kword:`NUAL`
  NO UPDATE ALPHA. This keyword should be followed by a line/lines specifying the element and type
  of the atomic parameter that should not be updated during the fitting. The example says that the
  oxygen type 2 atomic parameter should not be updated. ::

    NUAL
    8 2
    END

:kword:`NUKA`
  NO UPDATE KAPPA. This keyword should be followed by a line/lines specifying the element and type
  of the atomic parameter that should not be updated during the fitting. The example says that the
  oxygen type 2 atomic parameter should not be updated. ::

    NUKA
    8 2
    END

:kword:`NUAC`
  NO UPDATE CHARGE TRANSFER ALPHA. This keyword should be followed by a line/lines specifying the element and type
  of the atomic parameter that should not be updated during the fitting. The example says that the
  oxygen type 2 atomic parameter should not be updated.This only works for NEMO type 4. Check the NEMO keyword. ::

    NUAC
    8 2
    END

:kword:`NUKC`
  NO UPDATE CHARGE TRANSFER KAPPA. This keyword should be followed by a line/lines specifying the element and type
  of the atomic parameter that should not be updated during the fitting. The example says that the
  oxygen type 2 atomic parameter should not be updated.This only works for NEMO type 4. Check the NEMO keyword. ::

    NUKC
    8 2
    END

:kword:`NUSI`
  Not implemented. For future use.

:kword:`NUEP`
  Not implemented. For future use.

:kword:`NURE`
  NO UPDATE REPULSION FACTOR. This keyword should be followed by a line/lines specifying the element and type
  of the atomic parameter that should not be updated during the fitting. The example says that the
  oxygen type 2 atomic parameter should not be updated. This only works for NEMO type 1 and 3. Check the NEMO keyword. ::

    NUKC
    8 2
    END

:kword:`NUDI`
  NO UPDATE DISPERSION FACTOR. This keyword should be followed by a line/lines specifying the element and type
  of the atomic parameter that should not be updated during the fitting. The example says that the
  oxygen type 2 atomic parameter should not be updated. This only works for NEMO type 2 and 3. Check the NEMO keyword. ::

    NUKC
    8 2
    END

:kword:`GLOBal`
  The keyword should be followed by a line specifying the number of globalsteps.

:kword:`MACRo`
  The keyword should be followed by a line specifying the number of macrosteps.

:kword:`MICRo`
  The keyword should be followed by a line specifying the number of microsteps.

:kword:`TEMP`
  The keyword should be followed by a line specifying the temperature for the weighting procedure. See the keyword WEIG.

:kword:`SCFFit`
  By default the program tries to fit the second energy term in the NEMO file. Using this keyword the program uses the
  first energy term witch is a SCF type energy.

:kword:`CONVergence`
  The keyword should be followed by a line specifying the number for the convergence radii.

:kword:`RFACtor`
  The keyword should be followed by a line specifying the number for the scaling constant in the least square fit.

:kword:`WEIGht`
  The keyword should be followed by a line specifying the number of the weight type

  Optional WEIGht parameters:

  * **m=0**
    (Default) Weight=Min(2,Exp( -0.2*(E(dimer)-E(Monomer1)-E(Monomer2)) )

  * **m=1**
    Weight=exp(-(E(dimer)-E(Monomer1)-E(Monomer2))/kT)

:kword:`ERROr`
  The keyword should be followed by a line specifying the number of the error type

  Optional ERROr parameters:

  * **m=0**
    (Default) Error=Weight*( Exp( 0.15D0*(E(estimated)-E(reference)) )-1 )**2

  * **m=1**
    Error=Weight*(E(reference)-E(estimated))**2

:kword:`DISFactor`
  The keyword should be followed by a line specifying a scaling constant for the dispersion energy. (Default 1.0)

:kword:`LINEarsearch`
  The keyword can contain any of the keywords FORCe, SIMPlex, ITERation and CONVergence. It should also finnish by an END statement.

  .. Optional LINEar specific keywords:

:kword:`SIMPlex`
  Keyword for the simplex method.

:kword:`FORCe`
  Keyword for a steepest descent type method.

:kword:`ITERation`
  The keyword should be followed by a line specifying the number of interations.

:kword:`CONVergence`
  The keyword should be followed by a line specifying the number for the convergence.

Optional DIMEr specific keywords
................................

These keywords should begin by a DIMEr keyword and end with a END statement.

.. class:: keywordlist

:kword:`MOLEcules`
  The keyword should be followed by a line specifying a molecule by name exactly as they are named in the nemo file. All other molecular based keywords will be given to this molecule. That until a new molecule name is given with this keyword.

:kword:`METHod`
  Specifies the method to be used for the file to be opened. The program will find another method if the specified method cannot be found in the MPPROP file.

:kword:`MACRosteps`
  The keyword should be followed by a line specifying the number of macrosteps.

:kword:`MICRosteps`
  The keyword should be followed by a line specifying the number of microsteps.

:kword:`STARt`
  The keyword should be followed by a line specifying two numbers. The first number is search radii for coordinates and
  the second number is the search radii for the angles. In the first macrostep.

:kword:`RFACtor`
  The keyword should be followed by a line specifying the number of the scaling factor for the search radii each macrostep.

:kword:`CONVergence`
  The keyword should be followed by a line specifying the number for the convergence radii.

:kword:`DISFac`
  The keyword should be followed by a line specifying a scaling constant for the dispersion energy. (Default 1.0)

Optional POTSurf specific keywords
..................................

These keywords should begin by a POTSur keyword and end with a END statement.

.. class:: keywordlist

:kword:`MOLEcule`
  Specifies the start and the title of a new molecule. This means every keyword after this MOLEcule keyword will belong to the last specified MOLEcule.

:kword:`METHod`
  Specifies the method to be used for the file to be opened. The program will find another method if the specified method cannot be found in the MPPROP file.

:kword:`CLUSter`
  This keyword should be followed by a line that gives an integer number of witch cluster the lates molecule belongs to. Only the integer numbers 1 and 2 are valid for the PotSurf module.

:kword:`TROR`
  This keyword should be followed by a line that gives six numbers. The six numbers describes the translation in polar coordinates and the rotation in the three euler angles for the molecule given by the latest MOLEcule keyword. The sequence of the numbers are the following: R Theta Phi Alpha Beta Gamma (See Arfken for definitions)

:kword:`POTEntial`
  The keyword should be followed by one line specifying three numbers. The numbers gives the displacement vector in spherical poolar coordinates for the second cluster when calculating the potential energy. The numbers are given in the following order: R Theta Phi

:kword:`NPOInts`
  The keyword should be followed by one line specifying the number of points in the potential.

:kword:`TRANslation`
  The keyword should be followed by a line specifying up to five numbers. The first number specifies the type of potential coordinates.
  In order to visulize the potential curve one has to define a translation coordinate.
  The first column of the PotSurf file will consist of a coordinete specified by the iTrType parameter. The other parameters jTrType, kTrType ... are specified below.

  Optional TRANslation parameters:

  * **iTrType=0**
    The coordinate will be the length of the translation vector. (Default)

  * **iTrType=1**
    jTrType=coordinte (1=X,2=Y and 3=Z) index of kTrType=molecule given by the order of the apperence in the input section.

  * **iTrType=2**
    jTrType=Atom1 and kTrType=Atom2 on molecule=lTrType and mTrType respectively. The molecules are given by the order of the apperence in the input section. The potential coordinate will be the distance between Atom1 and Atom2. Note that if the potential coordinate is constant if the molcules belong to the same cluster.

:kword:`DISFactor`
  The keyword should be followed by a line specifying a scaling constant for the dispersion energy. (Default 1.0)

Optional SIMPar specific keywords
.................................

These keywords should begin by a DIMEr keyword and end with a END statement.

.. class:: keywordlist

:kword:`MOLEcules`
  Specifies the start and the title of a new molecule. This means that every keyword after this
  MOLEcule keyword will belong to the last specified MOLEcule.

:kword:`METHod`
  Specifies the method to be used for the file to be opened.
  The program will find another method if the specified method
  cannot be found in the MPPROP file.

:kword:`MOLSim`
  Tells the program to generate Molsim parameters and input files.

:kword:`EQUAlatoms`
  This keyword should be followed by a line specifying two atom numbers that should treated as equal.
  The atomic numbers are the numbers in sequence as they are found in the MPPROP file.
  For example, a water molecule in gasphase has the two hydrogen atoms equal by symmetry.
  They should thus be treated equally for the analysis in a simulation program. If the
  MPPROP file has the atoms in the sequence O H H the example below makes the two hydrogen equal ::

    EQUA
    2 3

:kword:`NUMBer`
  The keyword should by a line giving the number of latest molecule that will
  be used in the latter simulation. This information will be written in the MOLSIM file.

:kword:`DISFactor`
  The keyword should be followed by a line specifying a scaling constant for the dispersion energy. (Default 1.0)

Limitations
...........

The program package has no internal degrees of freedom.
The program cannot handle interactions including quadrupoles and higher.
The program cannot handle hyperpolarizabilities. For the time being we cannot handle more than two clusters.

.. Contacts
   ........

   It is hard to see what should be included in this manual, but if you have any questions or problems just send an email to daniel.hagberg@teokem.lu.se .
