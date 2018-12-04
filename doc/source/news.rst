New features and updates
========================

Below is presented a list of the major new features of |molcas|.
These features comprise a number of new codes and
introduction of new methods, but also considerable updates of many of the
programs in |molcas|. We keep some history, so that people who are using older
versions of |molcas| can get a feeling for what has happened on later versions

.. rubric:: New features in 8.0:

* General improvements:

  * includes major bug fixes;
  * enhanced performance;
  * better parallelization;
  * better support for the Intel, and GCC compilers;

* New codes and major updates:

  * MC-PDFT combines multiconfigurational wavefunctions with density functional theory to recover both static and dynamical correlation energy;
  * GASSCF allows for more flexibility in choosing the active space;
  * EMBQ is general purpose embedding technique;
  * FALCON is fragment-based approach for computing an electronic energy of the large systems;
  * GEO/HYPER module for constrained multi-fragment geometry optimisation in internal coordinates;
  * enhanced I/O via the Files In Memory technology;
  * SINGLE_ANISO code received several important updates:

    * CRYS: extraction of the parameters of the multiplet-specific crystal field for lanthanides;
    * UBAR: construction of the blocking barriers of single-molecule magnets;
    * ABCC: magnetic and anisotropy axes are given in the crystallographic :math:`abc` system;

* New features in existing codes:

  * Relativistic exact decoupling (X2C/BSS/infinite-order DKH);
  * Local X2C/BSS/infinite-order DKH;
  * RICD analytical gradients are available for the MBPT2, and CASSCF methods;
  * auto-segmentation in CD-based coupled cluster CHCC and CHT3 modules;
  * Orbital-free density embedding;
  * more robust and efficient SLAPAF module;
  * enhanced EMIL functional;

* Installation and tools:

  * first release of the Global Arrays free MOLCAS; a new parallel framework of MOLCAS requires only MPI-2 library;
  * better support for Mac OS X (including the both serial and parallel installations);

..
  .. rubric:: New features in 7.6:

  * Bug fixing release
  * Short guide for Molcas
  * GUI-ready release

  .. rubric:: New features in 7.4:

  * New codes and major updates:

    * There is a new set of coupled cluster codes added.
    * The M06 DFT functional have been implemented.
    * There are added constraints in :program:`slapaf`.
    * New method for transition state search and reaction coordinate analysis.

  * New features in existing codes:

    * There are improvements in the capabilities of the emil input.
    * It is now possible to specify the actual name of the orbital
      input files in modules :program:`SCF` and :program:`GRID_IT`.

  * Changes in usage of the package:

    * You can now get properties broken down by orbital contributions
      by setting environment variable.

  * Installation and tools

    * You can now tell |molcas| at configuration time to use an externally
      installed version of Global Arrays.
    * There are prebuilt versions of the GUI that can be installed in a very
      simple manner.
    * The default compiler on linux system is now gfortran.

  .. rubric:: New features in 7.2:

  * New codes and major updates:

    * pre-release version of GUI for input generation and |molcas| job submition (MING).
    * Module Seward has been split into Gateway (set up of molecular system)
      and Seward itself (computation of integrals).
    * Major improvements in runtime settings for the package, and new flags for |molcas| command
    * New manual for novice |molcas| users

  * Performance enhancements:

    * A new version of GA has been included.
    * Default integral thresholds are now changed to 1.0D-10.
    * RI code has been improved

  * New features in existing codes:

    * The exchange-hole dipole moments in :program:`LoProp` code
    * Better handling of supersymmetry in :program:`RASSCF` code
    * Localized natural orbitals in :program:`Localisation` code
    * BSSE calculations in :program:`SCF` code
    * A second finite nuclei charge distribution model, the so-called modified Gaussian charge distribution,
      has been implemented
    * Frequency calculations for :program:`MBPT2`
    * The :program:`ESPF` module can be used in order to compute electrostatic potential derived charges
    * Frozen Natural Orbital approach in :program:`CASPT2`
    * On-the-fly generation of RI auxiliary basis set
    * Flexible selection of orbitals in :program:`GRID_IT`
    * New features in GV code: visualization of molden files, selection of atomic groups, symmetry operations

  * Changes in usage of the package:

    * No shell scripts are needed to run |molcas|.
    * New EMIL commands for file handling
    * Control of the print level of the code

  * Installation and tools

    * New tools for memory and I/O profiling
    * New configuration files has been included

  .. rubric:: New features in 7.0:

  * New codes and major updates:

    * CHOLESKY --- a new approach to ab initio and first principle QM methods free
      from explicit two-electron integrals. SCF/DFT, RASSCF, RASSI and MP2 energy
      calculation can now be done with considerable improvement of performance
      and with controlled accuracy of the results.
    * The 1-center approximation of the Cholesky decomposition, 1-CCD
    * Resolution of Identity (RI)/ Density fitting (DF) scheme for SCF, DFT,
      CASSCF, RASSI and CASPT2
    * The :program:`CASPT2` module can be used in connection with Cholesky and RI/DF approximations,
      allowing for the treatment of larger systems
    * Update of guessorb code
    * Electrostatic potential fitted (ESPF) QM/MM interface for SCF, DFT,
      CASSCF, CASPT2, and CC. ESPF analytic gradients for SCF, DFT, and CASSCF.
    * Gradients for "pure" DFT for the 1-CCD, and RI/DF approximations
    * Scaled Opposite-Spin (SOS) and Scaled Spin Component (SCS) MP2 are implemented when
      using Cholesky or RI/DF approximation.
    * NEMO program: fitting of potential surfaces, energy optimizations, potential curves
      and simulation parameters.
    * interface to MOLSIM code
    * Major update for GUI code :program:`GV`, with a possibility to edit coordinates and
      visually select active space for RASSCF calculations.
    * A new program, :program:`EXPBAS`, has been introduced that allows expanding an
      orbital file from a small to a larger basis set.
    * Several different procedures for constructing localized orbitals have been
      implemented. Among them is one based on a Cholesky decomposition of the density
      matrix.

  * Performance enhancements:

    * Use of external blas libraries: lapack, GotoBLAS, Atlas, Intel MKL, ACML
    * New version of GA has been included.
    * Improved diagonalization routines and improved convergence in scf and rasscf
    * Some size limits in :program:`RASSCF` and :program:`CASPT2` have been increased or eliminated.
    * Automatic generation of starting orbitals for arbitrary valence and
      ECP basis sets.

  * New features in existing codes:

    * Natural orbitals for UHF calculations. Can, for example be used as
      starting orbitals for :program:`RASSCF`.
    * Natural Bond Order (NBO) based on the LoProp partitioning.
    * Arbitrary order Douglas--Kroll--Hess (DKH) transformation to include
      scalar relativistic effects.
    * Picture-change-corrected electric potential, electric field, and
      electric field gradient properties.
    * Automatic generation of rydberg orbitals in genano.
    * RASSI can compute g-tensors.
    * CASPT2 is able to run with Cholesky vectors instead of integrals.
    * Transverse constraint for geometry optimizations.
    * Numerical gradients for several methods.
    * Numerical IR intensities for Numerical Hessian.
    * Computation of charge capacitances for bonds using Loprop.
    * Localized exchange-hole dipole moments in Loprop.
    * Possibility to use loprop with user-defined densities.
    * Evaluation of transition density between two states.
    * Mulliken type multicenter multipole expansion and localized
      polarizablilites based on the uncoupled HF approach.
    * Several improvements and enhancements in the visualization program GV.
    * The ANO-RCC basis set is now complete covering all atoms H--Cm.
    * The GUESSORB facility is now included in :program:`SEWARD`, which automatically
      produces starting orbitals for arbitrary basis sets.

  * Changes in usage of the package:

    * Improvements in |molcas| input language.
    * |molcas| job can be submitted without shell scripts.
    * The programs are making extensive use of the runfile to simplify
      the input and eliminate unnecessary inputs.
    * automatic saving of output files (molden files, and orbital files)
    * The starting orbitals for :program:`RASSCF` can be taken from a number of sources
      (Guessorb, runfile, etc.), and this is done in a semi-intelligent
      way unless specified in user input.
    * simplified RASSCF input: number of
      orbitals, spin, etc can sometimes be deduced by the program from
      information available on the runfile or an orbital file.
      One can use CHARGE instead of the number of active electrons.
    * If used in multiple runs in one job, the RASSCF automatically
      selects suitable individual names for the JOBIPH files. The choice
      can be overridden by keyword input, but if not, it matches the
      default selection of JOBIPH names in :program:`RASSI`.
    * RASSI can use default selection of JOBIPH names, when used together with multiple
      RASSCF runs in one job.
    * RASSCF can use natural orbitals from a preceeding UHF calculation as input
      orbitals.

  * Installation and tools

    * improved installation procedure, with possibility to select compilers,
      BLAS libraries, and parallel environment.
    * Configuration files for new compilers, including gfortran, g95, SunStudio
    * Configuration files for OpenMP parallelization.
    * Tools for extracting information from RUNFILE and JOBIPH files.
