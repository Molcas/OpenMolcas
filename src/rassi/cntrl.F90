!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

module Cntrl

use Molcas, only: LenIn, MxAtom, MxOrb, MxRoot
use RASDim, only: mxTit
use Constants, only: Zero

! The parameters defined in module RASDim should be private
private mxTit
! The parameters defined in module Molcas should be private
private LenIn, MxAtom, MxOrb, MxRoot

integer, parameter :: MXJOB = 100, MXPROP = 30000
integer, parameter :: MXDISP = 500
real*8 :: PNUC(MXPROP) = Zero, PORIG(3,MXPROP) = Zero, CITHR
real*8 EMIN, ERFNUC, EPRTHR, EPRATHR, ALPHZ, BETAE
real*8 TSTART, TINCRE, BSTART, BINCRE, BANGRES
real*8 OSTHR_DIPR, OSTHR_QIPR
real*8 RSTHR, TOLERANCE
integer, dimension(MXPROP) :: ICOMP, ISOCMP, IPUSED, IPCODE
integer, dimension(MXJOB) :: NSTAT, ISTAT, NROOTS, NACTE, MLTPLT
integer, dimension(MXJOB) :: IRREP, NHOLE1, NELE3, NCONF, ISPACE
integer, dimension(MXJOB) :: NDET
integer NJOB, NSTATE, NPROP, NSOPR
integer NRNATO, NBINA, IBINA(2,MxRoot)
integer LSYM1, LSYM2, NCONF1, NCONF2
integer LCI1, LCI2, LCI3, LGAM1, LGAM2, LGAM3, LIPAIR
integer NTSTEP, NBSTEP, LOOPDIVIDE, LOOPMAX
integer DYSEXPSF, DYSEXPSO
integer OCAN, DCHO
character(len=16) :: OCAA(20)
! BP - SO Natural orbital information
! RF - SO Natural transition orbital information
! IFARGU           Do phase factor for SO-NTOs
! IFCURD           Do current density stuff
! Do_TMOM          Do (exact) Transition MOMents
! L_Eff            The order of the Lebedev solid angle integration
! Do_SK            Do k-vector directions of the exact
!                     Hamiltonian with a vector field.
! nQuad            Number of vectors in the isotropic quadrature
! Do_Pol           Specify a polarization vector direction
! DOCD             Regular circular dichroism - velocity and mixed gauge
! SaveDens         Save input-state transition densities in temp. file
logical IFCURD, Do_TMOM, Do_SK, Do_Pol, DOCD, Force_NON_AO_TDM, SaveDens, IFARGU
real*8 TDIPMIN, SOTHR_PRT, TMGr_thrs
integer NSOTHR_PRT, ISMGRD(MXDISP), LDISP(8), NDISP, NTDISP(MXDISP)
integer IFJ2, IFJZ
integer L_Eff, nQuad

! CITHR  - THRESHOLD FOR PRINTING CI COEFFICIENTS.
! ESHFT  - OPTIONAL ENERGY SHIFT OF EACH INPUT STATE.
! LROOT  - ORDERING NUMBER, ON ITS JOBIPH FILE, OF EACH INPUT STATE.
! NSTAT  - NR OF STATES TO BE PICKED UP FROM EACH JOBIPH.
! IRREP  - SYMMETRY OF THE WAVE FUNCTIONS ON EACH JOBIPH.
! NCONF  - SIZE OF CI ARRAYS ON EACH JOBIPH.
! ISPACE - Which determinant-CI space to use with each JOBIPH
! NJOB   - NR OF JOBIPH FILES TO BE USED.
! NSTATE - TOTAL NUMBER OF STATES.
! NPROP  - NR OF PROPERTIES TO COMPUTE MATRIX ELEMENTS FOR.
! NRNATO - NR OF EIGENSTATE TO COMPUTE NATURAL ORBITALS FOR.
! IPCODE - NUMERICAL CODE OF PROPERTY INTEGRALS TO USE.
!         (=ENTRY NUMBER INTO TABLE OF CONTENTS OF ONEINT FILE).
! LSYM1  - SYMMETRY OF CURRENTLY PROCESSED BRA STATE.
! LSYM2  - SYMMETRY OF CURRENTLY PROCESSED KET STATE.
! NCONF1, NCONF2, SIMILAR.
! THE REST ARE POINTERS TO DYNAMICALLY ALLOCATED ARRAYS:
! LCI1   - POINTER TO CI ARRAY OF CURRENTLY PROCESSED BRA STATE.
! LCMO1  - SIM., POINTER TO MO COEFFICIENT ARRAY.
! LTRA1  - SIM., TRNSFORMATION COEFFICIENT ARRAY.
! LCI2, LCMO2, LTRA2, AS ABOVE, BUT FOR KET STATE.
! LTUVX  - SIM., TWO-ELECTRON INTEGRALS.
! LGAM1  - POINTER TO ONE-ELECTRON TRANSITION DENSITY MATRIX.
! LGAM2  - SIM., TWO-ELECTRON MATRIX.
! LTDMAB - POINTER TO TRANSITION DENSITY MATRIX IN BION. BASIS.
! LTDMZZ - SIM., IN AO BASIS.
! iToc25 - Table-of-contents for the optional file TOFILE.
! ALPHZ - Value for alpha in DQV diabatization.
! BETAE - Value for beta in DQV diabatization.
character(len=8) PNAME(MXPROP), PTYPE(MXPROP), SOPRNM(MXPROP), SOPRTP(MXPROP), RASTYP(MXJOB)
character(len=128) JBNAME(MXJOB), MINAME(MXJOB)
! JBNAME - LOGICAL NAME OF EACH JOBIPH FILE.
! PNAME  - NAME OF EACH PROPERTY FOR WHICH MATRIX ELEMENTS ARE COMPUTED
! PTYPE  - TYPE NAME, ex. 'AntiSing' for an antihermitian, spin-singlet op.
! SOPRNM - LIST OF PROPERTY NAMES, LIKE PNAME, FOR MATRIX ELEMENTS OVER
!          SPIN-ORBIT STATES.
! SOPRTP - TYPE NAME, similar to PTYPE
! RassiT - Title of the Rassi-calculation.
logical PRDIPVEC, PRDIPCOM, PRSXY, PRORB, PRTRA
logical PRCI, CIH5, IFHAM, IFHEXT, IFHEFF, IFEJOB, IFHCOM
logical HAVE_HEFF, HAVE_DIAG, NOHAM
logical IFTRD1, IFTRD2, IFTDM, HOP, TRACK, ONLY_OVERLAPS
logical IFSHFT, IFHDIA, IFSO, IFTD2, NATO, RFpert, ToFile
logical BINA
logical PRXVR, PRXVE, PRXVS, PRMER, PRMEE, PRMES
logical IFGCAL, IFXCAL, IFMCAL, DQVD
logical DIPR, QIPR, QIALL
logical RSPR
logical DYSO, DYSEXPORT, TDYS, DCHS
logical QDPT2SC, QDPT2EV
logical PRRAW, PRWEIGHT
logical REDUCELOOP
logical SECOND_TIME, DoGSOR
logical RHODyn

! BP - Hyperfine tensor Flags
! IFACAL        TRUE to calculate hyperfine tensors
! IFACALFC      TRUE to calculate ONLY FC terms
! IFACALSD      TRUE to calculate ONLY SD terms
! K.Sharkas beg
! IFACALPSO     TRUE to calculate PSO terms individually
! IFACALFCON    TRUE to calculate FC  terms individually
! IFACALSDON    TRUE to calculate SD  terms individually
! IFACALFCSDON  TRUE to calculate FC +SD terms
! IFGTCALSA     TRUE to calculate single_aniso g-tensor in RASSI
! K.Sharkas end
logical IFACAL, IFACALFC, IFACALSD
logical IFACALFCON, IFACALSDON, IFACALPSO
logical IFACALFCSDON, IFVANVLECK, IFSONCINI
logical IFSONCIFC, IFGTCALSA, IFGTSHSA, IFATCALSA
integer NTS, NTP, NTF, MULTIP
real*8 TMINS, TMAXS, TMINP, TMAXP
! tjd- BMII: LPRPR set to .T. for easier parsable matrix output
! tjd- Yoni: LHAMI
logical LPRPR, LHAMI
real*8 TMINF, TMAXF

! BP - Testing flags
! NOSO      Disable SO contributions in the SONATORB and SODIAG code
logical NOSO

!nf
logical IfDCpl
!nf
! PRSXY  - PRINT MO OVERLAP MATRICES FOR INPUT JOBIPHS.
! PRORB  - PRINT INPUT ORBITALS.
! PRTRA  - PRINT TRANSFORMATION COEFFICIENTS.
! PRCI   - PRINT CI COEFFICIENTS OF INPUT STATES.
! CIH5   - Put CI coeffs and MO to HDF5
! IFHEXT - Spin-free Hamiltonian is taken from input.
! IFSHFT - Energy shifts of input states will be applied.
! IFHDIA - Diagonal H-matrix elements are taken from input.
! IFSO   - DO SPIN-ORBIT INTERACTION CALCULATION.
! IFTD2  - FLAG USED IN TRANS2 CALLS - CALCULATE 2-EL. TRANS.D.M.
!                              Rassi input...
! RFpert - This flag is used to signal a
!          reaction field calculation (perturbation approach).
! ToFile - Denotes if H-matrix and various one-electron matrices
!          are to be put on a file for subsequent programs.
! PRXVR, etc: Print expectation values for RasScf input states,
!          for (spin-free) eigenstates, and for SO states.
! PRMER, etc: Print matrix elements    for RasScf input states,
!          for (spin-free) eigenstates, and for SO states.
!nf
! IfDCpl - Flag for approximate derivative couplings
!nf
!IgorS 06-05-2009
! HOP    - Switch for Trajectory Surface Hopping Algorithm
! stknecht
! QDPT2SC - use SC effective Hamiltonian (rather than the PC one) from QD-NEVPT2
! QDPT2EV - use eigenvectors of effective Hamiltonian from QD-NEVPT2 to mix TDMs (in MPS-SI we do not use 'mixed MPS'
!           instead we mix the TDMs)
! NTO Calculation Section /// Jie Bao
logical IfNTO

! SONTO            Array of SO state pairs
integer, allocatable, public :: SONTO(:,:)
! SONTOSTATES      Number of state pairs to calculate
integer, public :: SONTOSTATES = 0

! SONAT            Array of SO state to compute
integer, allocatable, public :: SONAT(:)
! SONATNSTATE      Number of states to calculate
integer, public :: SONATNSTATE = 0

integer, allocatable, public :: SODIAG(:)
! SODIAGNSTATE     Number of states to diagonalize
integer, public :: SODIAGNSTATE = 0

real*8, allocatable, public :: RefEne(:), HEff(:,:)

integer, parameter :: MORSBITS = 8

! Note: MXATOM to be taken from module Molcas
integer NGROUP, IGROUP(8), NATOMS, COOR(3,MXATOM)
! Atom labels, 4 bytes each.
character(len=LenIn) ATLBL(MXATOM)

! TEMPORARY DATA FROM JOBIPHS
real*8 POTNU1
integer NACTE1, MPLET1, NSYM1, NFRO1(8), NISH1(8), NASH1(8), NDEL1(8), NBAS1(8), NRS11(8), NRS21(8), &
  NRS31(8), LROT1, NROOT1, IROOT1(mxRoot), NHOL11, NELE31
character(len=LenIn+8) NAME(mxOrb)
character(len=2) HEAD1(72)
character(len=4) TITLE1(18,mxTit)

! Cholesky/RI stuff
integer ALGO, Nscreen
real*8 dmpk

! This preserves the values of the variables for the
! trajectory surface hopping algorithm.
integer ISTATE1, ISTATE2, NCI1, NCI2, nHop
logical ChkHop, lHop
! ISTATE1 - Number of the current relaxed state
! ISTATE2 - Number of the state that interacts with the current state
! NCI1    - Configuration's number of state 1
! ChkHop  - If .TRUE. switches on the TSH algorithm
! lHop    - Is .TRUE. if nHop is set in the RunFile
! nHop    - Number of transitions (Hops) already occured

!----------------------------------------------------------------------*
!     Define files ( file names and unit numbers )                     *
!----------------------------------------------------------------------*

! LUIPH  - UNIT NUMBER OF JOBIPHS
! LUMCK  - UNIT NUMBER OF MCKINT FILES
! LUONE  - D:O, ONE-ELECTRON INTEGRAL FILE
! LUORD  - D:O, ORDERED TWO-ELECTRON INTEGRAL FILE
! IADR15 - TABLE OF CONTENTS, DISK ADDRESSES ON LUIPH.
! IDCMO  - Addresses to the CMO arrays on each JOBIPH
character(len=8) FnOne, FnIph, FnMck, FnOrd, FnTDM, FnExc, FnToM, FnEig
integer LUIPH, LUMCK, LUONE, LUORD, LUEIG
integer LUEXC, LUTDM, LUTOM, IDCMO(MXJOB), ITOC15(30)

end module Cntrl
