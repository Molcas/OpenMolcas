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
Module Cntrl
use Molcas, only: LenIn, MxAtom, MxOrb, MxRoot
use RASDim, only: mxTit
! The parameters defined in module RASDim should be private
Private mxTit
! The parameters defined in module Molcas should be private
Private LenIn, MxAtom, MxOrb, MxRoot

  INTEGER, PARAMETER :: MXJOB=100,MXPROP=30000
  INTEGER, PARAMETER :: MXDISP=500
  REAL*8  :: PNUC(MXPROP)=0.0D0,PORIG(3,MXPROP)=0.0D0,CITHR
  REAL*8  EMIN,ERFNUC,EPRTHR, EPRATHR,ALPHZ, BETAE
  REAL*8  TSTART,TINCRE,BSTART,BINCRE,BANGRES
  REAL*8  OSTHR_DIPR,OSTHR_QIPR
  REAL*8  RSTHR, TOLERANCE
  INTEGER, DIMENSION(MXPROP) :: ICOMP,ISOCMP,IPUSED,IPCODE
  INTEGER, DIMENSION(MXJOB) :: NSTAT,ISTAT,NROOTS,NACTE,MLTPLT
  INTEGER, DIMENSION(MXJOB) :: IRREP,NHOLE1,NELE3,NCONF,ISPACE
  INTEGER, DIMENSION(MXJOB) :: NDET
  INTEGER NJOB,NSTATE,NPROP,NSOPR
  INTEGER NRNATO,NBINA,IBINA(2,MxRoot)
  INTEGER LSYM1,LSYM2,NCONF1,NCONF2
  INTEGER LCI1,LCI2,LCI3,LGAM1,LGAM2,LGAM3,LIPAIR
  INTEGER NTSTEP,NBSTEP,LOOPDIVIDE,LOOPMAX
  INTEGER DYSEXPSF,DYSEXPSO
  INTEGER OCAN,DCHO
  CHARACTER(LEN=16) :: OCAA(20)
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
  LOGICAL IFCURD, Do_TMOM, Do_SK, Do_Pol, DOCD, Force_NON_AO_TDM, SaveDens, IFARGU
  REAL*8 TDIPMIN,SOTHR_PRT,TMGr_thrs
  INTEGER NSOTHR_PRT, ISMGRD(MXDISP), LDISP(8), NDISP, NTDISP(MXDISP)
  INTEGER IFJ2, IFJZ
  INTEGER L_Eff,nQuad

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
      CHARACTER(LEN=8) PNAME(MXPROP),PTYPE(MXPROP),SOPRNM(MXPROP),SOPRTP(MXPROP),RASTYP(MXJOB)
      CHARACTER(LEN=128) JBNAME(MXJOB),MINAME(MXJOB)
! JBNAME - LOGICAL NAME OF EACH JOBIPH FILE.
! PNAME  - NAME OF EACH PROPERTY FOR WHICH MATRIX ELEMENTS ARE COMPUTED
! PTYPE  - TYPE NAME, ex. 'AntiSing' for an antihermitian, spin-singlet op.
! SOPRNM - LIST OF PROPERTY NAMES, LIKE PNAME, FOR MATRIX ELEMENTS OVER
!          SPIN-ORBIT STATES.
! SOPRTP - TYPE NAME, similar to PTYPE
! RassiT - Title of the Rassi-calculation.
  LOGICAL PRDIPVEC,PRDIPCOM,PRSXY,PRORB,PRTRA
  LOGICAL PRCI,CIH5,IFHAM,IFHEXT,IFHEFF,IFEJOB,IFHCOM
  LOGICAL HAVE_HEFF,HAVE_DIAG,NOHAM
  LOGICAL IFTRD1,IFTRD2,IFTDM,HOP,TRACK,ONLY_OVERLAPS
  LOGICAL IFSHFT,IFHDIA,IFSO,IFTD2,NATO,RFpert,ToFile
  LOGICAL BINA
  LOGICAL PRXVR,PRXVE,PRXVS,PRMER,PRMEE,PRMES
  LOGICAL IFGCAL,IFXCAL,IFMCAL,DQVD
  LOGICAL DIPR,QIPR,QIALL
  LOGICAL RSPR
  LOGICAL DYSO,DYSEXPORT,TDYS,DCHS
  LOGICAL QDPT2SC, QDPT2EV
  LOGICAL PRRAW,PRWEIGHT
  LOGICAL REDUCELOOP
  LOGICAL SECOND_TIME,DoGSOR
  LOGICAL RHODyn

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
  LOGICAL IFACAL,IFACALFC,IFACALSD
  LOGICAL IFACALFCON,IFACALSDON,IFACALPSO
  LOGICAL IFACALFCSDON,IFVANVLECK,IFSONCINI
  LOGICAL IFSONCIFC,IFGTCALSA,IFGTSHSA,IFATCALSA
  INTEGER NTS,NTP,NTF,MULTIP
  REAL*8  TMINS,TMAXS,TMINP,TMAXP
! tjd- BMII: LPRPR set to .T. for easier parsable matrix output
! tjd- Yoni: LHAMI
  LOGICAL LPRPR,LHAMI
  REAL*8  TMINF,TMAXF

! BP - Testing flags
! NOSO      Disable SO contributions in the SONATORB and SODIAG code
  LOGICAL NOSO

!nf
  Logical IfDCpl
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
!      NTO Calculation Section /// Jie Bao
  Logical IfNTO

!    SONTO            Array of SO state pairs
Integer, Allocatable, Public:: SONTO(:,:)
!    SONTOSTATES      Number of state pairs to calculate
Integer, Public:: SONTOSTATES=0

!    SONAT            Array of SO state to compute
Integer, Allocatable, Public:: SONAT(:)
!    SONATNSTATE      Number of states to calculate
Integer, Public:: SONATNSTATE=0

Integer, Allocatable, Public:: SODIAG(:)
!    SODIAGNSTATE     Number of states to diagonalize
INTEGER, Public::  SODIAGNSTATE=0

Real*8, Allocatable, Public:: RefEne(:), HEff(:,:)


INTEGER, Parameter:: MORSBITS=8

! Note: MXATOM to be taken from module Molcas
INTEGER               NGROUP,IGROUP(8),NATOMS,COOR(3,MXATOM)
! Atom labels, 4 bytes each.
CHARACTER(LEN=LENIN) ATLBL(MXATOM)


! TEMPORARY DATA FROM JOBIPHS
REAL*8 POTNU1
Integer NACTE1,MPLET1,NSYM1,NFRO1(8),NISH1(8),NASH1(8),NDEL1(8),NBAS1(8),NRS11(8),NRS21(8), &
        NRS31(8),LROT1,NROOT1,IROOT1(mxRoot),NHOL11,NELE31
CHARACTER(LEN=LenIn+8) NAME(mxOrb)
CHARACTER(LEN=2) HEAD1(72)
CHARACTER(LEN=4) TITLE1(18,mxTit)

! Cholesky/RI stuff
Integer ALGO,Nscreen
Real*8 dmpk


! This preserves the values of the variables for the
! trajectory surface hopping algorithm.
INTEGER  ISTATE1,ISTATE2,NCI1,NCI2,nHop
LOGICAL  ChkHop,lHop
! ISTATE1 - Number of the current relaxed state
! ISTATE2 - Number of the state that interacts with the current state
! NCI1    - Configuration's number of state 1
! ChkHop  - If .TRUE. switches on the TSH algorithm
! lHop    - Is .TRUE. if nHop is set in the RunFile
! nHop    - Number of transitions (Hops) already occured


!----------------------------------------------------------------------*
!     Define files ( file names and unit numbers )                     *
!----------------------------------------------------------------------*
!
! LUIPH  - UNIT NUMBER OF JOBIPHS
! LUMCK  - UNIT NUMBER OF MCKINT FILES
! LUONE  - D:O, ONE-ELECTRON INTEGRAL FILE
! LUORD  - D:O, ORDERED TWO-ELECTRON INTEGRAL FILE
! IADR15 - TABLE OF CONTENTS, DISK ADDRESSES ON LUIPH.
! IDCMO  - Addresses to the CMO arrays on each JOBIPH
Character(LEN=8) FnOne,FnIph,FnMck,FnOrd,FnTDM,FnExc,FnToM,FnEig
INTEGER LUIPH, LUMCK, LUONE, LUORD, LUEIG
INTEGER LUEXC, LUTDM, LUTOM, IDCMO(MXJOB), ITOC15(30)

End Module Cntrl
!
