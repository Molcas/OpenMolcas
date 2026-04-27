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
use Definitions, only: wp, iwp

implicit none
private

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

! CITHR  - THRESHOLD FOR PRINTING CI COEFFICIENTS.
! NSTAT  - NR OF STATES TO BE PICKED UP FROM EACH JOBIPH.
! IRREP  - SYMMETRY OF THE WAVE FUNCTIONS ON EACH JOBIPH.
! NCONF  - SIZE OF CI ARRAYS ON EACH JOBIPH.
! NJOB   - NR OF JOBIPH FILES TO BE USED.
! NSTATE - TOTAL NUMBER OF STATES.
! NPROP  - NR OF PROPERTIES TO COMPUTE MATRIX ELEMENTS FOR.
! NRNATO - NR OF EIGENSTATE TO COMPUTE NATURAL ORBITALS FOR.
! LSYM1  - SYMMETRY OF CURRENTLY PROCESSED BRA STATE.
! LSYM2  - SYMMETRY OF CURRENTLY PROCESSED KET STATE.
! NCONF1, SIMILAR.
! THE REST ARE POINTERS TO DYNAMICALLY ALLOCATED ARRAYS:
! ALPHZ - Value for alpha in DQV diabatization.
! BETAE - Value for beta in DQV diabatization.

! JBNAME - LOGICAL NAME OF EACH JOBIPH FILE.
! PNAME  - NAME OF EACH PROPERTY FOR WHICH MATRIX ELEMENTS ARE COMPUTED
! PTYPE  - TYPE NAME, ex. 'AntiSing' for an antihermitian, spin-singlet op.
! SOPRNM - LIST OF PROPERTY NAMES, LIKE PNAME, FOR MATRIX ELEMENTS OVER
!          SPIN-ORBIT STATES.
! SOPRTP - TYPE NAME, similar to PTYPE

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

! tjd- BMII: LPRPR set to .T. for easier parsable matrix output
! tjd- Yoni: LHAMI

! BP - Testing flags
! NOSO      Disable SO contributions in the SONATORB and SODIAG code

! PRSXY  - PRINT MO OVERLAP MATRICES FOR INPUT JOBIPHS.
! PRORB  - PRINT INPUT ORBITALS.
! PRTRA  - PRINT TRANSFORMATION COEFFICIENTS.
! PRCI   - PRINT CI COEFFICIENTS OF INPUT STATES.
! CIH5   - Put CI coeffs and MO to HDF5
! IFHEXT - Spin-free Hamiltonian is taken from input.
! IFSHFT - Energy shifts of input states will be applied.
! IFHDIA - Diagonal H-matrix elements are taken from input.
! IFSO   - DO SPIN-ORBIT INTERACTION CALCULATION.

! Rassi input...
! RFpert - This flag is used to signal a
!          reaction field calculation (perturbation approach).
! ToFile - Denotes if H-matrix and various one-electron matrices
!          are to be put on a file for subsequent programs.
! PRXVR, etc: Print expectation values for RasScf input states,
!          for (spin-free) eigenstates, and for SO states.
! PRMER, etc: Print matrix elements    for RasScf input states,
!          for (spin-free) eigenstates, and for SO states.
! IfDCpl - Flag for approximate derivative couplings
!IgorS 06-05-2009
! HOP    - Switch for Trajectory Surface Hopping Algorithm
! stknecht
! QDPT2SC - use SC effective Hamiltonian (rather than the PC one) from QD-NEVPT2
! QDPT2EV - use eigenvectors of effective Hamiltonian from QD-NEVPT2 to mix TDMs (in MPS-SI we do not use 'mixed MPS'
!           instead we mix the TDMs)

! NTO Calculation Section /// Jie Bao
! SONTO            Array of SO state pairs
! SONTOSTATES      Number of state pairs to calculate
! SONAT            Array of SO state to compute
! SONATNSTATE      Number of states to calculate
! SODIAGNSTATE     Number of states to diagonalize

! TEMPORARY DATA FROM JOBIPHS
! POTNU1, NACTE1, MPLET1, NSYM1, NFRO1, NISH1, NASH1, NDEL1, NBAS1, NRS11, NRS21, NRS31, LROT1, NROOT1, IROOT1, NHOL11, NELE31,
! bNAME, HEAD1, TITLE1

! Cholesky/RI stuff
! ALGO, Nscreen, dmpk

! ISTATE1 - Number of the current relaxed state
! ISTATE2 - Number of the state that interacts with the current state
! NCI1    - Configuration's number of state 1
! ChkHop  - If .true. switches on the TSH algorithm

!----------------------------------------------------------------------*
!     Define files ( file names and unit numbers )                     *
!----------------------------------------------------------------------*

! LUIPH  - UNIT NUMBER OF JOBIPHS
! LUMCK  - UNIT NUMBER OF MCKINT FILES
! LUONE  - D:O, ONE-ELECTRON INTEGRAL FILE
! LUORD  - D:O, ORDERED TWO-ELECTRON INTEGRAL FILE
! IDCMO  - Addresses to the CMO arrays on each JOBIPH

integer(kind=iwp), parameter :: MORSBITS = 8, MXDISP = 500, MXJOB = 100, MXPROP = 30000

integer(kind=iwp) :: ALGO, DCHO, DYSEXPSF, DYSEXPSO, IBINA(2,MxRoot), ICOMP(MXPROP), IDCMO(MXJOB), IFJ2, IFJZ, IPUSED(MXPROP), &
                     IROOT1(mxRoot), IRREP(MXJOB), ISOCMP(MXPROP), ISTAT(MXJOB), ISTATE1, ISTATE2, ITOC15(30), L_Eff, LOOPDIVIDE, &
                     LOOPMAX, LROT1, LSYM1, LSYM2, LUEIG, LUEXC, LUIPH, LUMCK, LUONE, LUORD, LUTDM, LUTOM, MLTPLT(MXJOB), MPLET1, &
                     MULTIP, NACTE(MXJOB), NACTE1, NASH1(8), NATOMS, NBAS1(8), NBINA, NBSTEP, NCI1, NCI2, NCONF(MXJOB), NCONF1, &
                     NDEL1(8), NDET(MXJOB), NELE3(MXJOB), NELE31, NFRO1(8), NHOL11, NHOLE1(MXJOB), NISH1(8), NJOB, NPROP, nQuad, &
                     NRNATO, NROOT1, NROOTS(MXJOB), NRS11(8), NRS21(8), NRS31(8), Nscreen, NSOPR, NSOTHR_PRT, NSTAT(MXJOB), &
                     NSTATE, NSYM1, NTP, NTS, NTSTEP, OCAN, SODIAGNSTATE = 0, SONATNSTATE = 0, SONTOSTATES = 0
real(kind=wp) :: ALPHZ, BANGRES, BETAE, BINCRE, BSTART, CITHR, COOR(3,MXATOM), dmpk, EMIN, EPRATHR, EPRTHR, ERFNUC, OSTHR_DIPR, &
                 OSTHR_QIPR, PNUC(MXPROP) = Zero, PORIG(3,MXPROP) = Zero, RSTHR, SOTHR_PRT, TDIPMIN, TINCRE, TMAXP, TMAXS, &
                 TMGr_thrs, TMINP, TMINS, TOLERANCE, TSTART
logical(kind=iwp) :: BINA, ChkHop, CIH5, DCHS, DIPR, Do_Pol, Do_SK, Do_TMOM, DOCD, DoGSOR, DQVD, DYSEXPORT, DYSO, &
                     Force_NON_AO_TDM, HAVE_DIAG, HAVE_HEFF, HOP, IFACAL, IFACALFC, IFACALFCON, IFACALFCSDON, IFACALPSO, IFACALSD, &
                     IFACALSDON, IFARGU, IFATCALSA, IFCURD, IfDCpl, IFEJOB, IFGCAL, IFGTCALSA, IFGTSHSA, IFHAM, IFHCOM, IFHDIA, &
                     IFHEFF, IFHEXT, IFMCAL, IfNTO, IFSHFT, IFSO, IFSONCINI, IFTDM, IFTRD1, IFTRD2, IFVANVLECK, IFXCAL, LHAMI, &
                     LPRPR, NATO, NOHAM, NOSO, ONLY_OVERLAPS, PRCI, PRDIPCOM, PRDIPVEC, PRMEE, PRMER, PRMES, PRORB, PRRAW, PRSXY, &
                     PRTRA, PRWEIGHT, PRXVE, PRXVR, PRXVS, QDPT2EV, QDPT2SC, QIALL, QIPR, REDUCELOOP, RFpert, RHODyn, RSPR, &
                     SaveDens, SECOND_TIME, TDYS, ToFile, TRACK
character(len=LenIn+8) :: bNAME(mxOrb)
character(len=128) :: JBNAME(MXJOB), MINAME(MXJOB)
character(len=16) :: OCAA(20)
character(len=8) :: FnEig, FnToM, PNAME(MXPROP), PTYPE(MXPROP), RASTYP(MXJOB), SOPRNM(MXPROP), &
                    SOPRTP(MXPROP)
character(len=4) :: TITLE1(18,mxTit)
character(len=2) :: HEAD1(72)
integer(kind=iwp), allocatable :: SODIAG(:), SONAT(:), SONTO(:,:)
real(kind=wp), allocatable :: HEff(:,:), RefEne(:)

public :: ALGO, AlphZ, BAngRes, BetaE, BINA, BIncre, bNAME, BStart, ChkHop, CIH5, CITHR, Coor, DCHO, DCHS, DIPR, dmpk, Do_Pol, &
          Do_SK, DO_TMOM, DoCD, DOGSOR, DQVD, DYSEXPORT, DYSEXPSF, DYSEXPSO, DYSO, EMin, EPRATHR, EPRThr, ERFNuc, FnEig, FnTOM, &
          FORCE_NON_AO_TDM, HAVE_DIAG, HAVE_HEFF, HEAD1, HEff, HOP, IBINA, ICOMP, IDCMO, IFACAL, IFACALFC, IFACALFCON, &
          IFACALFCSDON, IFACALPSO, IFACALSD, IFACALSDON, IfArgu, IFATCALSA, IfCurd, IFDCPL, IFEJOB, IFGCAL, IFGTCALSA, IFGTSHSA, &
          IFHAM, IFHCOM, IFHDIA, IFHEFF, IFHEXT, IfJ2, IfJz, IFMCAL, IFNTO, IFSHFT, IFSO, IFSONCINI, IfTDM, IfTrD1, IFTRD2, &
          IfvanVleck, IFXCAL, IPUSED, IROOT1, IRREP, ISOCMP, ISTAT, ISTATE1, ISTATE2, iToc15, JBNAME, L_Eff, LHAMI, LOOPDIVIDE, &
          LOOPMAX, LPRPR, LROT1, lSym1, lSym2, LuEig, LuExc, LuIph, LuMck, LuOne, LuOrd, LuTDM, LUTOM, MINAME, MLTPLT, MORSBITS, &
          MPLET1, MULTIP, MXJOB, MXPROP, NACTE, NACTE1, NASH1, NATO, nAtoms, NBAS1, NBINA, NBSTep, nCI1, nCI2, NCONF, NCONF1, &
          NDEL1, NDET, NELE3, NELE31, NFRO1, NHOL11, NHOLE1, NISH1, NJOB, NOHAM, NOSO, NPROP, NQUAD, NrNATO, NROOT1, NROOTS, &
          NRS11, NRS21, NRS31, Nscreen, NSOPR, NSOThr_Prt, NSTAT, NSTATE, NSYM1, NTP, NTS, nTStep, OCAA, OCAN, ONLY_OVERLAPS, &
          OSThr_DiPr, OSThr_QIPR, PNAME, PNUC, PORIG, PRCI, PRDIPCOM, PrDipVec, PRMEE, PRMER, PRMES, PRORB, PrRaw, PRSXY, PRTRA, &
          PrWeight, PRXVE, PRXVR, PRXVS, PTYPE, QDPT2EV, QDPT2SC, QIAll, QIPR, RASTYP, REDUCELOOP, RefEne, RFPert, RhoDyn, RSPR, &
          RSThr, SAVEDENS, SECOND_TIME, SODIAG, SODIAGNSTATE, SONAT, SONATNSTATE, SONTO, SONTOSTATES, SOPRNM, SOPRTP, SOThr_Prt, &
          TDipMin, TDYS, TIncre, TITLE1, TMAXP, TMaxs, TMGR_Thrs, TMINP, TMins, ToFile, Tolerance, TRACK, TStart

end module Cntrl
