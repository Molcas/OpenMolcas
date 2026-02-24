!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) Per Ake Malmqvist                                      *
!***********************************************************************
!  RasScf_Init
!
!> @brief
!>   Initialize variables in commons, and set default values.
!>   Determine whether orbital files should be read, etc.
!> @author  P. &Aring;. Malmqvist
!>
!> @details
!> Sets values in the modules timers, rasscf_global and general_data.
!***********************************************************************

subroutine RasScf_Init()

use Fock_util_global, only: ALGO, Deco, DensityCheck, dmpk, DoCholesky, DoLocK, Estimate, Nscreen, Update
use casvb_global, only: ifvb
use Cholesky, only: ChFracMem, timings
use CMS, only: CMSGiveOpt, iCMSOpt
use UnixInfo, only: SuperName
use gas_data, only: IGSOCCX, NGAS, NGSSH
use timers, only: TimeAoMo, TimeCIOpt, TimeDavid, TimeDens, TimeFock, TimeHCSCE, TimeHDiag, TimeHSel, TimeInput, TimeOrb, &
                  TimePage, TimeRelax, TimeSigma, TimeTotal, TimeTrans, TimeWfn
use lucia_data, only: TDENSI, TSIGMA
use rasscf_global, only: CMSStartMat, CMSThreshold, CORESHIFT, Ener, ExFac, hRoots, iAlphaBeta, ICICH, ICICP, iCIonly, ICIRST, &
                         ICMSIterMax, ICMSIterMin, iCMSP, iExpand, IfCRPR, IfOrde, InOCalc, iOrbOnly, iOrbTyp, iOrdeM, iPCMRoot, &
                         iPhName, iPT2, iRLXRoot, IROOT, iRoot, irotPsi, iSave_Exp, iSPDen, iSupSM, itCore, ITMAX, iXMSP, ixSym, &
                         KSDFT, kTight, LowMS, LRoots, LvShft, MaxIt, MaxJT, MaxOrbOut, n_keep, NewFock, NonEq, NQUNE, NROOTS, &
                         OutFmt1, OutFmt2, PreThr, ProThr, PrwThr, Purify, QNSTEP, QNUPDT, RFPert, SXSel, ThFact, Thre, ThrEn, &
                         ThrSX, ThrTE, Title, TMin, Weight
use general_data, only: ISPIN, LOWDIN_ON, NACTEL, NALTER, NASH, NBAS, NDEL, NELEC3, NFRO, NHOLE1, NISH, NRS1, NRS2, NRS3, NRS3, &
                        NSEL, NSSH, STARTORBFILE, STSYM, SXDAMP
use spinfo, only: I_ELIMINATE_GAS_MOLCAS, ISPEED
use RASDim, only: MxCIIt, MxIter, MxSXIt
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: IPRGLB_IN, IPRLOC_IN(7)
integer(kind=iwp), external :: iPrintLevel
logical(kind=iwp), external :: Is_First_Iter

!----------------------------------------------------------------------*
! How was the program called?
!PAM 2009 Someone has put a number of possibilities here. Let it stand for now.
IfVB = 0
if (SuperName(1:6) == 'rasscf') then
  ICIRST = 0
  ! For geometry optimizations use the old CI coefficients.
  if (.not. Is_First_Iter()) ICIRST = 1
else if (SuperName(1:5) == 'casvb') then
  IfVB = 2
  ICIRST = 0
else if (SuperName(1:6) == 'loprop') then
  !ICIRST = 1 ! to be activated!
  ICIRST = 0
else if (SuperName(1:11) == 'last_energy') then
  ICIRST = 1
else if (SuperName(1:18) == 'numerical_gradient') then
  ICIRST = 1
else
  ICIRST = 0
end if

! Initialize print levels: Module output_ras
! Externally set default print level control. Should the program be silent?
IPRGLB_IN = iPrintLevel(-1)
IPRLOC_IN(:) = IPRGLB_IN
! Set print levels, and adjust them if needed:
call setprlev(IPRGLB_IN,IPRLOC_IN)

! Cholesky-related settings:
call DecideOnCholesky(DoCholesky)
ALGO = 1
DensityCheck = .false.
Deco = .true.
timings = .false.
DoLock = .true.
Nscreen = 10
dmpk = 1.0e-1_wp
Update = .true.
Estimate = .false.
!
#ifdef _MOLCAS_MPP_
ChFracMem = 0.3_wp
#else
ChFracMem = Zero
#endif

OutFmt1 = 'DEFAULT '
OutFmt2 = 'DEFAULT '

! Max nr of state-specific orbital files printed:
MAXORBOUT = 100
! Default title line:
TITLE(1) = '(No title given)'

! iteration control

! maximum number of RASSCF iterations
MAXIT = mxIter
! max number of super-CI iterations
ITMAX = mxSxIt
! max number of iterations in Davidson diagonalization
MAXJT = MXCIIT-2
! threshold for change in RASSCF energy
THRE = 1.0e-8_wp
!tbp, may 2013: no thre modification with Cholesky
!tbp if (DoCholesky) then
!tbp   call Get_dScalar('Cholesky Threshold',ThrCom)
!tbp   THRE = max(THRE,ThrCom)
!tbp end if
! threshold for max orbital rotation
!PAM2010 THRTE = 1.0e-4_wp
! PAM2010: Note: This is *not* a threshold that keeps rotation down
! between iterations in order to ensure proper function of the
! optimization -- it was intended as one of the thresholds that
! determine when the calculation has converged! As such, it is
! irrelevant! The relevant threshold is the max BLB.
THRTE = 1.0e-1_wp
! threshold for max BLB matrix element
! Note: If one changes the following value, please change it in
! fock_util/cho_LK_rassi.f and fock_util/cho_LK_rassi_x.f for consistency.
THRSX = 1.0e-4_wp
! Default damping in the SXCI orbital optimization
SXDAMP = 0.0002_wp
! Default thresholds used to determine convergence in CI
THREN = 1.0e-4_wp
THFACT = 1.0e-3_wp
! PAM 2017, Additional shift for douby occupied core states
! in order to compute core hole states. The core orbital is
! specified as one particular orbital in the input orbital set.
CORESHIFT = Zero
ITCORE = 0
IFCRPR = .false.
! PAM 2009, new default value for LVSHFT
! level shift parameter
LVSHFT = Half
! Quasi Newton update of the rotation matrix
NQUNE = 2
! only the CI calculation will be performed if iCIonly=1
iCIonly = 0
! only the orbitals from a JobIph to RasOrb if iOrbOnly=1
iOrbOnly = 0
! Default orbital type for RasOrb: Average orbitals
iOrbTyp = 1
! Root selection in the SXCI orbital optimization step.
! Values: LOWEST or HOMING.
SXSEL = 'LOWEST  '
! Choose to only expand or generate information for CI-vectors if INOCALC = 1
INOCALC = 0
! Save information on CI expansion if ISAVE_EXP = 1
ISAVE_EXP = 0
! Expand a smaller CI vector in a larger one if IEXPAND = 1
IEXPAND = 0

! wave function control bits

! new fock operator
NewFock = 1
! State used in response calculation
iPCMROOT = 1
! State to alaska
iRLXROOT = 0
! number of roots required in CI
NROOTS = 1
! number of roots actually used in CI-DAVIDSON
LROOTS = 1
! sequence numbers for roots in CI counted from
! lowest energy.
IROOT(1) = 1
iRoot(2:) = 0
! weights used for average energy calculations
WEIGHT(1) = One
WEIGHT(2:) = Zero
! iteration energies
ENER(:,:) = Zero

ICICH = 0
! if flag is active (ICICH=1) CI roots will be selected
! by maximum overlap with input CI function
! ICI(NROOTS,NREF)    CSF number for each root
! CCI(NROOTS,NREF)    corresponding CI coefficient
! maximum number is five csf's.

ISUPSM = 0
! make no use of supersymmetry
I_ELIMINATE_GAS_MOLCAS = 0
! Highly excited states are not default
hRoots = 0
! No hidden roots by default
n_keep = 0
! Number of kept vectors in Davidson chosen in ini_david by default
IORDEM = 0
! (SVC) do not force any ordering options
IFORDE = 1
! (SVC) use ordering of orbitals
! start CI Davidson with unit guess for CI vector
! use restart option if numerical gradients are computed.
PRWTHR = 0.05_wp
! threshold for printout of CI wave function

PROTHR = -One
! occupation threshold for printout of orbitals
! (The negative value serves to show if no user selection was made)
PRETHR = 999999.0_wp
! energy threshold for printout of orbitals

ICICP = 0
! no CI coupling (not active in this version)
NSEL = 200
! Default value for explicit Hamiltonian

TMIN = Zero
QNSTEP = 'SX'
QNUPDT = ' NO'
! Default value for tight parameter
KTIGHT = 0

! Default value for type of CASSCF (used for DFT)

KSDFT = 'SCF'
ExFac = One
!* Default orthonormalization of CMOs to be with
!* Gram-Schmidt
!Lowdin_ON = .false.
! PAM Jan 12 2010, on request, Lowdin ON has been made the default.
Lowdin_ON = .true.

! default for spin projection
LOWMS = 0
! default spin value (singlet)
ISPIN = 1
! default symmetry
STSYM = 1
! default number of active electrons
NACTEL = 0
! default maximum number of holes in RAS1
NHOLE1 = 0
! default maximum number of electrons in RAS3
NELEC3 = 0
! This run will not be the start for a CASPT2 calculation
IPT2 = 0
! This key will activate pertubational reaction field
! calculations.
RFpert = .false.
! Do compute the spin density matrix
ISPDEN = 1
! These keys will activate the calculation of the high
! frequency contribution to the reaction field
! ???
! This key controls if a non-equilibrium reaction field
! calculation is performed.
NonEq = .false.
! This initializes nr of input orbital swaps requested:
NAlter = 0
! set default values for orbitals

NFRO(:) = 0
NISH(:) = 0
NASH(:) = 0
NRS1(:) = 0
NRS2(:) = 0
NRS3(:) = 0
NSSH(:) = 0
NDEL(:) = 0
NBAS(:) = 0
! initialize occupation numbers for GAS

NGAS = 3
NGSSH = 0
IGSOCCX = 0
IXSYM(:) = 0
PURIFY = 'NO'

! Initial guess for jobiph name to use:
IPHNAME = 'JOBIPH'
! Initial guess for starting orbital file:
StartOrbFile = 'INPORB'
! Initialize alpha or beta orbitals (none):
iAlphaBeta = 0

! Initialize speed options (turn everything that's working on)
iSpeed(1:2) = 1
! The rest is at the present time just to allow testing
iSpeed(3:) = 0

TimeTotal = Zero
TimeInput = Zero
TimeWfn = Zero
TimeDens = Zero
TimeSigma = Zero
TimeHSel = Zero
TimeHDiag = Zero
TimeFock = Zero
TimeAoMo = Zero
TimeTrans = Zero
TimeCIOpt = Zero
TimeOrb = Zero
TimeDavid = Zero
TimePage = Zero
TimeHCSCE = Zero
TimeRelax = Zero

!SVC: lucia timers
tsigma(:) = Zero
tdensi(:) = Zero

! state rotation
iRotPsi = 0
iXMSP = 0
iCMSP = 0
ICMSIterMax = 100
ICMSIterMin = 5
CMSThreshold = 1.0e-8_wp
CMSStartMat = 'XMS'
iCMSOpt = 1
CMSGiveOpt = .false.

return

end subroutine RasScf_Init
