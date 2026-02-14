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
! Copyright (C) 1989, Per Ake Malmqvist                                *
!               1989, Bjorn O. Roos                                    *
!               1991,1993, Markus P. Fuelscher                         *
!               1991,1993, Jeppe Olsen                                 *
!               1998, Roland Lindh                                     *
!***********************************************************************

subroutine RASSCF(IRETURN)
!***********************************************************************
!                                                                      *
!           ######     #     #####   #####   #####  #######            *
!           #     #   # #   #     # #     # #     # #                  *
!           #     #  #   #  #       #       #       #                  *
!           ######  #     #  #####   #####  #       #####              *
!           #   #   #######       #       # #       #                  *
!           #    #  #     # #     # #     # #     # #                  *
!           #     # #     #  #####   #####   #####  #                  *
!                                                                      *
!                                                                      *
!                   A program for complete (CAS) and                   *
!                   restricted (RAS)SCF calculations                   *
!                                                                      *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M.P. Fuelscher and J. Olsen, P.Aa. Malmqvist and B.O. Roos       *
!     University of Lund, Sweden                                       *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history:                                                         *
!     MOLCAS version 1 by P.Aa. Malmqvist and B.O. Roos, 1989          *
!     MOLCAS version 2 by M.P. Fuelscher and J. Olsen  , 1991          *
!     MOLCAS version 3 by M.P. Fuelscher and J. Olsen  , 1993          *
!                                                                      *
!     Modified to process only unique symmetry blocks, R. Lindh,       *
!     March 1998.                                                      *
!                                                                      *
!***********************************************************************

use OneDat, only: sNoNuc, sNoOri
use Fock_util_global, only: ALGO, DoActive, DoCholesky
use write_orbital_files, only: OrbFiles, putOrbFile, write_orb_per_iter
use filesystem, only: copy_, real_path
use generic_CI, only: CI_solver_t
use fciqmc, only: DoNECI, fciqmc_solver_t, tGUGA_in
use fciqmc_read_RDM, only: dump_fciqmc_mats
use para_info, only: king
use fortran_strings, only: str
use spin_correlation, only: orb_range_p, orb_range_q, spin_correlation_driver
use CC_CI_mod, only: CC_CI_solver_t, Do_CC_CI
use fcidump, only: DumpOnly, make_fcidumps, transform
use orthonormalization, only: ON_scheme
use casvb_global, only: ifvb, invec_cvb
use OFembed, only: Do_OFemb, FMaux
use UnixInfo, only: ProgName
use rctfld_module, only: lRF
use Lucia_Interface, only: Lucia_Util
use wadr, only: CMO, D1A, D1I, DIAF, DMAT, DSPN, FA, FI, FockOcc, OccN, PA, PMAT, TUVX
use gugx, only: CIS, EXS, SGS
use gas_data, only: iDOGAS
use input_ras, only: KeyCION, KeyORBO, KeyORTH, KeySSCR, KeyTDM, KeyWRMA, LuInput
use raswfn, only: cre_raswfn, Wfn_FileID
use timers, only: TimeCIOpt, TimeInput, TimeOrb, TimeOutput, TimeRelax, TimeTotal, TimeTrans, TimeWfn
use rasscf_global, only: CBLBM, CMAX, Conv, DE, DOBLOCKDMRG, DoDMRG, DoFaro, DoFCIDump, ECAS, EMY, Ener, ESX, ExFac, FDIAG, HalfQ, &
                         iAdr15, iBLBM, ICICH, iCIOnly, iCIRST, iExpand, IfCrPr, InOCalc, IPCMROOT, iPr, iPT2, iRLXRoot, iRoot, &
                         iSave_Exp, iSymBB, ITER, ITERCI, ITERSX, JBLBM, KSDFT, KSDFT_Temp, l_casdft, lSquare, MaxIt, NAC, NACPAR, &
                         NACPR2, NewFock, nFint, no2m, NonEQ, nROOTS, PotNuc, QNSTEP, QNUPDT, ROTMax, Start_Vectors, SXShft, Thre, &
                         ThrSX, THRTE, TMin, Tot_Charge, VIA_DFT, Weight
use SplitCas_Data, only: DoSPlitCas, IterSplit, lRootSplit
use PrintLevel, only: DEBUG, TERSE, USUAL
use output_ras, only: IPRLOC, RC_CI, RC_SX
use general_data, only: CleanMask, CRPROJ, CRVec, INVEC, ISPIN, ITERFILE, JOBIPH, NALTER, NASH, NBAS, NCONF, NCRVEC, NDEL, NFRO, &
                        NISH, NRS1, NRS2, NRS3, NSYM, NTOT, NTOT1, NTOT2
use spinfo, only: DOBKAP
use DWSol, only: DWSol_final, DWSol_init, DWSolv
use Molcas, only: MxRoot
use RASDim, only: MxIter
#ifdef _DMRG_
use qcmaquis_interface, only: dmrg_energy, qcmaquis_interface_deinit, qcmaquis_interface_delete_chkp, &
                              qcmaquis_interface_prepare_hirdm_template, qcmaquis_param, TEMPLATE_4RDM, TEMPLATE_TRANSITION_3RDM
use qcmaquis_interface_mpssi, only: qcmaquis_mpssi_transform
use lucia_data, only: RF1, RF2
use rasscf_global, only: DoDelChk, DoMCPDFTDMRG, DoNEVPT2Prep, Twordm_qcm
#endif
#ifdef _FDE_
use Embedding_global, only: Eemb, embInt, embPot, embPotInBasis, embPotPath, embWriteEsp
#endif
#ifdef _HDF5_
use mh5, only: mh5_put_attr, mh5_put_dset
use csfbas, only: CONF
use lucia_data, only: CFTP, DStmp, Dtmp
use raswfn, only: wfn_energy, wfn_iter, wfn_transdens, wfn_transsdens
use rasscf_global, only: lRoots
use general_data, only: NACTEL, STSYM
#endif
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: IReturn
integer(kind=iwp) :: actual_iter, i, i_ROOT, iAd, iAd15, iBas, iComp, iFinal, iFlags, ihh, imm, Ind, IndT, IndType(56), iOff, &
                     iOpt, iPrLev, iRC, iRot, iShift, iss, iSyLbl, iSym, iTerm, j, kau, kRoot, LuOne, LuvvVec, mRoots, NoScr1, &
                     nTav, RC_RAS
real(kind=wp) :: CASDFT_E, CASDFT_FUNCT, DiffE, DiffETol, dum1, dum2, dum3, EAv, ECAS1, EVAC, ThMax, time0(2), time1(2), time2(2), &
                 time3(2), TMXTOT
logical(kind=iwp) :: DSCF, IfOpened, lOPTO, lTemp
character(len=80) :: Line, VecTyp
character(len=15) :: STLNE2
character(len=8) :: Label
character :: CTHRE, CTHRSX, CTHRTE
class(CI_solver_t), allocatable :: CI_solver
real(kind=wp), allocatable :: CMON(:), Dens(:), EDUM(:), Fock(:), folded_Fock(:), OCCX(:), orbital_E(:), PUVX(:), QMat(:), &
                              Scr1(:), Scr2(:), SMat(:), Tmp1(:), TmpD1S(:), TmpDMat(:), TmpDS(:)
#ifdef _HDF5_
integer(kind=iwp) :: iDX, jDisk, jRoot, kDisk
integer(kind=iwp), allocatable :: kcnf(:)
real(kind=wp), allocatable :: Tmp(:), VecL(:), VecR(:)
#endif
#ifdef _FDE_
integer(kind=iwp) :: iDummyEmb, iEmb, iUnit, nNuc
real(kind=wp), external :: EmbPotEneMODensities
#endif
#ifdef _DMRG_
integer(kind=iwp) :: maxBD, maxtrR
real(kind=wp) :: maxtrW
logical(kind=iwp) :: Do_ESPF
logical(kind=iwp), external :: PCM_On
#endif
integer(kind=iwp), external :: IsFreeUnit, isStructure
real(kind=wp), external :: Get_ExFac
#include "warnings.h"

! Set status line for monitor:
call StatusLine('RASSCF: ','Just started.')

! Set the return code(s)
ITERM = 0
IRETURN = _RC_ALL_IS_WELL_

! Set the HDF5 file id (a proper id will never be 0)
wfn_fileid = 0

! Set some Cholesky stuff
DoActive = .true.
lOPTO = .false.
! Initialise doDMRG if compiled without QCMaquis
#ifndef _DMRG_
DoDMRG = .false.
#endif

! Set variable IfVB to check if this is a VB job.
IfVB = 0
if (ProgName(1:5) == 'casvb') IfVB = 2
! Default option switches and values, and initial data.
THMAX = Zero
call RasScf_Init()
call Seward_Init()
! Open the one-olectron integral file:
LuOne = 77
LuOne = isFreeUnit(LuOne)
iRC = -1
iOpt = 0
call OpnOne(iRC,iOpt,'ONEINT',LuOne)
if (iRC /= 0) then
  write(u6,*) 'Error when trying to open the one-electron'
  write(u6,*) 'integral file.'
  call Quit(_RC_INTERNAL_ERROR_)
end if
call StatusLine('RASSCF: ','Read-in ONEINT')
if (IfVB == 2) go to 10

! Make a copy, upper-cased, left-adjusted, of the input between and including
! the '&RASSCF' and the 'End of input' markers, skipping all lines beginning
! with '*' or '!' or ' '  when left-adjusted, and replacing any rightmost
! substring beginning with '!' with blanks.
! That copy will be in file 'CleanInput', and its unit number is returned
! as LUInput in common (module file input_ras.F90) by the following call:
call cpinp(LUInput,iRc)
! If something wrong with input file:
if (iRc /= _RC_ALL_IS_WELL_) then
  call WarningMessage(2,'Input file is unusable.')
  write(u6,*) ' RASSCF Error: Could not make a clean copy of'
  write(u6,*) ' the input file. This is an unexpected bug.'
  IRETURN = _RC_INTERNAL_ERROR_
  goto 9990
end if

! Scan the input file for keywords:
call Scan_Inp(iRc)
! If something wrong with input file:
if (iRc /= _RC_ALL_IS_WELL_) then
  if (IPRLOC(1) >= TERSE) then
    call WarningMessage(2,'Scanning input file failed.')
    ! Calling again, now with iRc indicating an error, will echo the keywords:
    call Scan_Inp(iRc)
  end if
  IRETURN = _RC_INPUT_ERROR_
  goto 9990
end if

! Local print level in this routine:
IPRLEV = IPRLOC(1)

10 continue
! Open files
call OpnFls_RASSCF(DSCF,DoCholesky)

! Some preliminary input data:
call Rd1Int()
if (.not. DSCF) call Rd2Int_RASSCF()

! Printed program header:

! Process the input:
call StatusLine('RASSCF: ','Processing input')
call Proc_Inp(DSCF,lOPTO,iRc)
! If something goes wrong in proc_inp:
if (iRc /= _RC_ALL_IS_WELL_) then
  if (IPRLEV >= TERSE) then
    call WarningMessage(2,'Input processing failed.')
    write(u6,*) ' RASSCF Error: Proc_Inp failed unexpectedly.'
    write(u6,*) ' Check the output file for any previous messages'
    write(u6,*) ' that can help explain the failure.'
    write(u6,*) ' Here is a printing of the input file that'
    write(u6,*) ' was processed:'
    rewind(LUInput)
15  continue
    read(LuInput,'(A80)',end=16,err=16) Line
    write(u6,*) Line
    Go To 15
16  continue
  end if
  IRETURN = iRc
  goto 9990
end if
if (lRF) call DWSol_init(IPCMROOT,nRoots,NonEq)
!call DWSCF_init(1,nRoots)

! Local print level may have changed:
IPRLEV = IPRLOC(1)

call InpPri(lOpto)

! Note that CI_solver subclasses provide a cleanup procedure
! (C++ people might call it destructor). Hence the deallocation and
! cleanup is automatically performed, when it goes out of scope.
if (DoNECI) then
  allocate(CI_solver,source=fciqmc_solver_t(tGUGA_in))
else if (Do_CC_CI) then
  allocate(CI_solver,source=CC_CI_solver_t())
end if

! If this is not CASDFT make sure the DFT flag is unset

if (KSDFT(1:3) == 'SCF') then
  call Get_iScalar('System BitSwitch',iFlags)
  iFlags = iand(iFlags,not(2**6))
  call Put_iScalar('System BitSwitch',iFlags)
end if

! If the ORBONLY option was chosen, then Proc_Inp just generated
! orbitals from the JOBIPH file. Nothing more to do:
if (KeyORBO .or. (MAXIT == 0)) goto 9989
#ifdef _DMRG_
! delete old checkpoints, unless requested otherwise
! this flag is set in proc_inp
if (.not. DoDelChk) then
  do kroot=1,nroots
    call qcmaquis_interface_delete_chkp(iroot(kroot))
  end do
end if
#endif

! Allocate various matrices

call mma_allocate(FI,NTOT1,Label='FI')
call mma_allocate(FA,NTOT1,Label='FA')
call mma_allocate(D1I,NTOT2,Label='D1I')
call mma_allocate(D1A,NTOT2,Label='D1A')
call mma_allocate(OCCN,NTOT,Label='OccN')
call mma_allocate(CMO,NTOT2,Label='CMO')
call mma_allocate(DIAF,NTOT,Label='DIAF')
#ifdef _DMRG_
! Allocate RDMs for the reaction field reference root in QCMaquis calculations
if (doDMRG .and. PCM_On()) then
  call mma_allocate(RF1,NACPAR,Label='RF1')
  if (twordm_qcm) call mma_allocate(RF2,NACPR2,Label='RF2')
end if
#endif
FI(:) = Zero
FA(:) = Zero
DIAF(:) = Zero
ECAS1 = Zero
EVAC = Zero

if ((iCIRST == 1) .and. DumpOnly) then
  write(u6,*) 'ICIRST and DumpOnly flags are not compatible!'
  write(u6,*) 'Choose only one.'
  call Abend()
end if

if (DumpOnly) then
  write(u6,*) 'Dumping integrals.'
  write(u6,*) 'Nothing else will be done.'
end if

call mma_allocate(TUVX,NACPR2,Label='TUVX')
TUVX(:) = Zero
call mma_allocate(DSPN,NACPAR,Label='DSPN')
DSPN(:) = Zero
call mma_allocate(DMAT,NACPAR,Label='DMat')
DMAT(:) = Zero
call mma_allocate(PMAT,NACPR2,Label='PMat')
PMAT(:) = Zero
call mma_allocate(PA,NACPR2,Label='PA')
PA(:) = Zero
#ifdef _FDE_
! Embedding
iDummyEmb = 0
call Get_iScalar('embpot',iDummyEmb)
if (iDummyEmb == 1) embPot = .true.
if (embPot) call EmbPotRdRun()
if (embpot) then
  ! I have no idea why i need memory for x+4 entries
  ! and not just x...
  call mma_allocate(embInt,NTOT1+4,label='Emb')
  if (embPotInBasis) then
    ! If the potential is given in basis set representation it
    ! has not been calculated with a OneEl call and is just read
    ! from file here.
    iunit = isFreeUnit(1)
    call molcas_open(iunit,embPotPath)
    do iEmb=1,NTOT1
      read(iunit,*) embInt(iEmb)
    end do
  else
    ! Read in the embedding potential one-electron integrals
    Label = 'embpot  '
    iRC = -1
    iOpt = 0
    iComp = 1
    call RdOne(iRC,iOpt,Label,iComp,embInt,iSyLbl)
    if (iRC /= 0) then
      call WarningMessage(2,'Drv1El: Error reading ONEINT;Label='//Label)
      call Quit(_RC_IO_ERROR_READ_)
    end if
  end if
end if
#endif

! Get start orbitals

! Initialize OCCN array, to prevent false alarms later from
! automated detection of using uninitialized variables:
OccN(:) = Zero

! PAM03: Note that removal of linear dependence may change the nr
! of secondary/deleted orbitals, affecting some of the global
! variables: NSSH(),NDEL(),NORB(),NTOT3, etc etc
call ReadVc(CMO,OCCN,DMAT,DSPN,PMAT,PA,ON_scheme)
! TODO(Oskar): Add fourth argument OCC
! If the Occupation number is written properly as well.
if (KeyORTH) call putOrbFile(CMO=CMO(:),orbital_E=DIAF(:),iDoGAS=iDoGAS)
! Only now are such variables finally known.

! Allocate core space for dynamic storage of data

call ALLOC()

! Create job interphase on unit JOBIPH (FT15)

if (ifvb /= 2) then
  call CREIPH()
  call cre_raswfn()
end if
if (ifvb == 1) call casinfo2_cvb()

call Timing(dum1,dum2,time1(1),dum3)
TimeInput = time1(1)

!GG03 Aug 03
if (NAlter > 0) call Alter_MO(CMO)

! At this point all is ready to potentially dump MO integrals... just do it if required.
if (DumpOnly) goto 20
if (ifvb == 2) goto 20

if (dofcidump) then
  write(u6,*)
  write(u6,'(26X,A)') 'Dumping integrals on file FCIDUMP - nothing else to be done'
  goto 20
end if

!***********************************************************************
!
! Wave function section
!
!***********************************************************************

call StatusLine('RASSCF: ','Compute wave function.')
if ((IPRLEV >= 2) .and. (.not. lOPTO)) then
  write(u6,*)
  write(u6,'(6X,A)') repeat('*',120)
  write(u6,'(6X,A,118X,A)') '*','*'
  write(u6,'(6X,A,44X,A,45X,A)') '*','Wave function control section','*'
  write(u6,'(6X,A,118X,A)') '*','*'
  write(u6,'(6X,A)') repeat('*',120)
  write(u6,*)
end if

if ((IPRLEV >= 2) .and. (.not. lOPTO)) then
  if (ICIONLY == 0) then
    write(u6,*)
    if (doDMRG) then
      write(u6,'(41X,A)') 'DMRGSCF iterations: Energy and convergence statistics'
      write(u6,'(41X,A)') '-----------------------------------------------------'
    else
      write(u6,'(41X,A)') 'RASSCF iterations: Energy and convergence statistics'
      write(u6,'(41X,A)') '----------------------------------------------------'
    end if
    write(u6,*)
  else
    write(u6,*)
    if (doDMRG) then
      write(u6,'(41X,A)') 'DMRGCI only, no orbital optimization will be done.'
      write(u6,'(41X,A)') '--------------------------------------------------'
    else
      write(u6,'(41X,A)') 'CASCI only, no orbital optimization will be done.'
      write(u6,'(41X,A)') '-------------------------------------------------'
    end if
    write(u6,*)
  end if
# ifdef _DMRG_
  if (doDMRG) then
    write(u6,'(45x,a//,36x,a/,36x,a/,36x,a//,45x,a//,36x,a/,36x,a/,36x,a//,36x,a/,36x,a,a/,36x,a//)') &
      'Please cite for the QCMaquis-Molcas driver:', &
      'Freitag L.; Keller S.; Knecht S.; Lindh R.; Ma Y.; ', &
      'Stein C. J. and Reiher M., in preparation. (2018).', &
      '---------------------------------------------------------------', &
      'Please cite for the QCMaquis DMRG software:', &
      'S. Keller, M. Dolfi, M. Troyer, M. Reiher,', &
      'J. Chem. Phys. 143, 244118 (2015)', &
      '---------------------------------------------------------------'
  end if
# endif
  if (INOCALC == 1) then
    write(u6,*)
    write(u6,'(26X,A)') ' No calculation will be performed. Stopping in LUCIA'
  end if
  if (ISAVE_EXP == 1) then
    write(u6,*)
    write(u6,'(26X,A)') ' Information on the CI-vector will be written to ???.'
  end if
  if (IEXPAND == 1) then
    write(u6,*)
    write(u6,'(26X,A)') ' A shorter vector will be expanded in a longer'
  end if
  if (IPRLEV <= 3) then
    if (DoSplitCAS) then
      write(u6,'(6X,A)') 'Iter CI   SX   CI       SplitCAS       Energy    max ROT     max BLB   max BLB  Level Ln srch  '// &
                         'Step   QN   Walltime'
      write(u6,'(6X,A)') '    iter iter root      energy       change     param      element    value   shift minimum  '// &
                         'type update hh:mm:ss'
    else if (DoBKAP) then
      write(u6,'(6X,A)') 'Iter CI   SX   CI   RASSCF      CI    Energy    max ROT     max BLB   max BLB  Level Ln srch  '// &
                         'Step   QN   Walltime'
      write(u6,'(6X,A)') '    iter iter root  energy    energy  change     param      element    value   shift minimum  '// &
                         'type update hh:mm:ss'
    else if (DoDMRG .and. (ICIONLY == 0)) then
      write(u6,'(6X,A)') 'Iter num   Bond  DMRG max tr DMRG  SX      DMRGSCF       Energy    max ROT   max BLB     max BLB  '// &
                         'Level Ln srch  Step   QN     CPU Time'
      write(u6,'(6X,A)') '   sweeps/ dim  /root weight/root iter     energy        change     param    element      value   '// &
                         'shift minimum  type update   hh:mm:ss'
    else if (DoDMRG .and. (ICIONLY /= 0)) then

    else if (l_casdft) then

    else
      write(u6,'(6X,A)') 'Iter CI   SX   CI       RASSCF       Energy    max ROT     max BLB   max BLB  Level Ln srch  '// &
                         'Step   QN   Walltime'
      write(u6,'(6X,A)') '    iter iter root      energy       change     param      element    value   shift minimum  '// &
                         'type update hh:mm:ss'
    end if
  end if
end if
20 continue
!                                                                      *
!***********************************************************************
!                                                                      *
!     Start iterations
!                                                                      *
!***********************************************************************
!                                                                      *
Rc_CI = 0
Rc_SX = 0
ECAS = Zero
ROTMAX = Zero
ITER = 0
actual_iter = 0
IFINAL = 0
TMXTOT = Zero
call mma_allocate(FockOcc,nTot1,Label='FockOcc')
!                                                                      *
!***********************************************************************
!                                                                      *
!     Entry point for second and successive iterations
!                                                                      *
!***********************************************************************
!                                                                      *
1000 continue

if (l_casdft) then
  KSDFT_TEMP = KSDFT
  KSDFT = 'SCF'
  ExFac = One
else
  KSDFT_TEMP = KSDFT
  ExFac = Get_ExFac(KSDFT)
end if

ITER = ITER+1
write(STLNE2,'(A12,I3)') 'Iteration ',ITER
call StatusLine('RASSCF: ',STLNE2)
call Timing(dum1,dum2,time0(1),dum3)
#ifdef _DMRG_
! Leon 27/11/2017: Skip the first CI iteration if we're using
! DMRGCI and CIOnly.It's enabled only for DMRGCI with QCMaquis
! now, (to exclude potential side effects)
! but consider extending it to other cases!
call DecideOnESPF(Do_ESPF)
!write(u6,*) ' |rasscf> DecideOnESPF == ',Do_ESPF
if ((ITER == 1) .and. ((.not. (DoDMRG .and. (ICIONLY /= 0))) .or. lRf .or. domcpdftDMRG .or. Do_ESPF)) then
# else
if (ITER == 1) then
# endif
  !*********************************************************************
  !     ^   First iteration
  !*********************************************************************

  ! Print header to file containing informations on CI iterations.

  write(IterFile,'(A)') repeat('*',80)
  write(IterFile,'(15X,A)') 'RASSCF iteration: 1A'

  Start_Vectors = .true.
  lTemp = lRf

  ! Transform two-electron integrals and compute at the same time
  ! the Fock matrices FI and FA

  call Timing(dum1,dum2,time2(1),dum3)

  if ((.not. DoCholesky) .or. (ALGO == 1)) then
    call mma_allocate(PUVX,NFINT,Label='PUVX')
    PUVX(:) = Zero
  end if

  call Get_D1I_RASSCF(CMO,D1I)
  if (IPRLEV >= DEBUG) then
    write(u6,*)
    write(u6,*) ' D1I in AO basis in RASSCF'
    write(u6,*) ' ---------------------'
    write(u6,*)
    iOff = 1
    do iSym=1,nSym
      iBas = nBas(iSym)
      call wrtmat(D1I(ioff),iBas,iBas,iBas,iBas)
      iOff = iOff+iBas*iBas
    end do
  end if

  ! Compute D1A from CMO coefficients and, if CIREstart, old DMAT.
  if (iCIRST == 1) then

    call mma_allocate(TmpDMAT,NACPAR,Label='TmpDMAT')
    call dcopy_(NACPAR,DMAT,1,TmpDMAT,1)
    if (NASH(1) /= NAC) call DBLOCK(TmpDMAT)
    call Get_D1A_RASSCF(CMO,TmpDMAT,D1A)
    call mma_deallocate(TmpDMAT)

    DoActive = .true.

  else

    lRf = .false.
    if (.not. l_casdft) then
      KSDFT = 'SCF'
      ExFac = One
    end if
    D1A(:) = Zero

    DoActive = .false.

  end if

  IPR = 0
  if (IPRLOC(2) == 4) IPR = 5
  if (IPRLOC(2) == 5) IPR = 10

  if (IPRLEV >= DEBUG) then
    write(u6,*)
    write(u6,*) ' PUVX in rasscf bf first TRACTL2'
    write(u6,*) ' ---------------------'
    write(u6,*)
    call wrtmat(PUVX,1,nFint,1,nFint)

    write(u6,*)
    write(u6,*) ' ---------------------'
    write(u6,*)
    write(u6,*) ' D1A in AO basis in RASSCF bf TRACTL2 1'
    write(u6,*) ' ---------------------'
    write(u6,*)
    iOff = 1
    do iSym=1,nSym
      iBas = nBas(iSym)
      call wrtmat(D1A(ioff),iBas,iBas,iBas,iBas)
      iOff = iOff+iBas*iBas
    end do
  end if

  ! Transform two-electron integrals and compute the Fock matrices FI and FA
  ! FI and FA are output from TRACTL2...
  call TRACTL2(CMO,PUVX,TUVX,D1I,FI,D1A,FA,IPR,lSquare,ExFac)

  !write(u6,*) ' TUVX after TRACTL2'
  !write(u6,*) (UVX(ind),ind=1,NACPR2)
  ! Core shift applied to projection of WF with doubly occupied core
  if ((ITER == 1) .and. IfCRPR) call MkCRVEC(CMO,CRVEC)

  if (IPRLEV >= DEBUG) then
    write(u6,*)
    write(u6,*) ' D1A in AO basis in RASSCF af TRACTL2 1'
    write(u6,*) ' ---------------------'
    write(u6,*)
    iOff = 1
    do iSym=1,nSym
      iBas = nBas(iSym)
      call wrtmat(D1A(ioff),iBas,iBas,iBas,iBas)
      iOff = iOff+iBas*iBas
    end do

    write(u6,*)
    write(u6,*) ' ---------------------'
    write(u6,*)
    write(u6,*) ' PUVX in rasscf af first TRACTL2'
    write(u6,*) ' ---------------------'
    write(u6,*)
    call wrtmat(PUVX,1,nFint,1,nFint)

    write(u6,*)
    write(u6,*) ' ---------------------'

    write(u6,*) ' TUVX after TRACTL2'
    write(u6,*) (TUVX(ind),ind=1,NACPR2)
    write(u6,*)
    write(u6,*) ' ---------------------'
  end if

  if ((.not. DoCholesky) .or. (ALGO == 1)) call mma_deallocate(PUVX)

  call Timing(dum1,dum2,time2(2),dum3)
  TimeTrans = TimeTrans+time2(2)-time2(1)

  if (IPRLEV >= DEBUG) then
    write(u6,*)
    write(u6,*) ' CMO in RASSCF bf first call to CICTL'
    write(u6,*) ' ---------------------'
    write(u6,*)
    ioff = 1
    do iSym=1,nSym
      iBas = nBas(iSym)
      if (iBas /= 0) then
        write(u6,*) 'Sym =',iSym
        do i=1,iBas
          write(u6,*) (CMO(ioff+iBas*(i-1)+j),j=0,iBas-1)
        end do
        iOff = iOff+(iBas*iBas)
      end if
    end do
  end if

  ! Compute initial CI vectors and density matrices

  call Timing(dum1,dum2,time3(1),dum3)

  if (DumpOnly) then
    call mma_allocate(orbital_E,nTot)
    call mma_allocate(folded_Fock,nAcPar)
    call transform(iter, &
                   CMO=CMO(:), &
                   DIAF=DIAF(:), &
                   D1I_AO=D1I(:), &
                   D1A_AO=D1A(:), &
                   D1S_MO=DSPN(:), &
                   F_IN=FI(:), &
                   orbital_E=orbital_E, &
                   folded_Fock=folded_Fock)
    call make_fcidumps('FCIDUMP','H5FCIDUMP',orbital_E,folded_Fock,TUVX=tuvx(:),core_energy=EMY)
    call mma_deallocate(orbital_E)
    call mma_deallocate(folded_Fock)
    write(u6,*) 'FCIDMP file generated. Here for serving you!'
    goto 2010
  end if

  if (allocated(CI_solver)) then
    call CI_solver%run(actual_iter=actual_iter, &
                       ifinal=ifinal, &
                       iroot=iroot, &
                       weight=weight, &
                       CMO=CMO(:), &
                       DIAF=DIAF(:), &
                       D1I_AO=D1I(:), &
                       D1A_AO=D1A(:), &
                       TUVX=tuvx(:), &
                       F_IN=FI(:), &
                       D1S_MO=DSPN(:), &
                       DMAT=DMAT(:), &
                       PSMAT=pmat(:), &
                       PAMAT=pa(:))

# if defined (_ENABLE_BLOCK_DMRG_) || defined (_ENABLE_CHEMPS2_DMRG_) || defined (_ENABLE_DICE_SHCI_)
  else if (DoBlockDMRG) then
    call DMRGCTL(CMO,DMAT,DSPN,PMAT,PA,FI,D1I,D1A,TUVX,IFINAL,0)
# endif
  else
    call CICTL(CMO,DMAT,DSPN,PMAT,PA,FI,FA,D1I,D1A,TUVX,IFINAL)

    if (dofcidump) then
      write(u6,*) ' FCIDUMP file generated. This is the end...'
      goto 9990
    end if
#   ifdef _FDE_
    !Thomas Dresselhaus
    if (embpot) then
      !Eemb = DDot_(NACPAR,embInt,1,DMAT,1)
      !Eemb = embPotEne(D1I,D1A,embInt,CMO,nBasFunc,nFrozenOrbs,.true.)
      Eemb = embPotEneMODensities(D1I,D1A,embInt,nBas,nTot2,nSym)
      write(u6,*) 'Energy from embedding potential with the'
      write(u6,*) 'initial CI vectors: ',Eemb
    end if
    !!!!!!!!!!!!!!!!!!!
#   endif
    ! PAM 2015: Additional output line.
    if ((IPRLEV >= USUAL) .and. (.not. doDMRG)) write(u6,'(a,i4)') ' Nr of preliminary CI iterations:',ITERCI
  end if

  !.. dongxia testing jobiph
  !.. upt to here, jobiph are all zeros at iadr15(2)
  ! If CASVB job, go directly to return.
  invec_cvb = invec
  if (DSCF) NewFock = 1
  if (IfVB == 2) goto 9990

  EAV = Zero
  if (DoSplitCAS) then
    EAV = ENER(lRootSplit,ITER)
  else
    do KROOT=1,NROOTS
      EAV = EAV+ENER(IROOT(KROOT),ITER)*WEIGHT(KROOT)
    end do
  end if

  call Get_D1A_RASSCF(CMO,DMAT,D1A)

  if (IPRLEV >= DEBUG) then
    write(u6,*)
    write(u6,*) ' D1A in AO basis in RASSCF af Get_D1A_RASSCF'
    write(u6,*) ' ---------------------'
    write(u6,*)
    iOff = 1
    do iSym=1,nSym
      iBas = nBas(iSym)
      call wrtmat(D1A(ioff),iBas,iBas,iBas,iBas)
      iOff = iOff+iBas*iBas
    end do
  end if

  call Timing(dum1,dum2,time3(2),dum3)
  TimeCIOpt = TimeCIOpt+time3(2)-time3(1)
  lRf = lTemp

  if (.not. l_casdft) then
    KSDFT = KSDFT_TEMP
    ExFac = Get_ExFac(KSDFT)
  end if

  ! v GLM for MC-PDFT
  if ((KSDFT /= 'SCF') .and. (KSDFT /= 'PAM') .or. l_casdft) then
    if (IPRLEV >= DEBUG) then
      write(u6,*)
      write(u6,*) ' CMO in RASSCF bf call NATORB_RASSCF'
      write(u6,*) ' ---------------------'
      write(u6,*)
      ioff = 1
      do iSym=1,nSym
        iBas = nBas(iSym)
        if (iBas /= 0) then
          write(u6,*) 'Sym =',iSym
          do i=1,iBas
            write(u6,*) (CMO(ioff+iBas*(i-1)+j),j=0,iBas-1)
          end do
          iOff = iOff+(iBas*iBas)
        end if
      end do
    end if

    call mma_allocate(CMON,NTOT2,Label='CMON')
    call mma_allocate(OCCX,NTOT,Label='OCCX')
    noscr1 = max(nacpar,no2m)
    call mma_allocate(scr1,noscr1,Label='Scr1')
    call mma_allocate(scr2,NO2M,Label='Scr2')
    call mma_allocate(SMAT,NTOT1,Label='SMAT')
    call NATORB_RASSCF(CMO,scr1,scr2,SMAT,CMON,OCCX)
    call dCopy_(NTOT2,CMON,1,CMO,1)
    call Put_dArray('Last orbitals',CMO,ntot2)
    call mma_deallocate(scr1)
    call mma_deallocate(scr2)
    call mma_deallocate(SMAT)
    call mma_deallocate(OCCX)
    call mma_deallocate(CMON)
    if (IPRLEV >= DEBUG) then
      write(u6,*)
      write(u6,*) ' CMO in RASSCF af call NATORB_RASSCF & bf 2 CICTL'
      write(u6,*) ' ---------------------'
      write(u6,*)
      ioff = 1
      do iSym=1,nSym
        iBas = nBas(iSym)
        if (iBas /= 0) then
          write(u6,*) 'Sym =',iSym
          do i=1,iBas
            write(u6,*) (CMO(ioff+iBas*(i-1)+j),j=0,iBas-1)
          end do
          iOff = iOff+(iBas*iBas)
        end if
      end do
    end if
  end if
  !^ GLM End If for MC-PDFT

  EAV = Zero
  do KROOT=1,NROOTS
    EAV = EAV+ENER(IROOT(KROOT),ITER)*WEIGHT(KROOT)
  end do
  if (IPRLEV >= DEBUG) then
    if (l_casdft) then
      write(u6,*) 'EAV value in RASSCF after first call to CICTL:'
      write(u6,*) EAV
    end if
  end if
end if
!***********************************************************************
!^ End First iteration
!***********************************************************************

! Print header to file containing informations on CI iterations.

actual_iter = actual_iter+1
write(IterFile,*)
write(IterFile,'(A)') repeat('*',80)
if (Iter == 1) then
  write(IterFile,'(15X,A)') 'RASSCF iteration: 1B'
else
  write(IterFile,'(15X,A,I3)') 'RASSCF iteration: ',Iter
end if

if ((IPRLEV >= DEBUG) .and. l_casdft) then
  write(u6,*) repeat('*',70)
  write(u6,*) 'we are done withe first standard CAS-CI iteration  '
  write(u6,*) 'CI coeffs are known and mantained fix in next stage'
  write(u6,*) 'We are now going to remove exchange from FI and FA '
  write(u6,*) 'in TRACTL2 --> TRA_CTL2 --> TRADRV --> FTWO        '
  write(u6,*) 'the ExFac is going to be set to 0.0 as DT asked! '
  write(u6,*) 'FI and FA are going to change... '
  write(u6,*) 'Check with previous printout to see differences.\  '
  write(u6,*) repeat('*',70)
end if

if (l_casdft) then
  KSDFT = KSDFT_TEMP
  ExFac = Zero
end if

if (ICIONLY /= 0) IFINAL = 1

! Transform two-electron integrals and compute at the same time
! the Fock matrices FI and FA

call Timing(dum1,dum2,time2(1),dum3)
if (.not. DoCholesky .or. (ALGO == 1)) then
  call mma_allocate(PUVX,NFINT,Label='PUVX')
  PUVX(:) = Zero
end if
call Get_D1I_RASSCF(CMO,D1I)

DoActive = .true.

if (DoCholesky .and. (ALGO == 2)) then
  NTav = 0
  do iSym=1,nSym
    NTav = NTav+nBas(iSym)*nAsh(iSym)
  end do
  call mma_allocate(Qmat,NTav,Label='QMat')
  QMat(:) = Zero
end if

IPR = 0
if (IPRLOC(2) == 4) IPR = 5
if (IPRLOC(2) == 5) IPR = 10
if (IPRLEV >= DEBUG) then
  write(u6,*)
  write(u6,*) ' D1A in AO basis in RASSCF bf TRACTL2 2'
  write(u6,*) ' ---------------------'
  write(u6,*)
  iOff = 1
  do iSym=1,nSym
    iBas = nBas(iSym)
    call wrtmat(D1A(ioff),iBas,iBas,iBas,iBas)
    iOff = iOff+iBas*iBas
  end do
end if
call TRACTL2(CMO,PUVX,TUVX,D1I,FI,D1A,FA,IPR,lSquare,ExFac)

if (IPRLEV >= DEBUG) then
  write(u6,*)
  write(u6,*) ' D1A in AO basis in RASSCF af TRACTL2 2'
  write(u6,*) ' ---------------------'
  write(u6,*)
  iOff = 1
  do iSym=1,nSym
    iBas = nBas(iSym)
    call wrtmat(D1A(ioff),iBas,iBas,iBas,iBas)
    iOff = iOff+iBas*iBas
  end do
end if

if ((.not. DoCholesky) .or. (ALGO == 1)) call mma_deallocate(PUVX)

call Timing(dum1,dum2,time2(2),dum3)
TimeTrans = TimeTrans+time2(2)-time2(1)

! Compute the CI vectors and density matrices

call Timing(dum1,dum2,time3(1),dum3)
if (.not. l_casdft) then !the following is skipped in CASDFT-GLM

  if ((KSDFT /= 'SCF') .and. (KSDFT /= 'PAM')) call Put_dArray('Last orbitals',CMO,ntot2)

  if (allocated(CI_solver)) then
    ! The following is adapted from sxctl.f
    ! In addition to writing the last RasOrb to disk, the current
    ! orbitals have to be dumped *before* the CI step. The PERI
    ! keyword writes only the orbitals from the last iteration.
    iShift = 0
    do ISYM=1,NSYM
      IndT = 0
      IndType(1+iShift) = NFRO(ISYM)
      IndT = IndT+NFRO(ISYM)
      IndType(2+iShift) = NISH(ISYM)
      IndT = IndT+NISH(ISYM)
      IndType(3+iShift) = NRS1(ISYM)
      IndT = IndT+NRS1(ISYM)
      IndType(4+iShift) = NRS2(ISYM)
      IndT = IndT+NRS2(ISYM)
      IndType(5+iShift) = NRS3(ISYM)
      IndT = IndT+NRS3(ISYM)
      IndType(7+iShift) = NDEL(ISYM)
      IndT = IndT+NDEL(ISYM)
      IndType(6+iShift) = NBAS(ISYM)-IndT
      iShift = iShift+7
    end do
    call mma_allocate(EDUM,NTOT,Label='EDum')
    EDum(:) = Zero
    write(VecTyp,'(A)')
    VecTyp = '* RASSCF average (pseudo-natural) orbitals (Not final)'
    LuvvVec = 50
    LuvvVec = isfreeunit(LuvvVec)
    call WrVec('IterOrb',LuvvVec,'COE',NSYM,NBAS,NBAS,CMO(:),OCCN,EDUM,INDTYPE,VECTYP)
    call WrVec('IterOrb',LuvvVec,'AI',NSYM,NBAS,NBAS,CMO(:),OCCN,EDUM,INDTYPE,VECTYP)
    call mma_deallocate(EDUM)
    write(u6,*) 'MO coeffs for next iteration written to IterOrb.'

    call CI_solver%run(actual_iter=actual_iter, &
                       ifinal=ifinal, &
                       iroot=iroot, &
                       weight=weight, &
                       CMO=CMO(:), &
                       DIAF=DIAF(:), &
                       D1I_AO=D1I(:), &
                       D1A_AO=D1A(:), &
                       TUVX=tuvx(:), &
                       F_IN=FI(:), &
                       D1S_MO=DSPN(:), &
                       DMAT=DMAT(:), &
                       PSMAT=pmat(:), &
                       PAMAT=pa(:))
# if defined (_ENABLE_BLOCK_DMRG_) || defined (_ENABLE_CHEMPS2_DMRG_) || defined (_ENABLE_DICE_SHCI_)
  else if (DoBlockDMRG) then
    call DMRGCTL(CMO,DMAT,DSPN,PMAT,PA,FI,D1I,D1A,TUVX,IFINAL,1)
# endif
  else
    call CICTL(CMO,DMAT,DSPN,PMAT,PA,FI,FA,D1I,D1A,TUVX,IFINAL)
  end if

  ! call triprt('twxy',' ',TUVX,nAc*(nAc+1)/2)
  ! call triprt('P-mat 2',' ',PMAT,nAc*(nAc+1)/2)

  EAV = Zero
  if (DoSplitCAS) then
    EAV = ENER(lRootSplit,ITER)
  else
    do KROOT=1,NROOTS
      EAV = EAV+ENER(IROOT(KROOT),ITER)*WEIGHT(KROOT)
    end do
  end if

  if (IPRLEV >= DEBUG) then
    write(u6,*) 'EAV value in RASSCF after second call to CICTL:'
    write(u6,*) EAV

    write(u6,*) 'Printing matrices in RASSCF'
    write(u6,*)
    write(u6,*) ' CMO in RASSCF'
    write(u6,*) ' ---------------------'
    write(u6,*)
    ioff = 1
    do iSym=1,nSym
      iBas = nBas(iSym)
      if (iBas /= 0) then
        write(u6,*) 'Sym =',iSym
        do i=1,iBas
          write(u6,*) (CMO(ioff+iBas*(i-1)+j),j=0,iBas-1)
        end do
        iOff = iOff+(iBas*iBas)
      end if
    end do

    write(u6,*)
    write(u6,*) ' D1I in AO basis in RASSCF'
    write(u6,*) ' ---------------------'
    write(u6,*)
    iOff = 1
    do iSym=1,nSym
      iBas = nBas(iSym)
      call wrtmat(D1I(ioff),iBas,iBas,iBas,iBas)
      iOff = iOff+iBas*iBas
    end do
    write(u6,*)
    write(u6,*) 'Total Charge :',Tot_Charge

    call mma_allocate(Tmp1,nTot1,Label='Tmp1')
    iComp = 1
    iSyLbl = 1
    iRc = -1
    iOpt = ibset(ibset(0,sNoOri),sNoNuc)
    Label = 'OneHam'
    call RdOne(iRc,iOpt,Label,iComp,Tmp1,iSyLbl)
    if (iRc /= 0) then
      write(u6,*) 'SGFCIN: iRc from Call RdOne not 0'
#     ifdef _FDE_
      write(u6,*) 'Label = ',Label
#     endif
      write(u6,*) 'iRc = ',iRc
      call Abend()
    end if

    write(u6,*)
    write(u6,*) ' OneHam in AO basis in RASSCF'
    write(u6,*) ' ---------------------'
    write(u6,*)
    iOff = 1
    do iSym=1,nSym
      iBas = nBas(iSym)
      call TriPrt(' ','(5G17.11)',Tmp1(iOff),iBas)
      iOff = iOff+(iBas*iBas+iBas)/2
    end do

    call mma_deallocate(Tmp1)
    call Get_dScalar('PotNuc',potNuc)

    write(u6,*)
    write(u6,*) 'PotNuc :',PotNuc

    write(u6,*)
    write(u6,*) ' D1A in AO basis in RASSCF'
    write(u6,*) ' ---------------------'
    write(u6,*)
    iOff = 1
    do iSym=1,nSym
      iBas = nBas(iSym)
      call wrtmat(D1A(ioff),iBas,iBas,iBas,iBas)
      iOff = iOff+iBas*iBas
    end do
  end if
else
  call mma_allocate(FOCK,NACPAR,Label='Fock')
  ! To fix the DS bug... I forgot to transform it to the AO basis... Agrrrrhhh!
  if (iSpin == 1) then
    if (IPRLEV >= DEBUG) write(u6,*) 'running a singlet. DSPN set to zero!'
    DSPN(:) = Zero
  end if
  call mma_allocate(TmpDS,NACPAR,Label='TmpDS')
  call mma_allocate(TmpD1S,NTOT2,Label='TmpD1S')
  call dcopy_(NACPAR,DSPN,1,TmpDS,1)
  if (NASH(1) /= NAC) call DBLOCK(TmpDS)
  call Get_D1A_RASSCF(CMO,TmpDS,TmpD1S)
  call mma_deallocate(TmpDS)
  call CASDFT_terms(CMO,FOCK,FI,D1I,D1A,TmpD1S)
  call mma_deallocate(TmpD1S)
  call mma_deallocate(FOCK)
end if

!call TRIPRT('Averaged one-body density matrix, D, in RASSCF',' ',DMAT,NAC)
!call TRIPRT('Averaged one-body spin density matrix DS, RASSCF',' ',DSPN,NAC)
!call TRIPRT('Averaged two-body density matrix, P',' ',PMAT,NACPAR)
!call TRIPRT('Averaged antisym 2-body density matrix PA RASSCF',' ',PA,NACPAR)

call Get_D1A_RASSCF(CMO,DMAT,D1A)
if (IPRLEV >= DEBUG) then
  write(u6,*)
  write(u6,*) ' D1A in AO basis in RASSCF bf SXCTL'
  write(u6,*) ' ---------------------'
  write(u6,*)
  iOff = 1
  do iSym=1,nSym
    iBas = nBas(iSym)
    call wrtmat(D1A(ioff),iBas,iBas,iBas,iBas)
    iOff = iOff+iBas*iBas
  end do
end if
call Timing(dum1,dum2,time3(2),dum3)
TimeCIOpt = TimeCIOpt+time3(2)-time3(1)

!call rasscf_xml(Iter)
call rasscf_mcontrol(Iter)

! SX-section

call Timing(dum1,dum2,time2(1),dum3)

if (IPRLEV >= DEBUG) then
  write(u6,*) ' In RASSCF bf SXCTL'
  call TRIPRT('Averaged one-body density matrix, D, in RASSCF',' ',DMAT,NAC)
  call TRIPRT('Averaged two-body density matrix, P',' ',PMAT,NACPAR)
  call TRIPRT('Averaged antisym 2-body density matrix PA RASSCF',' ',PA,NACPAR)
end if
call SXCTL(CMO,OCCN,DMAT,PMAT,PA,FI,FA,D1A,THMAX,IFINAL)

if (IPRLEV >= DEBUG) then
  write(u6,*)
  write(u6,*) ' FI+FA in RASSCF after SXCTL'
  write(u6,*) ' ---------------------'
  write(u6,*)
  iOff = 1
  do iSym=1,nSym
    iBas = nBas(iSym)
    call TriPrt(' ',' ',FA(iOff),iBas)
    iOff = iOff+(iBas*iBas+iBas)/2
  end do
end if

!GLM write(u6,*) 'ECAS in RASSCF after call to SXCTL',ECAS
if (IPRLEV >= DEBUG) then
  write(u6,*)
  write(u6,*) ' D1A in AO basis in RASSCF af SXCTL'
  write(u6,*) ' ---------------------'
  write(u6,*)
  iOff = 1
  do iSym=1,nSym
    iBas = nBas(iSym)
    call wrtmat(D1A(ioff),iBas,iBas,iBas,iBas)
    iOff = iOff+iBas*iBas
  end do
end if
call Get_D1A_RASSCF(CMO,DMAT,D1A)
CASDFT_Funct = Zero
if ((KSDFT /= 'SCF') .and. (KSDFT /= 'PAM')) call Get_dScalar('CASDFT energy',CASDFT_Funct)

DE = (ECAS+CASDFT_Funct)-ECAS1
ECAS1 = ECAS+CASDFT_Funct

if (NAC == 0) EAV = ECAS

if ((KSDFT(1:3) /= 'SCF') .and. (KSDFT(1:3) /= 'PAM')) then
  if (nConf == 1) then
    EAV = ECAS
  else
    EAV = EAV-VIA_DFT-HALFQ
  end if
end if

if (ITER == 1) DE = Zero
call Timing(dum1,dum2,time2(2),dum3)
TimeOrb = TimeOrb+time2(2)-time2(1)
TMXTOT = max(TMXTOT,THMAX)

! Save energies and convergence parameters

CONV(1,ITER) = ECAS
CONV(2,ITER) = ESX
CONV(3,ITER) = CMAX
CONV(4,ITER) = DE
CONV(5,ITER) = CBLBM
CONV(6,ITER) = ROTMAX
IAD15 = IADR15(6)
call DDAFILE(JOBIPH,1,ENER,mxRoot*mxIter,IAD15)
call DDAFILE(JOBIPH,1,CONV,6*mxIter,IAD15)
#ifdef _HDF5_
call mh5_put_attr(wfn_iter,Iter)
call mh5_put_dset(wfn_energy,ENER(1,Iter))
#endif

! Print output of energies and convergence parameters

call Timing(dum1,dum2,time0(2),dum3)
time0(2) = time0(2)-time0(1)
! Character indicating unconvergence/convergence criterion fulfilled:
CTHRE = ' '
CTHRSX = ' '
CTHRTE = ' '
if (abs(DE) > THRE) CTHRE = '*'
if (abs(CBLBM) > THRSX) CTHRSX = '*'
if (abs(ROTMAX) > THRTE) CTHRTE = '*'
IROT = 0
if (NROOTS == 1) IROT = IROOT(1)
if ((IPRLEV >= 2) .and. (IPRLEV <= 3) .and. (.not. lOPTO)) then
  !----------------------------------
  ! Shift total energies (BOR 070411)
  if ((iter == 1) .and. (ICIONLY == 0)) then
    kau = int(ECAS*1.0e-3_wp)
    EVAC = 1.0e3_wp*real(kau,kind=wp)
    if (kau /= 0) then
      write(u6,'(6x,A,f23.2,A)') 'Total energies have been shifted. Add ',EVAC,' au'
    end if
  end if
  !----------------------------------
  ihh = int(time0(2)/3600)
  imm = int(time0(2)-ihh*3600)/60
  iss = int(time0(2)-ihh*3600-imm*60)
  if (DoSplitCAS) then
    write(u6,101) ITER,iterSplit,ITERSX,IROT,EAV,DE,CTHRE,ROTMAX,CTHRTE,IBLBM,JBLBM,ISYMBB,CBLBM,CTHRSX,SXSHFT,TMIN,QNSTEP,QNUPDT, &
                  ihh,':',imm,':',iss
  else if (DoBKAP) then
    write(u6,102) ITER,ITERCI,ITERSX,IROT,ECAS-EVAC+CASDFT_Funct,EAV,DE,CTHRE,ROTMAX,CTHRTE,IBLBM,JBLBM,ISYMBB,CBLBM,CTHRSX, &
                  SXSHFT,TMIN,QNSTEP,QNUPDT,ihh,':',imm,':',iss
  else
    if (doDMRG .and. KeyCION) then ! If DMRG only -- yma
      write(u6,'(/6X,a,F19.10)') 'DMRGCI energy              :',EAV+CASDFT_Funct
      write(u6,'(6X,a,I5,A1,I2.2,A1,I2.2/)') 'Total time spent (hh:mm:ss):        ',ihh,':',imm,':',iss
    else
      if (doDMRG) then
#       ifdef _DMRG_
        maxtrW = Zero
        maxtrR = -1
        maxBD = -1
        ! These dmrg variables are arrays of rank 1
        ITERCI = maxval(dmrg_energy%num_sweeps)
        IROT = maxloc(dmrg_energy%num_sweeps,1)
        maxtrW = maxval(dmrg_energy%max_truncW)
        maxtrR = maxloc(dmrg_energy%max_truncW,1)
        maxBD = maxval(dmrg_energy%bond_dim)
        write(u6,103) ITER,ITERCI,IROT,maxBD,maxtrW,maxtrR,ITERSX,ECAS-EVAC+CASDFT_Funct,DE,CTHRE,ROTMAX,CTHRTE,IBLBM,JBLBM, &
                      ISYMBB,CBLBM,CTHRSX,SXSHFT,TMIN,QNSTEP,QNUPDT,ihh,':',imm,':',iss
#       endif
      else
        write(u6,104) ITER,ITERCI,ITERSX,IROT,ECAS-EVAC+CASDFT_Funct,DE,CTHRE,ROTMAX,CTHRTE,IBLBM,JBLBM,ISYMBB,CBLBM,CTHRSX, &
                      SXSHFT,TMIN,QNSTEP,QNUPDT,ihh,':',imm,':',iss
      end if
    end if
  end if
else if (IPRLEV >= 4) then
  write(u6,'(6X,A,I4)') 'Energy statistics and convergence in iteration',ITER
  write(u6,'(6X,A,10X,F26.6)') 'Average of CI energies',EAV-EVAC
  write(u6,'(6X,A,F26.6)') 'Complete active space SCF energy',ECAS-EVAC
  write(u6,'(6X,A,17X,F26.6)') 'Super-CI energy',ESX
  write(u6,'(6X,A,12X,F26.6,3X,A1)') 'RASSCF energy change',DE,CTHRE
  write(u6,'(6X,A,I1,A1,2I3,F19.6)') 'Maximum BLB matrix element(sym=',ISYMBB,')',IBLBM,JBLBM,CBLBM
  write(u6,'(6X,A,10X,F26.6,3X,A1)') 'Max rotation parameter',ROTMAX,CTHRSX
  write(u6,'(6X,A,14X,F26.6,3X,A1)') 'Max rotation angle',THMAX,CTHRTE
  write(u6,'(6X,A,3X,F26.6)') 'Max change in MO coefficients',CMAX
end if
#ifdef _FDE_
! Embedding
if (embpot) then
  Eemb = embPotEneMODensities(D1I,D1A,embInt,nBas,nTot2,nSym)
  write(u6,*) 'E from embedding potential (<Psi|v_emb|Psi>): ',Eemb
end if
#endif
!GLM some additional printout for MC-PDFT

if (l_casdft) then
  CASDFT_E = ECAS-EVAC+CASDFT_Funct
  call print_mcpdft(CASDFT_E)
end if

! Compare RASSCF energy and average CI energy

DIFFE = abs((ECAS-EAV)/ECAS)
if (.not. DoSplitCAS) then
  if (DoBKAP) then
    if (iter == 1) then
      if ((DIFFE > 1.0e-10_wp) .and. (NROOTS == 1)) then
        write(u6,'(6X,A)') repeat('=',120)
        call WarningMessage(2,'Rasscf and CI energies will differ.')
        write(u6,*) 'This is the price you pay by the diagonal approximation over the BB block in the SplitCAS method.'
        write(u6,*) 'The RASSCF energy might also diverge!'
        write(u6,'(A)') repeat('#',80)
      end if
    end if
  else if (l_casdft) then
    write(u6,'(6X,A)') repeat('=',80)
    write(u6,'(10X,A)') 'This is a POST-SCF correction using a modified  Hamiltonian.'
    write(u6,'(10X,A)') 'The RASSCF energy has been corrected and it will differ from'
    write(u6,'(10X,A)') 'the preceding CI energy.'
    write(u6,'(6X,A)') repeat('=',80)
  else
    if (doDMRG) then

#     ifdef _DMRG_DEBUGPRINT_
      write(u6,*) 'DMRG-SCF energy    ',ECAS
      write(u6,*) 'DMRG sweeped energy',EAV
#     endif

      if (.not. KeyCION) then
        if ((DIFFE > 1.0e-6_wp) .and. (NROOTS == 1)) then
          write(u6,'(6X,A)') repeat('=',120)
          call WarningMessage(2,'DMRGSCF and DMRG energies differ.')
          write(u6,'(6X,A,I11)') 'iteration           ',ITER
          write(u6,'(6X,A,F22.10)') 'DMRGSCF energy      ',ECAS
          write(u6,'(6X,A,F22.10)') 'DMRG energy         ',EAV
          write(u6,'(6X,A,F22.10)') 'relative difference ',DIFFE
          write(u6,*) 'About this difference:'
          write(u6,*) '1) If possible, consider a larger M value'
          write(u6,*) '2) Severe convergence problems. Maybe active'
          write(u6,*) '   space is unsuitable for this system?'
          write(u6,'(6X,A)') repeat('=',120)
          if ((DIFFE > 5.0e-4_wp) .and. (NROOTS == 1)) then
            write(u6,*)
            write(u6,*) 'Warning : '
            write(u6,*) ' Relative difference is near unacceptable'
            write(u6,*) ' If possible, consider a larger M value'
            write(u6,*)
          end if
        end if
      end if
    else
      DIFFETol = 1.0e-10_wp
#     ifdef _ENABLE_DICE_SHCI_
      if (DoBlockDMRG) DIFFETol = 1.0e-8_wp
#     endif
      if ((DIFFE > DIFFETol) .and. (NROOTS == 1)) then
        write(u6,'(6X,A)') repeat('=',120)
        call WarningMessage(2,'Rasscf and CI energies differ.')
        write(u6,'(6X,A,I11)') 'iteration           ',ITER
        write(u6,'(6X,A,F22.10)') 'RASSCF energy       ',ECAS
        write(u6,'(6X,A,F22.10)') 'CI energy           ',EAV
        write(u6,'(6X,A,F22.10)') 'relative difference ',DIFFE
        write(u6,*) 'Severe convergence problems. Maybe the active'
        write(u6,*) '   space is unsuitable for this system?'
        write(u6,'(6X,A)') repeat('=',120)
        if ((DIFFE > 1.0e-4_wp) .and. (NROOTS == 1) .and. (.not. l_casdft)) then
          write(u6,*)
          write(u6,'(6X,A)') 'The program has to stop !!!'
          write(u6,*)
          write(u6,'(6X,A)') repeat('=',120)
          write(u6,*)
          ITERM = 99
          goto 2000
        end if
      end if
    end if
  end if
else
  !write(u6,'(6X,A,F22.10)') 'Split-RASSCF energy    ',ECAS
  if (DIFFE > 5.0e-3_wp) then
    write(u6,'(6X,A)') repeat('*',120)
    write(u6,'(6X,A)') 'The Split-RASSCF and Split-CI energies differ !!!'
    !write(u6,'(6X,A,I11)') 'iteration           ',ITER
    write(u6,'(6X,A,F22.10)') 'Split-RASSCF energy    ',ECAS
    write(u6,'(6X,A,F22.10)') 'Split-CI energy        ',EAV
    write(u6,'(6X,A,F22.10)') 'Relative difference    ',DIFFE
    write(u6,*) '     Smaller is the difference more realiable is the result.'
    write(u6,*) '     To make the difference smaller try to select a bigger AA block or use firstOrder keyword.'
    write(u6,'(6X,A)') repeat('*',120)
  end if
end if

if (write_orb_per_iter .and. king()) then
  call copy_(real_path('RASORB'),real_path('ITERORB.'//str(actual_iter)))
# ifdef _HDF5_
  call copy_(real_path('RASWFN'),real_path('RASWFN.'//str(actual_iter)))
# endif
end if

! Convergence check:
! check is done on largest BLB matrix
! element (CBLBM), on difference in
! average energy, DE and on maximum value of non-
! diagonal rotation matrix element.
!
!***********************************************************************
!***********************************************************************
! IF CIONLY calculation the convergence is skipped and goes to line 2000
!***********************************************************************
!***********************************************************************

if (IFINAL == 1) goto 2000
if (DE > One) then
  call StatusLine('RASSCF: ','No convergence.')
  write(u6,*)
  write(u6,'(6X,A)') repeat('=',120)
  call WarningMessage(2,'Rasscf energy diverges.')
  write(u6,'(6X,A,I11)') 'iteration           ',ITER
  write(u6,'(6X,A,F22.10)') 'RASSCF energy       ',ECAS
  write(u6,'(6X,A,F22.10)') 'energy difference   ',DE
  write(u6,'(6X,A)') repeat('=',120)
  write(u6,*)
  write(u6,'(6X,A)') '!!! The program was forced to stop !!!'
  write(u6,*)
  write(u6,'(6X,A)') repeat('=',120)
  write(u6,*)
  ITERM = 99
  goto 2000
end if
if (ITER < MAXIT) then
  write(STLNE2,'(A12,I3)') 'Iteration ',ITER
  call StatusLine('RASSCF converged: ',STLNE2)
  if (abs(DE) > THRE) GO TO 1000
  if (abs(CBLBM) > THRSX) GO TO 1000
  if (abs(ROTMAX) > THRTE) GO TO 1000
  if ((ITER <= 3) .and. (ICIONLY == 0)) GO TO 1000   ! 3->0 checking
  if (IPRLEV >= TERSE) write(u6,'(6X,A,I3,A)') 'Convergence after',ITER,' iterations'
  if (IPRLEV >= DEBUG) then
    write(u6,*)
    write(u6,*) ' CMO in RASSCF after convergence printout'
    write(u6,*) ' ---------------------'
    write(u6,*)
    ioff = 1
    do iSym=1,nSym
      iBas = nBas(iSym)
      if (iBas /= 0) then
        write(u6,*) 'Sym =',iSym
        do i=1,iBas
          write(u6,*) (CMO(ioff+iBas*(i-1)+j),j=0,iBas-1)
        end do
        iOff = iOff+(iBas*iBas)
      end if
    end do
  end if
  IFINAL = 1
  call Add_Info('RASSCF_ITER',[real(ITER,kind=wp)],1,8)
  !call Add_Info('RASSCF_THMX',TMXTOT,1,5)
  goto 1000
else
  if (IPRLEV >= TERSE) write(u6,'(6X,A,I3,A)') 'No convergence after',ITER,' iterations'
  write(STLNE2,'(A12,I3)') 'Iteration ',ITER
  call StatusLine('RASSCF max iter: ',STLNE2)
  IFINAL = 1
  ITERM = 16
  goto 1000
end if
!                                                                      *
!***********************************************************************
!                                                                      *
!     Compute Final CI vectors
!                                                                      *
!***********************************************************************
!                                                                      *
2000 IFINAL = 2
ICICH = 0

if (KeyWRMA) then
  call dump_fciqmc_mats(dmat=DMAT(:),dspn=DSPN(:),psmat=pmat(:),pamat=pa(:))
end if

!***********************************************************************
!*****************           Closing up MC-PDFT      *******************
!***********************************************************************

! Clean-close as much as you can the CASDFT stuff...
if (l_casdft) goto 2010

!* IPT2 = 1 for OUTO, CANOnical keyword...
if (IPT2 == 1) then
  IAD = IADR15(9)
  call DDAFILE(JOBIPH,2,CMO,NTOT2,IAD)
else
  IAD = IADR15(2)
  call DDAFILE(JOBIPH,2,CMO,NTOT2,IAD)
end if
if (IPRLEV >= DEBUG) then
  write(u6,*)
  write(u6,*) ' CMO in RASSCF after DDAFILE'
  write(u6,*) ' ---------------------'
  write(u6,*)
  ioff = 1
  do iSym=1,nSym
    iBas = nBas(iSym)
    if (iBas /= 0) then
      write(u6,*) 'Sym =',iSym
      do i=1,iBas
        write(u6,*) (CMO(ioff+iBas*(i-1)+j),j=0,iBas-1)
      end do
      iOff = iOff+(iBas*iBas)
    end if
  end do
end if
if (NROOTS > 1) then
  call StatusLine('RASSCF: ','Compute final CI vectors')
else
  call StatusLine('RASSCF: ','Compute final CI vector')
end if

! Transform two-electron integrals

call Timing(dum1,dum2,time2(1),dum3)
if ((.not. DoCholesky) .or. (ALGO == 1)) then
  call mma_allocate(PUVX,NFINT,Label='PUVX')
  PUVX(:) = Zero
end if

call Get_D1I_RASSCF(CMO,D1I)

IPR = 0
if (IPRLOC(2) == 4) IPR = 5
if (IPRLOC(2) == 5) IPR = 10
call TRACTL2(CMO,PUVX,TUVX,D1I,FI,D1A,FA,IPR,lSquare,ExFac)

if ((.not. DoCholesky) .or. (ALGO == 1)) call mma_deallocate(PUVX)

call Timing(dum1,dum2,time2(2),dum3)
TimeTrans = TimeTrans+time2(2)-time2(1)

! CI-section (to obtain final wave function)
! 1st and 2nd order density matrix in MO basis
! for the gradients.

! Print header to file containing informations on CI iterations.

write(IterFile,*)
write(IterFile,'(A)') repeat('*',80)
write(IterFile,'(15X,A)') 'RASSCF iteration: Final'

call Timing(dum1,dum2,time3(1),dum3)

if (allocated(CI_solver)) then
  call CI_solver%run(actual_iter=actual_iter, &
                     ifinal=ifinal, &
                     iroot=iroot, &
                     weight=weight, &
                     CMO=CMO(:), &
                     DIAF=DIAF(:), &
                     D1I_AO=D1I(:), &
                     D1A_AO=D1A(:), &
                     TUVX=tuvx(:), &
                     F_IN=FI(:), &
                     D1S_MO=DSPN(:), &
                     DMAT=DMAT(:), &
                     PSMAT=pmat(:), &
                     PAMAT=pa(:))

#if defined (_ENABLE_BLOCK_DMRG_) || defined (_ENABLE_CHEMPS2_DMRG_) || defined (_ENABLE_DICE_SHCI_)
else if (DoBlockDMRG) then
  call DMRGCTL(CMO,DMAT,DSPN,PMAT,PA,FI,D1I,D1A,TUVX,IFINAL,1)
#endif
#ifdef _DMRG_
! Leon 27/11/2017: Skip the final CI iteration if we're using DMRGCI
! and CIOnly. It's enabled only for DMRGCI with QCMaquis now
! (to exclude potential side effects)
! but consider extending it to other cases!
else if (doDMRG .and. (ICIONLY /= 0)) then
  continue
#endif
else
  call CICTL(CMO,DMAT,DSPN,PMAT,PA,FI,FA,D1I,D1A,TUVX,IFINAL)
end if
if (lRF .and. ((iPCMRoot <= 0) .or. (DWSolv%DWZeta /= Zero))) then
  IAD15 = IADR15(6)
  call DDAFILE(JOBIPH,1,ENER,mxRoot*mxIter,IAD15)
end if

EAV = Zero
if (DoSplitCAS) then
  EAV = ENER(lRootSplit,ITER)
else
  do KROOT=1,NROOTS
    EAV = EAV+ENER(IROOT(KROOT),ITER)*WEIGHT(KROOT)
  end do
end if
if (NAC == 0) EAV = ECAS
if (NCRVEC > 0) then
! Core shift has been used
  call mma_deallocate(CRVEC)
  call mma_deallocate(CRPROJ)
end if
call Timing(dum1,dum2,time3(2),dum3)
TimeCIOpt = TimeCIOpt+time3(2)-time3(1)
TimeWfn = TimeWfn+time3(2)-time1(1)
time1(1) = time3(2)

! Calculation of natural orbitals. These orbitals are stored on
! JOBIPH in IADR15(12), followed by the occupation numbers.
! Usage: Only for one-electron properties
! Note: in an average calculation natural orbitals are obtained
! for all roots. Stored sequentially at IADR15(12) on JOBPIH.
!PAM2009 This is done in NATORB. Note that it places the resulting
! natural orbitals and occ nos into JOBIPH. Thus the arguments nr
! 5 and 6 are simply those for the last root that NATORB processed,
! and they should not be reused for anything after return from NATORB.
! NATORB args: Arg1 is current CMO coeffs, used in CI;
!  all the rest should be regarded as scratch.

call StatusLine('RASSCF: ','Compute natural orbitals')
IPR = 0
if (IPRLOC(6) == 4) IPR = 5
if (IPRLOC(6) == 5) IPR = 10
call mma_allocate(CMON,NTOT2,Label='CMON')
call mma_allocate(OCCX,NTOT,Label='OCCX')
call mma_allocate(scr1,max(NACPAR,NO2M),Label='scr1')
call mma_allocate(scr2,NO2M,Label='scr2')
call mma_allocate(SMAT,NTOT1,Label='SMAT')
!PAM2009 NATORB args: Arg1 is current CMO coeffs, used in CI;
!  all the rest should be regarded as scratch.
call NATORB_RASSCF(CMO,scr1,scr2,SMAT,CMON,OCCX)
call mma_deallocate(scr1)
call mma_deallocate(scr2)
!PAM2009 Deallocate CMON, OCCX.
call mma_deallocate(CMON)
call mma_deallocate(OCCX)

! Compute transition density matrices
if (KeyTDM) then
# ifdef _HDF5_
  call mma_allocate(Tmp,NConf,Label='Tmp')
  call mma_allocate(VecL,NConf,Label='VecL')
  call mma_allocate(VecR,NConf,Label='VecR')
  call mma_allocate(kcnf,NACTEL,Label='kcnf')
  call mma_allocate(Dtmp,NAC*NAC,Label='Dtmp')
  call mma_allocate(DStmp,NAC*NAC,Label='DStmp')
  jDisk = IADR15(4)
  call DDafile(JOBIPH,2,Tmp,nConf,jDisk)
  do jRoot=2,lRoots
    ! Read and reorder the left CI vector
    call DDafile(JOBIPH,2,Tmp,nConf,jDisk)
    call Reord2(NAC,NACTEL,STSYM,1,CONF,CFTP,Tmp,VecL,kcnf)
    kDisk = IADR15(4)
    do kRoot=1,jRoot-1
      ! Read and reorder the right CI vector
      call DDafile(JOBIPH,2,Tmp,nConf,kDisk)
      call Reord2(NAC,NACTEL,STSYM,1,CONF,CFTP,Tmp,VecR,kcnf)
      ! Compute TDM and store in h5 file
      call Lucia_Util('Densi',CI_Vector=VecL(:),RVec=VecR(:))
      idx = (jRoot-2)*(jRoot-1)/2+kRoot
      call mh5_put_dset(wfn_transdens,Dtmp(1:NAC*NAC),[NAC,NAC,1],[0,0,idx-1])
      if (iSpin > 1) call mh5_put_dset(wfn_transsdens,DStmp(1:NAC**2),[NAC,NAC,1],[0,0,idx-1])
    end do
  end do
  call mma_deallocate(TMP)
  call mma_deallocate(VecL)
  call mma_deallocate(VecR)
  call mma_deallocate(kcnf)
  call mma_deallocate(Dtmp)
  call mma_deallocate(DStmp)
# else
  call WarningMessage(1,'HDF5 support disabled, TDM keyword ignored.')
# endif
end if

if (KeySSCR) then
  call spin_correlation_driver(orb_range_p,orb_range_q,iroot)
  call mma_deallocate(orb_range_p)
  call mma_deallocate(orb_range_q)
end if

!***********************************************************************
! Export all information relevant to geometry optimizations.
! Save also the reaction field operator.
call Timing(dum1,dum2,time2(1),dum3)
if (iRlxRoot == 0) iRlxRoot = iRoot(1)

! Replace average occ Fock with occ Fock for state iRlxRoot
! and densities with the densities for state iRLXRoot
!write(u6,*) 'I am in RASSCF before call to PutRlx!'
if (ITERM /= 99) then
  call mma_allocate(Dens,nTot1,Label='Dens')
  call PutRlx(DMAT,DSPN,PMAT,Dens,CMO)
  call Export1(IFINAL,CMO,DMAT,PMAT,Dens,FockOcc)
  call mma_deallocate(Dens)
end if
call Timing(dum1,dum2,time2(2),dum3)
TimeRelax = time2(2)-time2(1)
!***********************************************************************

EMY = EMY+CASDFT_Funct
!                                                                      *
!***********************************************************************
!                                                                      *
! Output section
!                                                                      *
!***********************************************************************
!                                                                      *

call StatusLine('RASSCF: ','Printing results')
if ((IPRLEV >= USUAL) .and. (.not. lOPTO)) then
  write(u6,*)
  write(u6,'(6X,A)') repeat('*',120)
  write(u6,'(6X,A,118X,A)') '*','*'
  write(u6,'(6X,A,52X,A,53X,A)') '*','Final results','*'
  write(u6,'(6X,A,118X,A)') '*','*'
  write(u6,'(6X,A)') repeat('*',120)
  write(u6,*)
end if
#ifdef _FDE_
! Embedding
if (embpot) then
  Eemb = embPotEneMODensities(D1I,D1A,embInt,nBas,nTot2,nSym)
  write(u6,*) 'Final energy from embedding potential: ',Eemb
  ! Write out ESP on grid if requested
  if (embWriteEsp) then
    call Get_iScalar('Unique atoms',nNuc)
    call embPotOutputMODensities(nNuc,nSym,D1I,D1A,nBas,nTot2)
  end if
end if
#endif

if (ITERM /= 99) then
  if (.not. DoSplitCAS) then
    call OUTCTL(CMO,OCCN,SMAT,lOPTO)
  else
    call OUTCTLSplit(CMO,OCCN,SMAT,lOPTO)
  end if
end if

call mma_deallocate(SMAT)

! Write information for MOLDEN

!  i_root=0 gives state- and spin-averaged natural orbitals
!  i_root>0 gives natural spin orbitals for that root
mroots = nroots
if ((iSpin == 1) .and. (nroots == 1)) mroots = 0
do i_root=0,mroots
  call Interf(i_root,FDIAG,1,0)
end do

! Create output orbital files:
call OrbFiles(JOBIPH,IPRLEV)

!***********************************************************************
!*****************           Closing up RASSCF       *******************
!***********************************************************************

2010 continue

if (DoCholesky .and. (ALGO == 2)) call mma_deallocate(Qmat)

! Release  some memory allocations
call mma_deallocate(DIAF)
call mma_deallocate(FockOcc)
call mma_deallocate(FI)
call mma_deallocate(FA)
call mma_deallocate(D1I)
call mma_deallocate(D1A)
call mma_deallocate(OccN)
call mma_deallocate(CMO)
#ifdef _DMRG_
! Free RDMs for the reaction field reference root in QCMaquis calculations
if (doDMRG .and. PCM_On()) then
  call mma_deallocate(RF1)
  if (twordm_qcm) call mma_deallocate(RF2)
end if
#endif
if (lRF) call DWSol_final()
!call DWSCF_final()

! deallocating TUVX memory...
call mma_deallocate(TUVX)
call mma_deallocate(DSPN)
call mma_deallocate(DMAT)
call mma_deallocate(PMAT)
call mma_deallocate(PA)
!Leon: The velociraptor comes! xkcd.com/292/
9989 continue

! release SEWARD
call ClsSew()
! ClsSew is needed for releasing memory used by integral_util, rys... which is allocated when MC-PDFT run is performed.

! Finalize Cholesky information if initialized
if (DoCholesky) then
  call Cho_X_Final(irc)
  if (irc /= 0) then
    write(u6,*) 'RASSCF: Cho_X_Final fails with return code ',irc
    write(u6,*) ' Try to recover. Calculation continues.'
  end if
  if (Do_OFemb) then
    call mma_deallocate(FMaux)
    call OFE_print(EAV)
  end if
end if

!do i=1,NTOT2
!  write(u6,*) 'A,I',D1A(i),D1I(i)
!end do
call mma_deallocate(CleanMask,safe='*')

! Skip Lucia stuff if NECI or BLOCK-DMRG is on
if (.not. any([allocated(CI_solver),DumpOnly,doDMRG,doBlockDMRG])) call Lucia_Util('CLOSE')

call StatusLine('RASSCF: ','Finished.')
if (IPRLEV >= 2) write(u6,*)
if (ifvb == 1) call make_close_rvb()
!vv call to grid is moved up, in order to call clssew safely..
!if (iCIonly == 0) call Grid_driver(-1,'RASSCF','RASORB',iR)

call Timing(dum1,dum2,time1(2),dum3)
TimeTotal = time1(2)
TimeOutput = TimeOutput+time1(2)-time1(1)
if (IPRLEV >= 3) then
  call PrtTim()
  call FastIO('STATUS')
end if
call ClsFls_RASSCF()

! Rc_RAS  =  0 : The RASSCF wave function is converged
!         = 16 : The RASSCF wave function is not(!) converged
!         = 99 : The RASSCF energy is divergent or
!                the CI and SX energies differ
Rc_RAS = ITERM
Rc_RAS = max(RC_RAS,Rc_CI)
Rc_RAS = max(RC_RAS,Rc_SX)
if (Rc_Ras == 0) then
  ireturn = _RC_ALL_IS_WELL_
else if (Rc_Ras == 16) then
  ireturn = _RC_NOT_CONVERGED_
else
  call WarningMessage(2,'Something is wrong: Did CI fail?')
  ireturn = _RC_GENERAL_ERROR_
end if

if (Do_OFemb) then
  if (isStructure() == 1) then
    if (iReturn /= _RC_ALL_IS_WELL_) call WarningMessage(1,'RASSCF: non-zero return code.')
    iReturn = _RC_CONTINUE_LOOP_
    call Check_FThaw(iReturn)
  end if
end if

if (.not. (iDoGas .or. doDMRG .or. doBlockDMRG .or. allocated(CI_solver) .or. DumpOnly)) call MKGUGA_FREE(SGS,CIS,EXS)

if (DoFaro) then
  call faroald_free()
  call citrans_free()
end if

if (allocated(CI_solver)) then
  call CI_solver%cleanup()
  deallocate(CI_solver)
end if

! DMRG: Save results for other use
! ==========================================================
if (doDMRG) then
# ifdef _DMRG_
  !Leon: Generate 4-RDM evaluation templates for NEVPT2

  ! In the new interface the EvRDM keyword will be ignored.
  ! Instead, NEVPT2Prep will always generate the template.
  ! RDM evaluation will now happen in the NEVPT2 module
  ! where NEVPT2 either can attempt to compute it directly
  ! with the new interface or read from QCMaquis HDF5 result
  ! file.
  if (DoNEVPT2Prep) then
    if (MAXIT == 0) write(u6,*) ' --- DMRG-SCF iterations are skipped, only QCMaquis input for higher-order RDMs will be generated.'
    if (NACTEL > 3) then ! Ignore 4-RDM if we have <4 electrons
      do i=1,NROOTS
        write(u6,'(a)') 'Writing 4-RDM QCMaquis template for state '//str(i)
        call qcmaquis_interface_prepare_hirdm_template(filename='meas-4rdm.'//str(i-1)//'.in',state=i-1,tpl=TEMPLATE_4RDM)
        call qcmaquis_mpssi_transform(trim(qcmaquis_param%workdir)//'/'//trim(qcmaquis_param%project_name),i)
      end do
    else
      write(u6,*) 'Skipping 4-RDM QCMaquis template generation since we have less than 4 electrons.'
    end if
    ! Generate 3-TDM templates
    if (NACTEL > 2) then ! but only if we have more than 3 el.
      do i=1,NROOTS
        do j=i+1,NROOTS
          write(u6,'(a)') 'Writing 3-TDM QCMaquis template for states '//str(i)//' and '//str(j)
          call qcmaquis_interface_prepare_hirdm_template(filename='meas-3tdm.'//str(i-1)//'.'//str(j-1)//'.in', &
                                                         state=i-1, &
                                                         state_j=j-1, &
                                                         tpl=TEMPLATE_TRANSITION_3RDM)
        end do
      end do
    else
      write(u6,*) 'Skipping 3-RDM QCMaquis template generation since we have less than 3 electrons.'
    end if
  end if
  ! is it really needed in the times of Fortran 2008?
  call qcmaquis_interface_deinit()
# endif
end if
! ==========================================================
! Exit

9990 continue
! Close the one-electron integral file:
iRC = -1
iOpt = 0
call ClsOne(iRC,iOpt)
if (iRC /= 0) then
  write(u6,*) 'Error when trying to close the one-electron'
  write(u6,*) 'integral file.'
  call Quit(_RC_INTERNAL_ERROR_)
end if

if (IfVB /= 2) then
  do I=10,99
    inquire(unit=I,opened=IfOpened)
    if (IfOpened .and. (I /= 19)) close(I)
  end do
  close(LUInput)
end if

return

101 format(6X,I3,I4,I5,I5,F15.8,ES12.2,A1,ES10.2,A1,2I4,I2,ES10.2,A1,F6.2,F7.2,4X,A2,3X,A3,I5,A1,I2.2,A1,I2.2)
102 format(3X,I3,I4,I2,I2,F15.8,F15.8,ES12.2,A1,ES10.2,A1,2I4,I2,ES10.2,A1,F6.2,F7.2,4X,A2,3X,A3,I5,A1,I2.2,A1,I2.2)
#ifdef _DMRG_
103 format(6X,I3,I3,I4,I7,ES12.2,I4,I5,F15.8,ES12.2,A1,ES9.2,A1,2I4,I2,ES10.2,A1,F6.2,F7.2,4X,A2,3X,A3,I7,A1,I2.2,A1,I2.2)
#endif
104 format(6X,I3,I4,I5,I5,F15.8,ES12.2,A1,ES10.2,A1,2I4,I2,ES10.2,A1,F6.2,F7.2,4X,A2,3X,A3,I5,A1,I2.2,A1,I2.2)

end subroutine rasscf
