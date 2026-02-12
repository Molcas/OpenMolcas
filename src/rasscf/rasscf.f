************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1989, Per Ake Malmqvist                                *
*               1989, Bjorn O. Roos                                    *
*               1991,1993, Markus P. Fuelscher                         *
*               1991,1993, Jeppe Olsen                                 *
*               1998, Roland Lindh                                     *
************************************************************************
      SUBROUTINE RASSCF(IRETURN)
************************************************************************
*                                                                      *
*           ######     #     #####   #####   #####  #######            *
*           #     #   # #   #     # #     # #     # #                  *
*           #     #  #   #  #       #       #       #                  *
*           ######  #     #  #####   #####  #       #####              *
*           #   #   #######       #       # #       #                  *
*           #    #  #     # #     # #     # #     # #                  *
*           #     # #     #  #####   #####   #####  #                  *
*                                                                      *
*                                                                      *
*                   A program for complete (CAS) and                   *
*                   restricted (RAS)SCF calculations                   *
*                                                                      *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.P. Fuelscher and J. Olsen, P.Aa. Malmqvist and B.O. Roos       *
*     University of Lund, Sweden                                       *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history:                                                         *
*     MOLCAS version 1 by P.Aa. Malmqvist and B.O. Roos, 1989          *
*     MOLCAS version 2 by M.P. Fuelscher and J. Olsen  , 1991          *
*     MOLCAS version 3 by M.P. Fuelscher and J. Olsen  , 1993          *
*                                                                      *
*     Modified to process only unique symmetry blocks, R. Lindh,       *
*     March 1998.                                                      *
*                                                                      *
************************************************************************

#ifdef _DMRG_
!     module dependencies
      use qcmaquis_interface, only: qcmaquis_interface_delete_chkp,
     &  qcmaquis_interface_prepare_hirdm_template,
     &  qcmaquis_interface_deinit, qcmaquis_param,
     &  TEMPLATE_4RDM, TEMPLATE_TRANSITION_3RDM, dmrg_energy
      use qcmaquis_interface_mpssi, only: qcmaquis_mpssi_transform
      use lucia_data, only: RF1, RF2
      use rasscf_global, only: DoNEVPT2Prep, DoDelChk
      use rasscf_global, only: Twordm_qcm, DoMCPDFTDMRG
#endif
      use OneDat, only: sNoNuc, sNoOri
      use Fock_util_global, only: ALGO, DoActive, DoCholesky
      use write_orbital_files, only : OrbFiles, putOrbFile,
     &  write_orb_per_iter
      use filesystem, only: copy_, real_path
      use generic_CI, only: CI_solver_t
      use fciqmc, only: DoNECI, fciqmc_solver_t, tGUGA_in
      use fciqmc_read_RDM, only: dump_fciqmc_mats
      use para_info, only: king
      use fortran_strings, only: str
      use spin_correlation, only: spin_correlation_driver,
     &    orb_range_p, orb_range_q
      use CC_CI_mod, only: Do_CC_CI, CC_CI_solver_t
      use fcidump, only : make_fcidumps, transform, DumpOnly
      use orthonormalization, only : ON_scheme
      use casvb_global, only: ifvb, invec_cvb
#ifdef _FDE_
      use Embedding_global, only: Eemb, embInt, embPot, embPotInBasis,
     &    embPotPath, embWriteEsp
#endif
#ifdef _HDF5_
      use mh5, only: mh5_put_attr, mh5_put_dset
      use csfbas, only: CONF
      use lucia_data, only: CFTP, DStmp, Dtmp
      use raswfn, only: wfn_iter, wfn_energy, wfn_transdens,
     &                  wfn_transsdens
      use rasscf_global, only: lRoots
      use general_data, only: NACTEL,STSYM
#endif
      use OFembed, only: Do_OFemb, FMaux
      use UnixInfo, only: ProgName
      use stdalloc, only: mma_allocate, mma_deallocate
      use rctfld_module, only: lRF
      use Lucia_Interface, only: Lucia_Util
      use wadr, only: DMAT, PMAT, PA, FockOcc, TUVX, FI, FA, DSPN,
     &                D1I, D1A, OccN, CMO, DIAF
      use sxci
      use gugx, only: SGS, CIS, EXS
      use general_data, only: CRVec, CleanMask, CRPROJ
      use gas_data, only: iDOGAS
      use input_ras, only: KeyORBO, KeyORTH, KeyCION, KeyWRMA, KeyTDM,
     &                     KeySSCR, LuInput
      use raswfn, only: cre_raswfn, Wfn_FileID
      use timers, only: TimeCIOpt, TimeInput, TimeOrb, TimeOutput,
     &                  TimeRelax, TimeTotal, TimeTrans, TimeWfn
      use rasscf_global, only: KSDFT, CBLBM, CMAX, DE, DOBLOCKDMRG,
     &                         DoFaro, DoFCIDump,               ECAS,
     &                         ESX, ExFac, FDIAG, HalfQ, iBLBM,
     &                         ICICH, iCIOnly, iExpand, IfCrPr,
     &                         InOCalc, iPr, iPT2, iRLXRoot, iSave_Exp,
     &                         iSymBB, ITER, ITERCI, ITERSX, JBLBM,
     &                         l_casdft,         lSquare, MaxIt, NAC,
     &                         NACPAR, NACPR2, NewFock, nFint, no2m,
     &                         NonEQ, nROOTS, PotNuc, QNSTEP, QNUPDT,
     &                         ROTMax, Start_Vectors, SXShft, Thre,
     &                         ThrSX, THRTE, TMin, Tot_Charge, EMY,
     &                         VIA_DFT, iRoot, Weight, iAdr15, Ener,
     &                         Conv, DoDMRG, iCIRST, KSDFT_Temp
      use SplitCas_Data, only: DoSPlitCas,IterSplit,lRootSplit
      use PrintLevel, only: DEBUG,USUAL,TERSE
      use output_ras, only: LF,IPRLOC,RC_CI,RC_SX
      use general_data, only: NALTER,ITERFILE,NSYM,INVEC,ISPIN,NCONF,
     &                        NCRVEC,JOBIPH,NASH,NBAS,NDEL,NFRO,
     &                        NISH,NRS1,NRS2,NRS3,NTOT,NTOT1,NTOT2
      use spinfo, only: DOBKAP
      use rasscf_global, only: IPCMROOT
      use DWSol, only: DWSolv, DWSol_final, DWSol_init
      use Molcas, only: MxRoot
      use RASDim, only: MxIter

      Implicit None

#include "warnings.h"

      Integer IReturn, RC_RAS
      Logical DSCF
      Logical lTemp, lOPTO
      Character(LEN=80) Line
      Character(LEN=8) Label
      Character(LEN=1) CTHRE, CTHRSX, CTHRTE
      Logical IfOpened
      Real*8 ECAS1, EVAC
#ifdef _DMRG_
      Logical Do_ESPF
      ! function defined in misc_util/pcm_on.f
      Logical, External :: PCM_On
#endif
      Real*8 CASDFT_E, CASDFT_FUNCT,
     &       DiffE, DiffETol, dum1, dum2, dum3, EAv, ThMax, TMXTOT,
     &       time0(2), time1(2), time2(2), time3(2)
      Real*8, External:: Get_ExFac
      Integer i, i_ROOT, iAd, iAd15, iBas,      iComp, iFinal, ihh,
     &        imm, Ind, IndT, iOff, iOpt, iPrLev, iRC, iRot, iShift,
     &        iss, iSyLbl, iSym, iTerm, j,               kau,
     &        kRoot, LuOne, LuvvVec, mRoots, nTav, iFlags, NoScr1
      Integer, External:: IsFreeUnit
#ifdef _HDF5_
      Integer iDX, jDisk, jRoot, kDisk
#endif
#ifdef _FDE_
      Integer iDummyEmb, iEmb, iUnit, nNuc
      Real*8, External:: EmbPotEneMODensities
#endif

* --------- FCIDUMP stuff:
      real*8, allocatable :: orbital_E(:), folded_Fock(:)
* --------- End FCIDUMP stuff:
* --------- CI-solver class
        class(CI_solver_t), allocatable :: CI_solver

! actual_iter starts at 0, so iter 1A == 0, 1B == 1, 2 == 2, 3 == 3 and so on
      integer :: actual_iter


      Character(LEN=15) STLNE2
      External RasScf_Init
      External Scan_Inp
      External Proc_Inp
#ifdef _DMRG_
      integer :: maxtrR
      integer :: maxBD
      real*8  :: maxtrW
#endif
      Integer IndType(56)
      Character(len=80) ::  VecTyp

#ifdef _HDF5_
      Real*8, Allocatable:: VecL(:), VecR(:), Tmp(:)
      Integer, Allocatable:: kcnf(:)
#endif
      Real*8, Allocatable:: Dens(:), PUVX(:), TmpDMat(:), CMON(:),
     &                      OCCX(:), Scr1(:), Scr2(:), SMat(:),
     &                      QMat(:), EDUM(:), Tmp1(:), Fock(:),
     &                      TmpDS(:), TmpD1S(:)
      Integer, External :: isStructure

* Set status line for monitor:
      Call StatusLine('RASSCF: ','Just started.')

* Set the return code(s)
      ITERM  = 0
      IRETURN=_RC_ALL_IS_WELL_

* Set the HDF5 file id (a proper id will never be 0)
      wfn_fileid = 0

* Set some Cholesky stuff
      DoActive=.true.
      lOPTO=.False.
* Initialise doDMRG if compiled without QCMaquis
#ifndef _DMRG_
      DoDMRG = .false.
#endif

* Set variable IfVB to check if this is a VB job.
      IfVB=0
      If (ProgName(1:5).eq.'casvb') IfVB=2
* Default option switches and values, and initial data.
      THMAX=0.0d0
      Call RasScf_Init()
      Call Seward_Init()
* Open the one-olectron integral file:
       LuOne=77
       LuOne=isFreeUnit(LuOne)
       iRC=-1
       iOpt=0
       Call OpnOne(iRC,iOpt,'ONEINT',LuOne)
       If (iRC.ne.0) Then
        Write (6,*) 'Error when trying to open the one-electron'
        Write (6,*) 'integral file.'
        Call Quit(_RC_INTERNAL_ERROR_)
       End If
       Call StatusLine('RASSCF: ','Read-in ONEINT')
      If (IfVB.eq.2) go to 10

*
* Make a copy, upper-cased, left-adjusted, of the input between and including
* the '&RASSCF' and the 'End of input' markers, skipping all lines beginning
* with '*' or '!' or ' '  when left-adjusted, and replacing any rightmost
* substring beginning with '!' with blanks.
* That copy will be in file 'CleanInput', and its unit number is returned
* as LUInput in common (module file input_ras.F90) by the following call:
      Call cpinp(LUInput,iRc)
* If something wrong with input file:
      If (iRc.ne._RC_ALL_IS_WELL_) Then
       Call WarningMessage(2,'Input file is unusable.')
       Write(6,*)' RASSCF Error: Could not make a clean copy of'
       Write(6,*)' the input file. This is an unexpected bug.'
       IRETURN=_RC_INTERNAL_ERROR_
       GOTO 9990
      End If

* Scan the input file for keywords:
      Call Scan_Inp(iRc)
* If something wrong with input file:
      If (iRc.ne._RC_ALL_IS_WELL_) Then
       If (IPRLOC(1).GE.TERSE) Then
        Call WarningMessage(2,'Scanning input file failed.')
* Calling again, now with iRc indicating an error, will echo the keywords:
        Call Scan_Inp(iRc)
       End If
       IRETURN=_RC_INPUT_ERROR_
       GOTO 9990
      End If

* Local print level in this routine:
      IPRLEV=IPRLOC(1)
*
10    CONTINUE
* Open files
      Call OpnFls_RASSCF(DSCF,DoCholesky)

* Some preliminary input data:
      Call Rd1Int()
      If ( .not.DSCF ) Call Rd2Int_RASSCF()

* Printed program header:

* Process the input:
      Call StatusLine('RASSCF: ','Processing input')
      Call Proc_Inp(DSCF,lOPTO,iRc)
* If something goes wrong in proc_inp:
      If (iRc.ne._RC_ALL_IS_WELL_) Then
       If (IPRLEV.ge.TERSE) Then
        Call WarningMessage(2,'Input processing failed.')
        Write(6,*)' RASSCF Error: Proc_Inp failed unexpectedly.'
        Write(6,*)' Check the output file for any previous messages'
        Write(6,*)' that can help explain the failure.'
        Write(6,*)' Here is a printing of the input file that'
        Write(6,*)' was processed:'
        Rewind(LUInput)
  15    Continue
        Read(LuInput,'(A80)',End=16,Err=16) Line
        Write(6,*) Line
        Go To 15
  16    Continue
       End If
       IRETURN=iRc
       GOTO 9990
      End If
      if (lRF) call DWSol_init(IPCMROOT,nRoots,NonEq)
*     call DWSCF_init(1,nRoots)


* Local print level may have changed:
      IPRLEV=IPRLOC(1)


      Call InpPri(lOpto)

* Note that CI_solver subclasses provide a cleanup procedure
* (C++ people might call it destructor). Hence the deallocation and
* cleanup is automatically performed, when it goes out of scope.
      if (DoNECI) then
        allocate(CI_solver, source=fciqmc_solver_t(tGUGA_in))
      else if (Do_CC_CI) then
        allocate(CI_solver, source=CC_CI_solver_t())
      end if

*
* If this is not CASDFT make sure the DFT flag is unset
*
      If (KSDFT(1:3).eq.'SCF') Then
        Call Get_iScalar('System BitSwitch',iFlags)
        iFlags=iAnd(iFlags,Not(2**6))
        Call Put_iScalar('System BitSwitch',iFlags)
      End If

* If the ORBONLY option was chosen, then Proc_Inp just generated
*  orbitals from the JOBIPH file. Nothing more to do:
      IF((KeyORBO).or.(MAXIT.eq.0)) GOTO 9989
#ifdef _DMRG_
      ! delete old checkpoints, unless requested otherwise
      ! this flag is set in proc_inp
      if (.not.DoDelChk) then
        do kroot=1,nroots
          call qcmaquis_interface_delete_chkp(iroot(kroot))
        end do
      end if
#endif
*
* Allocate various matrices
*
      Call mma_allocate(FI,NTOT1,Label='FI')
      Call mma_allocate(FA,NTOT1,Label='FA')
      Call mma_allocate(D1I,NTOT2,Label='D1I')
      Call mma_allocate(D1A,NTOT2,Label='D1A')
      Call mma_allocate(OCCN,NTOT,Label='OccN')
      Call mma_allocate(CMO,NTOT2,Label='CMO')
      Call mma_allocate(DIAF,NTOT,Label='DIAF')
#ifdef _DMRG_
* Allocate RDMs for the reaction field reference root in QCMaquis calculations
      if (doDMRG.and.PCM_On()) then
        Call mma_allocate(RF1,NACPAR,Label='RF1')
        if (twordm_qcm) then
          Call mma_allocate(RF2,NACPR2,Label='RF2')
        end if
      end if
#endif
      FI(:)=0.0D0
      FA(:)=0.0D0
      DIAF(:)=0.0D0
      ECAS1=0.0D0
      EVAC=0.0D0
*
      If (iCIRST.eq.1.and.DumpOnly) then
        write(6,*) 'ICIRST and DumpOnly flags are not compatible!'
        write(6,*) 'Choose only one.'
        Call Abend
      end if

      If(DumpOnly) then
        write(6,*) 'Dumping integrals.'
        write(6,*) 'Nothing else will be done.'
      end if

        Call mma_allocate(TUVX,NACPR2,Label='TUVX')
        TUVX(:)=0.0D0
        Call mma_allocate(DSPN,NACPAR,Label='DSPN')
        DSPN(:)=0.0D0
        Call mma_allocate(DMAT,NACPAR,Label='DMat')
        DMAT(:)=0.0D0
        Call mma_allocate(PMAT,NACPR2,Label='PMat')
        PMAT(:)=0.0D0
        Call mma_allocate(PA,NACPR2,Label='PA')
        PA(:)=0.0D0
#ifdef _FDE_
      ! Embedding
      iDummyEmb=0
      Call Get_iScalar('embpot', iDummyEmb)
      if (iDummyEmb.eq.1) embPot=.true.
      if (embPot) then
       Call EmbPotRdRun
      end if
      if (embpot) then
       ! I have no idea why i need memory for x+4 entries
       ! and not just x...
       call mma_allocate(embInt,NTOT1+4,label='Emb')
       if (embPotInBasis) then
        ! If the potential is given in basis set representation it
        ! has not been calculated with a OneEl call and is just read
        ! from file here.
        iunit = isFreeUnit(1)
        call molcas_open(iunit, embPotPath)
        do iEmb=1, NTOT1
         read(iunit,*) embInt(iEmb)
        end do
       else
        ! Read in the embedding potential one-electron integrals
        Label='embpot  '
        iRC=-1
        iOpt=0
        iComp=1
        Call RdOne(iRC,iOpt,Label,iComp,embInt,iSyLbl)
        If (iRC.ne.0) then
         Call WarningMessage(2,
     &                'Drv1El: Error reading ONEINT;'
     &              //'Label='//Label)
         Call Quit(_RC_IO_ERROR_READ_)
        End If
       end if
      end if
#endif

*
* Get start orbitals

* Initialize OCCN array, to prevent false alarms later from
* automated detection of using uninitialized variables:
      OccN(:)=0.0D0

* PAM03: Note that removal of linear dependence may change the nr
* of secondary/deleted orbitals, affecting some of the global
* variables: NSSH(),NDEL(),NORB(),NTOT3, etc etc
      Call ReadVc(CMO,OCCN,DMAT,DSPN,PMAT,PA,ON_scheme)
      if (KeyORTH) then
! TODO(Oskar): Add fourth argument OCC
!   If the Occupation number is written properly as well.
        call putOrbFile(CMO=CMO(:),
     &                  orbital_E=DIAF(:),
     &                  iDoGAS=iDoGAS)
      end if
* Only now are such variables finally known.

*
* Allocate core space for dynamic storage of data
*
      CALL ALLOC()
*
* Create job interphase on unit JOBIPH (FT15)
*

      if(ifvb.ne.2) then
        CALL CREIPH()
        call cre_raswfn()
      end if
      if(ifvb.eq.1)call casinfo2_cvb()

      Call Timing(dum1,dum2,time1(1),dum3)
      TimeInput = time1(1)

CGG03 Aug 03
      If(NAlter.gt.0) Call Alter_MO(CMO)

c At this point all is ready to potentially dump MO integrals... just do it if required.
      If(DumpOnly) goto 20
      if(ifvb.eq.2)goto 20

      if(dofcidump)then
        Write(LF,*)
        Write(LF,'(26X,A)')
     &  'Dumping integrals on file FCIDUMP - nothing else to be done'
        goto 20
      end if

************************************************************************
*
* Wave function section
*
************************************************************************

      Call StatusLine('RASSCF: ','Compute wave function.')
      If ( IPRLEV.GE.2 .AND..NOT.lOPTO) then
       Write(LF,*)
       Write(LF,'(6X,A)') repeat('*',120)
       Write(LF,'(6X,A,118X,A)') '*','*'
       Write(LF,'(6X,A,44X,A,45X,A)')
     &      '*','Wave function control section','*'
       Write(LF,'(6X,A,118X,A)') '*','*'
       Write(LF,'(6X,A)') repeat('*',120)
       Write(LF,*)
      End If

      If ( IPRLEV.GE.2 .AND..NOT.lOPTO) then
       IF (ICIONLY.eq.0) THEN
        Write(LF,*)
        if(doDMRG)then
          Write(LF,'(41X,A)')
     &         'DMRGSCF iterations: Energy and convergence statistics'
          Write(LF,'(41X,A)')
     &         '-----------------------------------------------------'
        else
          Write(LF,'(41X,A)')
     &         'RASSCF iterations: Energy and convergence statistics'
          Write(LF,'(41X,A)')
     &         '----------------------------------------------------'
        end if
        Write(LF,*)
       ELSE
        Write(LF,*)
        if(doDMRG)then
          Write(LF,'(41X,A)')
     &        'DMRGCI only, no orbital optimization will be done.'
          Write(LF,'(41X,A)')
     &        '--------------------------------------------------'
        else
          Write(LF,'(41X,A)')
     &        'CASCI only, no orbital optimization will be done.'
          Write(LF,'(41X,A)')
     &        '-------------------------------------------------'
        end if
        Write(LF,*)
       END IF
#ifdef _DMRG_
       if(doDMRG)then
          Write(LF,'(45x,a//,36x,a/,36x,a/,36x,a//,45x,a//,'//
     &               '36x,a/,36x,a/,36x,a//,36x,a/,36x,a,a/,36x,a//)')
     &        'Please cite for the QCMaquis-Molcas driver:',
     &        'Freitag L.; Keller S.; Knecht S.; Lindh R.; Ma Y.; ',
     &        'Stein C. J. and Reiher M., in preparation. (2018).',
     &        '------------------------------------------------'//
     &        '---------------',
     &        'Please cite for the QCMaquis DMRG software:',
     &        'S. Keller, M. Dolfi, M. Troyer, M. Reiher,',
     &        'J. Chem. Phys. 143, 244118 (2015)',
     &        '------------------------------------------------'//
     &        '---------------'
       end if
#endif
       IF(INOCALC.EQ.1) THEN
         WRITE(LF,*)
         WRITE(LF,'(26X,A)')
     &      ' No calculation will be performed. Stopping in LUCIA'
       END IF
       IF(ISAVE_EXP.EQ.1) THEN
         WRITE(LF,*)
         WRITE(LF,'(26X,A)')
     &      ' Information on the CI-vector will be written to ???.'
       END IF
       IF(IEXPAND.EQ.1) THEN
         WRITE(LF,*)
         WRITE(LF,'(26X,A)')
     &      ' A shorter vector will be expanded in a longer'
       END IF
       If ( IPRLEV.LE.3) then
        If(DoSplitCAS) then
         Write(LF,'(6X,A)')
     &         'Iter CI   SX   CI       SplitCAS       Energy    '//
     &         'max ROT     max BLB   max BLB  Level Ln srch  Step '//
     &         '  QN   Walltime'
         Write(LF,'(6X,A)')
     &         '    iter iter root      energy       change    '//
     &         ' param      element    value   shift minimum  type '//
     &         'update hh:mm:ss'
        else if (DoBKAP) then
         Write(LF,'(6X,A)')
     &         'Iter CI   SX   CI   RASSCF      CI    Energy    '//
     &         'max ROT     max BLB   max BLB  Level Ln srch  Step '//
     &         '  QN   Walltime'
         Write(LF,'(6X,A)')
     &         '    iter iter root  energy    energy  change    '//
     &         ' param      element    value   shift minimum  type '//
     &         'update hh:mm:ss'
        else if (DoDMRG .and. ICIONLY == 0)then
         Write(LF,'(6X,A)')
     &         'Iter num   Bond  DMRG max tr DMRG  SX'//
     &         '      DMRGSCF       Energy    '//
     &         'max ROT   max BLB     max BLB  Level Ln srch  Step '//
     &         '  QN     CPU Time'
         Write(LF,'(6X,A)')
     &         '   sweeps/ dim  /root weight/root iter'//
     &         '     energy        change    '//
     &         ' param    element      value   shift minimum  type '//
     &         'update   hh:mm:ss'
        else if (DoDMRG .and. ICIONLY /= 0)then

        else if( l_casdft ) then

        else
         Write(LF,'(6X,A)')
     &         'Iter CI   SX   CI       RASSCF       Energy    '//
     &         'max ROT     max BLB   max BLB  Level Ln srch  Step '//
     &         '  QN   Walltime'
         Write(LF,'(6X,A)')
     &         '    iter iter root      energy       change    '//
     &         ' param      element    value   shift minimum  type '//
     &         'update hh:mm:ss'
        end if
       End If
      End If
20    continue
*                                                                      *
************************************************************************
*                                                                      *
*     Start iterations
*                                                                      *
************************************************************************
*                                                                      *
      Rc_CI  = 0
      Rc_SX  = 0
      ECAS   = 0.0d0
      ROTMAX = 0.0d0
      ITER   = 0
      actual_iter = 0
      IFINAL = 0
      TMXTOT = 0.0D0
      Call mma_allocate(FockOcc,nTot1,Label='FockOcc')
*                                                                      *
************************************************************************
*                                                                      *
*     Entry point for second and successive iterations
*                                                                      *
************************************************************************
*                                                                      *
 1000 CONTINUE

      if( l_casdft ) then
        KSDFT_TEMP=KSDFT
        KSDFT='SCF'
        ExFac=1.0D0
      else
         KSDFT_TEMP=KSDFT
         ExFac=Get_ExFac(KSDFT)
      end if

      ITER=ITER+1
      Write(STLNE2,'(A12,I3)')'Iteration ',ITER
      Call StatusLine('RASSCF: ',STLNE2)
      Call Timing(dum1,dum2,time0(1),dum3)
#ifdef _DMRG_
      ! Leon 27/11/2017: Skip the first CI iteration if we're using
      ! DMRGCI and CIOnly.It's enabled only for DMRGCI with QCMaquis
      ! now, (to exclude potential side effects)
      ! but consider extending it to other cases!
      call DecideOnESPF(Do_ESPF)
      !write(LF,*) ' |rasscf> DecideOnESPF == ',Do_ESPF
      If (( ITER.EQ.1 ).and.((.not.(DoDMRG.and.(ICIONLY.NE.0))).or.lRf
     &    .or.domcpdftDMRG.or.Do_ESPF))THEN
#else
      If ( ITER.EQ.1 ) THEN
#endif
************************************************************************
*     ^   First iteration
************************************************************************
*

* Print header to file containing informations on CI iterations.
*
        Write(IterFile,'(20A4)') ('****',i=1,20)
        Write(IterFile,'(15X,A)') 'RASSCF iteration: 1A'
*
        Start_Vectors=.True.
        lTemp = lRf
*
* Transform two-electron integrals and compute at the same time
* the Fock matrices FI and FA
*
        Call Timing(dum1,dum2,time2(1),dum3)

        If (.not.DoCholesky .or. ALGO.eq.1) Then
           Call mma_allocate(PUVX,NFINT,Label='PUVX')
           PUVX(:)=0.0D0
        EndIf

        Call Get_D1I_RASSCF(CMO,D1I)
        If ( IPRLEV.ge.DEBUG ) then
          Write(LF,*)
          Write(LF,*) ' D1I in AO basis in RASSCF'
          Write(LF,*) ' ---------------------'
          Write(LF,*)
          iOff=1
          Do iSym = 1,nSym
            iBas = nBas(iSym)
             call wrtmat(D1I(ioff),iBas,iBas, iBas, iBas)
            iOff = iOff + iBas*iBas
          End Do
        End If

* Compute D1A from CMO coefficients and, if CIREstart, old DMAT.
        If (iCIRST.eq.1) Then

           Call mma_allocate(TmpDMAT,NACPAR,Label='TmpDMAT')
           call dcopy_(NACPAR,DMAT,1,TmpDMAT,1)
           If (NASH(1).ne.NAC) Call DBLOCK(TmpDMAT)
           Call Get_D1A_RASSCF(CMO,TmpDMAT,D1A)
           Call mma_deallocate(TmpDMAT)

           DoActive = .true.

        Else

           lRf = .false.
           IF( .not.l_casdft )  then
             KSDFT='SCF'
             ExFac=1.0D0
           end IF
           D1A(:)=0.0D0

           DoActive = .false.

        End If

        IPR=0
        IF(IPRLOC(2).EQ.4) IPR=5
        IF(IPRLOC(2).EQ.5) IPR=10

        If ( IPRLEV.ge.DEBUG ) then
          Write(LF,*)
          Write(LF,*) ' PUVX in rasscf bf first TRACTL2'
          Write(LF,*) ' ---------------------'
          Write(LF,*)
          call wrtmat(PUVX,1,nFint, 1, nFint)

          Write(LF,*)
          Write(LF,*) ' ---------------------'
          Write(LF,*)
          Write(LF,*) ' D1A in AO basis in RASSCF bf TRACTL2 1'
          Write(LF,*) ' ---------------------'
          Write(LF,*)
          iOff=1
          Do iSym = 1,nSym
            iBas = nBas(iSym)
            call wrtmat(D1A(ioff),iBas,iBas, iBas, iBas)
            iOff = iOff + iBas*iBas
          End Do
        end if

*
* Transform two-electron integrals and compute the Fock matrices FI and FA
* FI and FA are output from TRACTL2...
        CALL TRACTL2(CMO,PUVX,TUVX,D1I,
     &               FI,D1A,FA,IPR,lSquare,ExFac)

c         Write(6,*) ' TUVX after TRACTL2'
c         write(6,*) (UVX(ind),ind=1,NACPR2)
        IF (ITER.eq.1 .and. IfCRPR) Then
* Core shift applied to projection of WF with doubly occupied core
          Call MkCRVEC(CMO,CRVEC)
        END IF

        If ( IPRLEV.ge.DEBUG ) then
          Write(LF,*)
          Write(LF,*) ' D1A in AO basis in RASSCF af TRACTL2 1'
          Write(LF,*) ' ---------------------'
          Write(LF,*)
          iOff=1
          Do iSym = 1,nSym
            iBas = nBas(iSym)
            call wrtmat(D1A(ioff),iBas,iBas, iBas, iBas)
            iOff = iOff + iBas*iBas
          End Do

          Write(LF,*)
          Write(LF,*) ' ---------------------'
          Write(LF,*)
          Write(LF,*) ' PUVX in rasscf af first TRACTL2'
          Write(LF,*) ' ---------------------'
          Write(LF,*)
          call wrtmat(PUVX,1,nFint, 1, nFint)

          Write(LF,*)
          Write(LF,*) ' ---------------------'

          Write(6,*) ' TUVX after TRACTL2'
          write(6,*) (TUVX(ind),ind=1,NACPR2)
          Write(LF,*)
          Write(LF,*) ' ---------------------'
        end if

        If (.not.DoCholesky .or. ALGO.eq.1) Then
          Call mma_deallocate(PUVX)
        EndIf

        Call Timing(dum1,dum2,time2(2),dum3)
        TimeTrans = TimeTrans + time2(2) - time2(1)

        If ( IPRLEV.ge.DEBUG ) then
          Write(LF,*)
          Write(LF,*) ' CMO in RASSCF bf first call to CICTL'
          Write(LF,*) ' ---------------------'
          Write(LF,*)
          ioff=1
          Do iSym = 1,nSym
            iBas = nBas(iSym)
            if(iBas.ne.0) then
              write(6,*) 'Sym =', iSym
              do i= 1,iBas
                write(6,*) (CMO(ioff+iBas*(i-1)+j),j=0,iBas-1)
              end do
              iOff = iOff + (iBas*iBas)
            end if
          End Do
        End If

*
* Compute initial CI vectors and density matrices
*
        Call Timing(dum1,dum2,time3(1),dum3)

        if (DumpOnly) then
          call mma_allocate(orbital_E, nTot)
          call mma_allocate(folded_Fock, nAcPar)
          call transform(iter,
     &                   CMO=CMO(:),
     &                   DIAF=DIAF(:),
     &                   D1I_AO=D1I(:),
     &                   D1A_AO=D1A(:),
     &                   D1S_MO=DSPN(:),
     &                   F_IN=FI(:),
     &                   orbital_E=orbital_E,
     &                   folded_Fock=folded_Fock)
          call make_fcidumps('FCIDUMP', 'H5FCIDUMP',
     &      orbital_E, folded_Fock,
     &      TUVX=tuvx(:), core_energy=EMY)
          call mma_deallocate(orbital_E)
          call mma_deallocate(folded_Fock)
          write(6,*) "FCIDMP file generated. Here for serving you!"
          goto 2010
        end if

        if (allocated(CI_solver)) then
          call CI_solver%run(actual_iter=actual_iter,
     &                    ifinal=ifinal,
     &                    iroot=iroot,
     &                    weight=weight,
     &                    CMO=CMO(:),
     &                    DIAF=DIAF(:),
     &                    D1I_AO=D1I(:),
     &                    D1A_AO=D1A(:),
     &                    TUVX=tuvx(:),
     &                    F_IN=FI(:),
     &                    D1S_MO=DSPN(:),
     &                    DMAT=DMAT(:),
     &                    PSMAT=pmat(:),
     &                    PAMAT=pa(:))

#if defined (_ENABLE_BLOCK_DMRG_) || defined (_ENABLE_CHEMPS2_DMRG_) || defined (_ENABLE_DICE_SHCI_)
        else If(DoBlockDMRG) then
          CALL DMRGCTL(CMO,DMAT,DSPN,PMAT,PA,FI,D1I,D1A,TUVX,IFINAL,0)
#endif
        else
          CALL CICTL(CMO,DMAT,DSPN,PMAT,PA,FI,FA,D1I,D1A,TUVX,IFINAL)

          if(dofcidump)then
           write(LF,*) " FCIDUMP file generated. This is the end..."
           goto 9990
          end if
#ifdef _FDE_
          !Thomas Dresselhaus
          if (embpot) then
            !Eemb=DDot_(NACPAR,embInt,1,DMAT,1)
!           Eemb=embPotEne(D1I, D1A, embInt,
!    &                     CMO, nBasFunc, nFrozenOrbs, .true.)
            Eemb=embPotEneMODensities(D1I, D1A,
     &            embInt, nBas, nTot2, nSym)
            Write(LF,*) "Energy from embedding potential with the"
            Write(LF,*) "initial CI vectors: ", Eemb
          end if
        !!!!!!!!!!!!!!!!!!!
#endif
* PAM 2015: Additional output line.
          If ( IPRLEV.ge.USUAL .and. .not. doDMRG) Then
             write(6,'(a,i4)')' Nr of preliminary CI iterations:',ITERCI
          End If
        end if

c.. dongxia testing jobiph
c.. upt to here, jobiph are all zeros at iadr15(2)
* If CASVB job, go directly to return.
        invec_cvb=invec
        If (DSCF) NewFock=1
        If (IfVB.eq.2) GoTo 9990

        EAV=0.0d0
        If(DoSplitCAS) then
          EAV = ENER(lRootSplit,ITER)
        Else
          DO KROOT=1,NROOTS
            EAV=EAV+ENER(IROOT(KROOT),ITER)*WEIGHT(KROOT)
          END DO
        End If

        Call Get_D1A_RASSCF(CMO,DMAT,D1A)

        If ( IPRLEV.ge.DEBUG ) then
          Write(LF,*)
          Write(LF,*) ' D1A in AO basis in RASSCF af Get_D1A_RASSCF'
          Write(LF,*) ' ---------------------'
          Write(LF,*)
          iOff=1
          Do iSym = 1,nSym
            iBas = nBas(iSym)
            call wrtmat(D1A(ioff),iBas,iBas, iBas, iBas)
            iOff = iOff + iBas*iBas
          End Do
        end if

        Call Timing(dum1,dum2,time3(2),dum3)
        TimeCIOpt = TimeCIOpt + time3(2) - time3(1)
        lRf = lTemp

        IF( .not.l_casdft ) then
          KSDFT=KSDFT_TEMP
          ExFac=Get_ExFac(KSDFT)
        end IF

*     v GLM for MC-PDFT
      If(KSDFT.ne.'SCF'.and.KSDFT.ne.'PAM'.or.l_casdft) Then
        If ( IPRLEV.ge.DEBUG ) then
          Write(LF,*)
          Write(LF,*) ' CMO in RASSCF bf call NATORB_RASSCF'
          Write(LF,*) ' ---------------------'
          Write(LF,*)
          ioff=1
          Do iSym = 1,nSym
            iBas = nBas(iSym)
            if(iBas.ne.0) then
              write(6,*) 'Sym =', iSym
              do i= 1,iBas
                write(6,*) (CMO(ioff+iBas*(i-1)+j),j=0,iBas-1)
              end do
              iOff = iOff + (iBas*iBas)
            end if
          End Do
        End If

        Call mma_allocate(CMON,NTOT2,Label='CMON')
        Call mma_allocate(OCCX,NTOT,Label='OCCX')
        noscr1=max(nacpar,no2m)
        Call mma_allocate(scr1,noscr1,Label='Scr1')
        Call mma_allocate(scr2,NO2M,Label='Scr2')
        Call mma_allocate(SMAT,NTOT1,Label='SMAT')
        CALL NATORB_RASSCF(CMO,scr1,scr2,SMAT,CMON,OCCX)
        Call dCopy_(NTOT2,CMON,1,CMO,1)
        Call Put_dArray('Last orbitals',CMO,ntot2)
        Call mma_deallocate(scr1)
        Call mma_deallocate(scr2)
        Call mma_deallocate(SMAT)
        Call mma_deallocate(OCCX)
        Call mma_deallocate(CMON)
        If ( IPRLEV.ge.DEBUG ) then
         Write(LF,*)
         Write(LF,*) ' CMO in RASSCF af call NATORB_RASSCF & bf 2 CICTL'
         Write(LF,*) ' ---------------------'
         Write(LF,*)
         ioff=1
         Do iSym = 1,nSym
           iBas = nBas(iSym)
           if(iBas.ne.0) then
             write(6,*) 'Sym =', iSym
             do i= 1,iBas
               write(6,*) (CMO(ioff+iBas*(i-1)+j),j=0,iBas-1)
             end do
             iOff = iOff + (iBas*iBas)
           end if
         End Do
        End If
      End If
*     ^ GLM End If for MC-PDFT

        EAV=0.0d0
        DO KROOT=1,NROOTS
         EAV=EAV+ENER(IROOT(KROOT),ITER)*WEIGHT(KROOT)
        END DO
        If ( IPRLEV.ge.DEBUG ) then
          IF( l_casdft ) then
            write(6,*) 'EAV value in RASSCF after first call to CICTL:'
            write(6,*) EAV
          END if
        End if
      END IF
************************************************************************
*     ^ End First iteration
************************************************************************
*
* Print header to file containing informations on CI iterations.
*
      actual_iter = actual_iter + 1
      Write(IterFile,*)
      Write(IterFile,'(20A4)') ('****',i=1,20)
      IF (Iter .Eq. 1) Then
         Write(IterFile,'(15X,A)') 'RASSCF iteration: 1B'
      Else
         Write(IterFile,'(15X,A,I3)') 'RASSCF iteration: ',Iter
      End If
*
      If ( IPRLEV.ge.DEBUG .and. l_casdft) then
        write(6,*) repeat('*',70)
        write(6,*) 'we are done withe first standard CAS-CI iteration  '
        write(6,*) 'CI coeffs are known and mantained fix in next stage'
        write(6,*) 'We are now going to remove exchange from FI and FA '
        write(6,*) 'in TRACTL2 --> TRA_CTL2 --> TRADRV --> FTWO        '
        write(6,*) 'the ExFac is going to be set to 0.0d0 as DT asked! '
        write(6,*) 'FI and FA are going to change... '
        write(6,*) 'Check with previous printout to see differences.\  '
        write(6,*) repeat('*',70)
      End if


      IF( l_casdft ) then
        KSDFT=KSDFT_TEMP
        ExFac=0.0d0
      end IF


      IF(ICIONLY.NE.0) IFINAL=1
*
* Transform two-electron integrals and compute at the same time
* the Fock matrices FI and FA
*
      Call Timing(dum1,dum2,time2(1),dum3)
      If (.not.DoCholesky .or. ALGO.eq.1) Then
         Call mma_allocate(PUVX,NFINT,Label='PUVX')
         PUVX(:)=0.0D0
      EndIf
      Call Get_D1I_RASSCF(CMO,D1I)

      DoActive = .true.

      If (DoCholesky.and.ALGO.eq.2) Then
         NTav=0
         do iSym=1,nSym
            NTav = NTav + nBas(iSym)*nAsh(iSym)
         end do
         Call mma_allocate(Qmat,NTav,Label='QMat')
         QMat(:)=0.0D0
      EndIf

       IPR=0
       IF(IPRLOC(2).EQ.4) IPR=5
       IF(IPRLOC(2).EQ.5) IPR=10
       If ( IPRLEV.ge.DEBUG ) then
         Write(LF,*)
         Write(LF,*) ' D1A in AO basis in RASSCF bf TRACTL2 2'
         Write(LF,*) ' ---------------------'
         Write(LF,*)
         iOff=1
         Do iSym = 1,nSym
          iBas = nBas(iSym)
          call wrtmat(D1A(ioff),iBas,iBas, iBas, iBas)
          iOff = iOff + iBas*iBas
         End Do
       end if
       CALL TRACTL2(CMO,PUVX,TUVX,D1I,
     &              FI,D1A,FA,IPR,lSquare,ExFac)

      If ( IPRLEV.ge.DEBUG ) then
         Write(LF,*)
         Write(LF,*) ' D1A in AO basis in RASSCF af TRACTL2 2'
         Write(LF,*) ' ---------------------'
         Write(LF,*)
         iOff=1
         Do iSym = 1,nSym
          iBas = nBas(iSym)
          call wrtmat(D1A(ioff),iBas,iBas, iBas, iBas)
          iOff = iOff + iBas*iBas
         End Do
      end if

      If (.not.DoCholesky .or. ALGO.eq.1) Then
         Call mma_deallocate(PUVX)
      EndIf

      Call Timing(dum1,dum2,time2(2),dum3)
      TimeTrans = TimeTrans + time2(2) - time2(1)

*
* Compute the CI vectors and density matrices
*
      Call Timing(dum1,dum2,time3(1),dum3)
      IF (.not. l_casdft) THEN !the following is skipped in CASDFT-GLM

        If(KSDFT.ne.'SCF'.and.KSDFT.ne.'PAM') Then
          Call Put_dArray('Last orbitals',CMO,ntot2)
        End If

        if (allocated(CI_solver)) then
         ! The following is adapted from sxctl.f
         ! In addition to writing the last RasOrb to disk, the current
         ! orbitals have to be dumped *before* the CI step. The PERI
         ! keyword writes only the orbitals from the last iteration.
         iShift=0
         DO ISYM=1,NSYM
           IndT=0
           IndType(1+iShift)= NFRO(ISYM)
           IndT=IndT+NFRO(ISYM)
           IndType(2+iShift)= NISH(ISYM)
           IndT=IndT+NISH(ISYM)
           IndType(3+iShift)= NRS1(ISYM)
           IndT=IndT+NRS1(ISYM)
           IndType(4+iShift)= NRS2(ISYM)
           IndT=IndT+NRS2(ISYM)
           IndType(5+iShift)= NRS3(ISYM)
           IndT=IndT+NRS3(ISYM)
           IndType(7+iShift)= NDEL(ISYM)
           IndT=IndT+NDEL(ISYM)
           IndType(6+iShift)= NBAS(ISYM)-IndT
           iShift=iShift+7
         EndDo
         call mma_allocate(EDUM,NTOT,Label='EDum')
         EDum(:)=0.0D0
         Write(VecTyp,'(A)')
         VecTyp='* RASSCF average (pseudo-natural) orbitals (Not final)'
         LuvvVec=50
         LuvvVec=isfreeunit(LuvvVec)
         call WrVec('IterOrb',LuvvVec,'COE',NSYM,NBAS,
     &               NBAS, CMO(:), OCCN,
     &               EDUM, INDTYPE,VECTYP)
         call WrVec('IterOrb',LuvvVec,'AI',NSYM,NBAS,
     &               NBAS, CMO(:), OCCN,
     &               EDUM, INDTYPE,VECTYP)
         call mma_deallocate(EDUM)
         write(6,*) "MO coeffs for next iteration written to IterOrb."

          call CI_solver%run(actual_iter=actual_iter,
     &                    ifinal=ifinal,
     &                    iroot=iroot,
     &                    weight=weight,
     &                    CMO=CMO(:),
     &                    DIAF=DIAF(:),
     &                    D1I_AO=D1I(:),
     &                    D1A_AO=D1A(:),
     &                    TUVX=tuvx(:),
     &                    F_IN=FI(:),
     &                    D1S_MO=DSPN(:),
     &                    DMAT=DMAT(:),
     &                    PSMAT=pmat(:),
     &                    PAMAT=pa(:))
#if defined (_ENABLE_BLOCK_DMRG_) || defined (_ENABLE_CHEMPS2_DMRG_) || defined (_ENABLE_DICE_SHCI_)
        else If(DoBlockDMRG) Then
            CALL DMRGCTL(CMO,DMAT,DSPN,PMAT,PA,FI,D1I,D1A,TUVX,IFINAL,1)
#endif
        else
          CALL CICTL(CMO,DMAT,DSPN,PMAT,PA,FI,FA,D1I,D1A,TUVX,IFINAL)
        end if

c      call triprt('twxy',' ',TUVX,nAc*(nAc+1)/2)
c      call triprt('P-mat 2',' ',PMAT,nAc*(nAc+1)/2)

        EAV=0.0d0
        If (DoSplitCAS) Then
          EAV = ENER(lRootSplit,ITER)
        Else
         DO KROOT=1,NROOTS
           EAV=EAV+ENER(IROOT(KROOT),ITER)*WEIGHT(KROOT)
         END DO
        End If

        IF (IPRLEV.ge.DEBUG) THEN
          write(6,*) 'EAV value in RASSCF after second call to CICTL:'
          write(6,*) EAV

          write(6,*) 'Printing matrices in RASSCF'
          Write(LF,*)
          Write(LF,*) ' CMO in RASSCF'
          Write(LF,*) ' ---------------------'
          Write(LF,*)
          ioff=1
          Do iSym = 1,nSym
            iBas = nBas(iSym)
            if(iBas.ne.0) then
              write(6,*) 'Sym =', iSym
              do i= 1,iBas
                write(6,*) (CMO(ioff+iBas*(i-1)+j),j=0,iBas-1)
              end do
              iOff = iOff + (iBas*iBas)
            end if
          End Do

          Write(LF,*)
          Write(LF,*) ' D1I in AO basis in RASSCF'
          Write(LF,*) ' ---------------------'
          Write(LF,*)
          iOff=1
          Do iSym = 1,nSym
            iBas = nBas(iSym)
            call wrtmat(D1I(ioff),iBas,iBas, iBas, iBas)
            iOff = iOff + iBas*iBas
          End Do
          write(6,*)
          write(6,*) 'Total Charge :', Tot_Charge

          Call mma_allocate(Tmp1,nTot1,Label='Tmp1')
          iComp  =  1
          iSyLbl =  1
          iRc    = -1
          iOpt   =  ibset(ibset(0,sNoOri),sNoNuc)
          Label  = 'OneHam'
          Call RdOne(iRc,iOpt,Label,iComp,Tmp1,iSyLbl)
          If ( iRc.ne.0 ) then
           Write(LF,*) 'SGFCIN: iRc from Call RdOne not 0'
#ifdef _FDE_
           Write(LF,*) 'Label = ',Label
#endif
           Write(LF,*) 'iRc = ',iRc
           Call Abend
          End if

          Write(LF,*)
          Write(LF,*) ' OneHam in AO basis in RASSCF'
          Write(LF,*) ' ---------------------'
          Write(LF,*)
          iOff=1
          Do iSym = 1,nSym
            iBas = nBas(iSym)
            Call TriPrt(' ','(5G17.11)',Tmp1(iOff),iBas)
            iOff = iOff + (iBas*iBas+iBas)/2
          End Do

          Call mma_deallocate(Tmp1)
          Call Get_dScalar('PotNuc',potNuc)

          write(6,*)
          write(6,*) 'PotNuc :', PotNuc

          Write(LF,*)
          Write(LF,*) ' D1A in AO basis in RASSCF'
          Write(LF,*) ' ---------------------'
          Write(LF,*)
          iOff=1
          Do iSym = 1,nSym
            iBas = nBas(iSym)
            call wrtmat(D1A(ioff),iBas,iBas, iBas, iBas)
            iOff = iOff + iBas*iBas
          End Do
        End if
      else
        CALL mma_allocate(FOCK,NACPAR,Label='Fock')
* To fix the DS bug... I forgot to transform it to the AO basis... Agrrrrhhh!
        if(iSpin.eq.1) then
          If ( IPRLEV.ge.DEBUG ) then
            write(6,*) 'running a singlet. DSPN set to zero!'
          end if
          DSPN(:)=0.0D0
        end if
        CALL mma_allocate(TmpDS ,NACPAR,Label='TmpDS')
        CALL mma_allocate(TmpD1S,NTOT2,Label='TmpD1S')
        Call dcopy_(NACPAR,DSPN,1,TmpDS,1)
        IF ( NASH(1).NE.NAC ) then
          CALL DBLOCK(TmpDS)
        end if
        Call Get_D1A_RASSCF(CMO,TmpDS,TmpD1S)
        CALL mma_deallocate(TmpDS)
        call CASDFT_terms(CMO,FOCK,FI,D1I,D1A,TmpD1S)
        CALL mma_deallocate(TmpD1S)
        CALL mma_deallocate(FOCK)
      end if

c        CALL TRIPRT('Averaged one-body density matrix, D, in RASSCF',
c     &              ' ',DMAT,NAC)
c        CALL TRIPRT('Averaged one-body spin density matrix DS, RASSCF',
c     &              ' ',DSPN,NAC)
c        CALL TRIPRT('Averaged two-body density matrix, P',
c     &              ' ',PMAT,NACPAR)
c        CALL TRIPRT('Averaged antisym 2-body density matrix PA RASSCF',
c     &              ' ',PA,NACPAR)

      Call Get_D1A_RASSCF(CMO,DMAT,D1A)
       If ( IPRLEV.ge.DEBUG ) then
        Write(LF,*)
        Write(LF,*) ' D1A in AO basis in RASSCF bf SXCTL'
        Write(LF,*) ' ---------------------'
        Write(LF,*)
        iOff=1
        Do iSym = 1,nSym
         iBas = nBas(iSym)
         call wrtmat(D1A(ioff),iBas,iBas, iBas, iBas)
         iOff = iOff + iBas*iBas
        End Do
       End If
      Call Timing(dum1,dum2,time3(2),dum3)
      TimeCIOpt = TimeCIOpt + time3(2) - time3(1)

*
c      Call rasscf_xml(Iter)
      Call rasscf_mcontrol(Iter)
*
* SX-section
*
      Call Timing(dum1,dum2,time2(1),dum3)

      If ( IPRLEV.ge.DEBUG ) then
       Write(LF,*) ' In RASSCF bf SXCTL'
       CALL TRIPRT('Averaged one-body density matrix, D, in RASSCF',
     &             ' ',DMAT,NAC)
       CALL TRIPRT('Averaged two-body density matrix, P',
     &             ' ',PMAT,NACPAR)
       CALL TRIPRT('Averaged antisym 2-body density matrix PA RASSCF',
     &             ' ',PA,NACPAR)
      end if
      CALL SXCTL(CMO,OCCN,DMAT,PMAT,PA,FI,FA,D1A,THMAX,IFINAL)


      If ( IPRLEV.ge.DEBUG ) then
       Write(LF,*)
       Write(LF,*) ' FI+FA in RASSCF after SXCTL'
       Write(LF,*) ' ---------------------'
       Write(LF,*)
       iOff=1
       Do iSym = 1,nSym
        iBas = nBas(iSym)
        Call TriPrt(' ',' ',FA(iOff),iBas)
        iOff = iOff + (iBas*iBas+iBas)/2
       End Do
      End If

cGLM   write(6,*) 'ECAS in RASSCF after call to SXCTL', ECAS
      If ( IPRLEV.ge.DEBUG ) then
        Write(LF,*)
        Write(LF,*) ' D1A in AO basis in RASSCF af SXCTL'
        Write(LF,*) ' ---------------------'
        Write(LF,*)
        iOff=1
        Do iSym = 1,nSym
         iBas = nBas(iSym)
         call wrtmat(D1A(ioff),iBas,iBas, iBas, iBas)
         iOff = iOff + iBas*iBas
        End Do
      end if
      Call Get_D1A_RASSCF(CMO,DMAT,D1A)
      CASDFT_Funct=0.0d0
      If (KSDFT.ne.'SCF'.and.KSDFT.ne.'PAM') then
        Call Get_dScalar('CASDFT energy',CASDFT_Funct)
cGLM        write(6,*) 'CASDFT energy :', CASDFT_Funct

      end IF
      DE=(ECAS+CASDFT_Funct)-ECAS1
      ECAS1=ECAS+CASDFT_Funct

      IF(NAC.EQ.0)   EAV=ECAS

      If (KSDFT(1:3).ne.'SCF'.and.KSDFT(1:3).ne.'PAM') Then
        IF (nConf.EQ.1) Then
          EAV = ECAS
        Else
          EAV = EAV - VIA_DFT - HALFQ
        End If
      End IF

      IF (ITER.EQ.1) DE=0.0d0
      Call Timing(dum1,dum2,time2(2),dum3)
      TimeOrb = TimeOrb + time2(2) - time2(1)
      TMXTOT=MAX(TMXTOT,THMAX)
*
* Save energies and convergence parameters
*
      CONV(1,ITER)=ECAS
      CONV(2,ITER)=ESX
      CONV(3,ITER)=CMAX
      CONV(4,ITER)=DE
      CONV(5,ITER)=CBLBM
      CONV(6,ITER)=ROTMAX
      IAD15 = IADR15(6)
      CALL DDAFILE(JOBIPH,1,ENER,mxRoot*mxIter,IAD15)
      CALL DDAFILE(JOBIPH,1,CONV,6*mxIter,IAD15)
#ifdef _HDF5_
      call mh5_put_attr(wfn_iter, Iter)
      call mh5_put_dset(wfn_energy, ENER(1,Iter))
#endif

*
* Print output of energies and convergence parameters
*
      Call Timing(dum1,dum2,time0(2),dum3)
      time0(2) = time0(2) - time0(1)
* Character indicating unconvergence/convergence criterion fulfilled:
      CTHRE=' '
      CTHRSX=' '
      CTHRTE=' '
      IF(ABS(DE).GT.THRE) CTHRE='*'
      IF(ABS(CBLBM).GT.THRSX) CTHRSX='*'
      IF(ABS(ROTMAX).GT.THRTE) CTHRTE='*'
      IROT=0
      IF(NROOTS.EQ.1) IROT=IROOT(1)
      IF ( IPRLEV.GE.2 .and. IPRLEV.LE.3 .AND..NOT.lOPTO) THEN
*----------------------------------
* Shift total energies (BOR 070411)
        if(iter.eq.1 .and. ICIONLY == 0)then
         kau=int(ECAS/1000.d0)
         EVAC=1000.d0*DBLE(kau)
         if(kau.ne.0) THEN
          write(6,'(6x,A,f23.2,A)')
     &       'Total energies have been shifted. Add ', EVAC,' au'
         end if
        endif
*----------------------------------
        ihh=int(time0(2)/3600)
        imm=int(time0(2)-ihh*3600)/60
        iss=int(time0(2)-ihh*3600-imm*60)
        if (DoSplitCAS) then
         Write(LF,'(6X,I3,I4,I5,I5,F15.8,ES12.2,A1,ES10.2,A1,2I4,I2,'//
     &            'ES10.2,A1,F6.2,F7.2,4X,A2,3X,A3,I5,A1,I2.2,A1,I2.2)')
     &        ITER,iterSplit,
     &        ITERSX,IROT,EAV,
     &        DE,CTHRE,ROTMAX,CTHRTE,IBLBM,JBLBM,ISYMBB,CBLBM,CTHRSX,
     &        SXSHFT,TMIN,QNSTEP,QNUPDT,ihh,':',imm,':',iss
        else if(DoBKAP) then
      Write(LF,'(3X,I3,I4,I2,I2,F15.8,F15.8,ES12.2,A1,ES10.2,A1,2I4,'//
     &         'I2,ES10.2,A1,F6.2,F7.2,4X,A2,3X,A3,I5,A1,I2.2,A1,I2.2)')
     &        ITER,ITERCI,
     &        ITERSX,IROT,ECAS-EVAC+CASDFT_Funct,EAV,DE,CTHRE,
     &        ROTMAX,CTHRTE,IBLBM,JBLBM,ISYMBB,CBLBM,CTHRSX,
     &        SXSHFT,TMIN,QNSTEP,QNUPDT,ihh,':',imm,':',iss

        else
          if(doDMRG.and.KeyCION)then ! If DMRG only -- yma
            Write(LF,'(/6X,a,F19.10)')
     &      'DMRGCI energy              :',EAV+CASDFT_Funct
            Write(LF,'(6X,a,I5,A1,I2.2,A1,I2.2/)')
     &      'Total time spent (hh:mm:ss):        ',ihh,':',imm,':',iss
          else
            if(doDMRG)then
#ifdef _DMRG_
              maxtrW = 0.0d0
              maxtrR = -1
              maxBD = -1
! These dmrg variables are arrays of rank 1
              ITERCI = MAXVAL(dmrg_energy%num_sweeps)
              IROT   = MAXLOC(dmrg_energy%num_sweeps, 1)
              maxtrW = MAXVAL(dmrg_energy%max_truncW)
              maxtrR = MAXLOC(dmrg_energy%max_truncW, 1)
              maxBD   = MAXVAL(dmrg_energy%bond_dim)
         Write(LF,'(6X,I3,I3,I4,I7,ES12.2,I4,I5,F15.8,ES12.2,A1,'//
     &            'ES9.2,A1,2I4,I2,ES10.2,A1,F6.2,F7.2,4X,A2,3X,A3,'//
     &            'I7,A1,I2.2,A1,I2.2)')
     &        ITER,ITERCI,IROT,maxBD,maxtrW,maxtrR,
     &        ITERSX,ECAS-EVAC+CASDFT_Funct,DE,CTHRE,
     &        ROTMAX,CTHRTE,IBLBM,JBLBM,ISYMBB,CBLBM,CTHRSX,
     &        SXSHFT,TMIN,QNSTEP,QNUPDT,ihh,':',imm,':',iss
#endif
            else
            Write(LF,'(6X,I3,I4,I5,I5,F15.8,ES12.2,A1,ES10.2,A1,2I4,'//
     &               'I2,ES10.2,A1,F6.2,F7.2,4X,A2,3X,A3,I5,A1,I2.2,'//
     &               'A1,I2.2)')
     &          ITER,ITERCI,
     &          ITERSX,IROT,ECAS-EVAC+CASDFT_Funct,DE,CTHRE,
     &          ROTMAX,CTHRTE,IBLBM,JBLBM,ISYMBB,CBLBM,CTHRSX,
     &          SXSHFT,TMIN,QNSTEP,QNUPDT,ihh,':',imm,':',iss
            end if
          end if
        end if
      ELSE IF ( IPRLEV.GE.4) THEN
        Write(LF,'(6X,A,I4)') 'Energy statistics and convergence'//
     &                                          ' in iteration',ITER
        Write(LF,'(6X,A,10X,F26.6)') 'Average of CI energies',EAV-EVAC
        Write(LF,'(6X,A,F26.6)') 'Complete active space SCF energy',
     &                                                       ECAS-EVAC
        Write(LF,'(6X,A,17X,F26.6)') 'Super-CI energy',ESX
        Write(LF,'(6X,A,12X,F26.6,3X,A1)') 'RASSCF energy change',
     &                                                        DE,CTHRE
        Write(LF,'(6X,A,I1,A1,2I3,F19.6)') 'Maximum BLB matrix'
     &                  //' element(sym=',ISYMBB,')',IBLBM,JBLBM,CBLBM
        Write(LF,'(6X,A,10X,F26.6,3X,A1)') 'Max rotation parameter',
     &                                               ROTMAX,CTHRSX
        Write(LF,'(6X,A,14X,F26.6,3X,A1)') 'Max rotation angle',
     &                                               THMAX,CTHRTE
        Write(LF,'(6X,A,3X,F26.6)') 'Max change in MO coefficients',CMAX
      END IF
#ifdef _FDE_
      ! Embedding
      if (embpot) then
       Eemb=embPotEneMODensities(D1I, D1A, embInt,
     &                nBas, nTot2, nSym)
       Write(LF,*)"E from embedding potential (<Psi|v_emb|Psi>): ",Eemb
      end if
#endif
      Call XFlush(6)
cGLM some additional printout for MC-PDFT

      If( l_casdft ) then
             CASDFT_E = ECAS-EVAC+CASDFT_Funct
             call print_mcpdft(CASDFT_E)
      end if
*
* Compare RASSCF energy and average CI energy
*
        DIFFE = ABS((ECAS-EAV)/ECAS)
      if(.not.DoSplitCAS) then
        if(DoBKAP) then
         if(iter.eq.1) then
          IF (DIFFE.GT.1.D-10 .AND. NROOTS.EQ.1) THEN
            Write(LF,'(6X,A)') repeat('=',120)
            Call WarningMessage(2,'Rasscf and CI energies will differ.')
            Write(LF,*)'This is the price you pay by the diagonal '
     &   //'approximation over the BB block in the SplitCAS method.'
            Write(LF,*)'The RASSCF energy might also diverge!'
            Write(LF,'(A)') repeat('#',80)
          END IF
         end if
        else if( l_casdft ) then
          Write(LF,'(6X,A)') repeat('=',80)
          Write(LF,'(10X,A)') 'This is a POST-SCF correction using a '
     & //'modified  Hamiltonian.'
          write(LF,'(10X,A)') 'The RASSCF energy has been corrected and'
     & //' it will differ from'
          write(LF,'(10X,A)') 'the preceding CI energy.'
          Write(LF,'(6X,A)') repeat('=',80)
        else
          IF(doDMRG)then

#ifdef _DMRG_DEBUGPRINT_
            write(lf,*) "DMRG-SCF energy    ",ECAS
            write(lf,*) "DMRG sweeped energy",EAV
#endif

            if(KeyCION)then
            else
              IF (DIFFE.GT.1.D-6 .AND. NROOTS.EQ.1) THEN
                Write(LF,'(6X,A)') repeat('=',120)
              Call WarningMessage(2,'DMRGSCF and DMRG energies differ.')
                Write(LF,'(6X,A,I11)')    'iteration           ',ITER
                Write(LF,'(6X,A,F22.10)') 'DMRGSCF energy      ',ECAS
                Write(LF,'(6X,A,F22.10)') 'DMRG energy         ',EAV
                Write(LF,'(6X,A,F22.10)') 'relative difference ',DIFFE
                Write(LF,*)'About this difference:'
                Write(LF,*)'1) If possible, consider a larger M value'
               Write(LF,*)'2) Severe convergence problems. Maybe active'
                Write(LF,*)'   space is unsuitable for this system?'
                Write(LF,'(6X,A)') repeat('=',120)
                IF (DIFFE.GT.5.D-04 .AND. NROOTS.EQ.1) THEN
                  Write(LF,*)
                  Write(LF,*)"Warning : "
                  Write(LF,*)" Relative difference is near unacceptable"
                  Write(LF,*)" If possible, consider a larger M value"
                  Write(LF,*)
                end if
              end if
            end if
          else
            DIFFETol = 1.D-10
#ifdef _ENABLE_DICE_SHCI_
            if (DoBlockDMRG) DIFFETol = 1.D-8
#endif
            IF (DIFFE.GT.DIFFETol .AND. NROOTS.EQ.1) THEN
              Write(LF,'(6X,A)') repeat('=',120)
              Call WarningMessage(2,'Rasscf and CI energies differ.')
              Write(LF,'(6X,A,I11)')    'iteration           ',ITER
              Write(LF,'(6X,A,F22.10)') 'RASSCF energy       ',ECAS
              Write(LF,'(6X,A,F22.10)') 'CI energy           ',EAV
              Write(LF,'(6X,A,F22.10)') 'relative difference ',DIFFE
              Write(LF,*)'Severe convergence problems. Maybe the active'
              Write(LF,*)'   space is unsuitable for this system?'
              Write(LF,'(6X,A)') repeat('=',120)
              IF(DIFFE.GT.1.D-04.AND.NROOTS.EQ.1.AND. .not.l_casdft)THEN
                Write(LF,*)
                Write(LF,'(6X,A)') 'The program has to stop !!!'
                Write(LF,*)
                Write(LF,'(6X,A)') repeat('=',120)
                Write(LF,*)
                ITERM=99
                GOTO 2000
              END IF
            end if
          END IF
        end if
      else
*          Write(LF,'(6X,A,F22.10)') 'Split-RASSCF energy    ',ECAS
        IF (DIFFE.GT.5.0D-03) THEN
          Write(LF,'(6X,A)') repeat('*',120)
          Write(LF,'(6X,A)') 'The Split-RASSCF and Split-CI '//
     &                       'energies differ !!!'
*          Write(LF,'(6X,A,I11)')    'iteration           ',ITER
          Write(LF,'(6X,A,F22.10)') 'Split-RASSCF energy    ',ECAS
          Write(LF,'(6X,A,F22.10)') 'Split-CI energy        ',EAV
          Write(LF,'(6X,A,F22.10)') 'Relative difference    ', DIFFE
          Write(LF,*)'     Smaller is the difference more realiable is',
     &               ' the result.'
          Write(LF,*)'     To make the difference smaller try',
     &               ' to select a bigger AA block or use firstOrder ',
     &               ' keyword.'
          Write(LF,'(6X,A)') repeat('*',120)
        END IF
      end if

      if (write_orb_per_iter .and. king()) then
        call copy_(real_path('RASORB'),
     &             real_path('ITERORB.'//str(actual_iter)))
#ifdef _HDF5_
        call copy_(real_path('RASWFN'),
     &             real_path('RASWFN.'//str(actual_iter)))

#endif
      end if

*
* Convergence check:
* check is done on largest BLB matrix
* element (CBLBM), on difference in
* average energy, DE and on maximum value of non-
* diagonal rotation matrix element.
*
************************************************************************
************************************************************************
* IF CIONLY calculation the convergence is skipped and goes to line 2000
************************************************************************
************************************************************************

      IF (IFINAL.EQ.1) GOTO 2000
      IF (DE.GT.1.0D0) THEN
        Call StatusLine('RASSCF: ','No convergence.')
        Write(LF,*)
        Write(LF,'(6X,A)') repeat('=',120)
        Call WarningMessage(2,'Rasscf energy diverges.')
        Write(LF,'(6X,A,I11)')    'iteration           ',ITER
        Write(LF,'(6X,A,F22.10)') 'RASSCF energy       ',ECAS
        Write(LF,'(6X,A,F22.10)') 'energy difference   ',DE
        Write(LF,'(6X,A)') repeat('=',120)
        Write(LF,*)
        Write(LF,'(6X,A)') '!!! The program was forced to stop !!!'
        Write(LF,*)
        Write(LF,'(6X,A)') repeat('=',120)
        Write(LF,*)
        ITERM=99
        GOTO 2000
      ENDIF
      IF(ITER.LT.MAXIT) THEN
        Write(STLNE2,'(A12,I3)')'Iteration ',ITER
        Call StatusLine('RASSCF converged: ',STLNE2)
        IF(ABS(DE).GT.THRE) GO TO 1000
        IF(ABS(CBLBM).GT.THRSX) GO TO 1000
        IF(ABS(ROTMAX).GT.THRTE) GO TO 1000
        IF(ITER.LE.3.AND.ICIONLY.EQ.0) GO TO 1000   ! 3->0 checking
        IF(IPRLEV.ge.TERSE) Write(LF,'(6X,A,I3,A)')
     &        'Convergence after',ITER,' iterations'
        If ( IPRLEV.ge.DEBUG ) then
         Write(LF,*)
         Write(LF,*) ' CMO in RASSCF after convergence printout'
         Write(LF,*) ' ---------------------'
         Write(LF,*)
         ioff=1
         Do iSym = 1,nSym
          iBas = nBas(iSym)
          if(iBas.ne.0) then
            write(6,*) 'Sym =', iSym
            do i= 1,iBas
              write(6,*) (CMO(ioff+iBas*(i-1)+j),j=0,iBas-1)
            end do
            iOff = iOff + (iBas*iBas)
          end if
         End Do
        End If
        IFINAL=1
        Call Add_Info('RASSCF_ITER',[DBLE(ITER)],1,8)
*    Call Add_Info('RASSCF_THMX',TMXTOT,1,5)
        GOTO 1000
      ELSE
        IF(IPRLEV.ge.TERSE) Write(LF,'(6X,A,I3,A)')
     &        'No convergence after',ITER,' iterations'
        Write(STLNE2,'(A12,I3)')'Iteration ',ITER
        Call StatusLine('RASSCF max iter: ',STLNE2)
        IFINAL=1
        ITERM=16
        GOTO 1000
      ENDIF
*                                                                      *
************************************************************************
*                                                                      *
*     Compute Final CI vectors
*                                                                      *
************************************************************************
*                                                                      *
 2000 IFINAL=2
      ICICH=0

      if (KeyWRMA) then
        call dump_fciqmc_mats(dmat=DMAT(:),
     &                        dspn=DSPN(:),
     &                        psmat=pmat(:),
     &                        pamat=pa(:))
      end if

************************************************************************
******************           Closing up MC-PDFT      *******************
************************************************************************


c Clean-close as much as you can the CASDFT stuff...
      if( l_casdft ) goto 2010

** IPT2 = 1 for OUTO, CANOnical keyword...
      IF(IPT2.EQ.1) THEN
        IAD=IADR15(9)
        CALL DDAFILE(JOBIPH,2,CMO,NTOT2,IAD)
      ELSE
        IAD=IADR15(2)
        CALL DDAFILE(JOBIPH,2,CMO,NTOT2,IAD)
      ENDIF
        If ( IPRLEV.ge.DEBUG ) then
         Write(LF,*)
         Write(LF,*) ' CMO in RASSCF after DDAFILE'
         Write(LF,*) ' ---------------------'
         Write(LF,*)
         ioff=1
         Do iSym = 1,nSym
          iBas = nBas(iSym)
          if(iBas.ne.0) then
            write(6,*) 'Sym =', iSym
            do i= 1,iBas
              write(6,*) (CMO(ioff+iBas*(i-1)+j),j=0,iBas-1)
            end do
            iOff = iOff + (iBas*iBas)
          end if
         End Do
        End If
      IF (NROOTS.GT.1) THEN
       Call StatusLine('RASSCF: ','Compute final CI vectors')
      ELSE
       Call StatusLine('RASSCF: ','Compute final CI vector')
      END IF
*
* Transform two-electron integrals
*
      Call Timing(dum1,dum2,time2(1),dum3)
      If (.not.DoCholesky .or. ALGO.eq.1) Then
         Call mma_allocate(PUVX,NFINT,Label='PUVX')
         PUVX(:)=0.0D0
      EndIf

      Call Get_D1I_RASSCF(CMO,D1I)

       IPR=0
       IF(IPRLOC(2).EQ.4) IPR=5
       IF(IPRLOC(2).EQ.5) IPR=10
       CALL TRACTL2(CMO,PUVX,TUVX,D1I,
     &              FI,D1A,FA,IPR,lSquare,ExFac)
*
       If (.not.DoCholesky .or. ALGO.eq.1) Then
          Call mma_deallocate(PUVX)
       EndIf

      Call Timing(dum1,dum2,time2(2),dum3)
      TimeTrans = TimeTrans + time2(2) - time2(1)
*
* CI-section (to obtain final wave function)
* 1st and 2nd order density matrix in MO basis
* for the gradients.
*
*
* Print header to file containing informations on CI iterations.
*
       Write(IterFile,*)
       Write(IterFile,'(20A4)') ('****',i=1,20)
       Write(IterFile,'(15X,A)') 'RASSCF iteration: Final'
*
      Call Timing(dum1,dum2,time3(1),dum3)

      if (allocated(CI_solver)) then
          call CI_solver%run(actual_iter=actual_iter,
     &                    ifinal=ifinal,
     &                    iroot=iroot,
     &                    weight=weight,
     &                    CMO=CMO(:),
     &                    DIAF=DIAF(:),
     &                    D1I_AO=D1I(:),
     &                    D1A_AO=D1A(:),
     &                    TUVX=tuvx(:),
     &                    F_IN=FI(:),
     &                    D1S_MO=DSPN(:),
     &                    DMAT=DMAT(:),
     &                    PSMAT=pmat(:),
     &                    PAMAT=pa(:))

#if defined (_ENABLE_BLOCK_DMRG_) || defined (_ENABLE_CHEMPS2_DMRG_) || defined (_ENABLE_DICE_SHCI_)
      else If(DoBlockDMRG) Then
        CALL DMRGCTL(CMO,DMAT,DSPN,PMAT,PA,FI,D1I,D1A,TUVX,IFINAL,1)
#endif
#ifdef _DMRG_
! Leon 27/11/2017: Skip the final CI iteration if we're using DMRGCI
! and CIOnly. It's enabled only for DMRGCI with QCMaquis now
! (to exclude potential side effects)
! but consider extending it to other cases!
      else if (doDMRG .and. ICIONLY/=0) then
        continue
#endif
      else
        CALL CICTL(CMO,DMAT,DSPN,PMAT,PA,FI,FA,D1I,D1A,TUVX,IFINAL)
      end if
      if (lRF .and. (iPCMRoot<=0 .or. DWSolv%DWZeta/=0.0d+00)) then
        IAD15 = IADR15(6)
        CALL DDAFILE(JOBIPH,1,ENER,mxRoot*mxIter,IAD15)
      end if

      EAV=0.0d0
      If(DoSplitCAS) then
        EAV = ENER(lRootSplit,ITER)
      else
       DO KROOT=1,NROOTS
         EAV=EAV+ENER(IROOT(KROOT),ITER)*WEIGHT(KROOT)
       END DO
      end if
      IF(NAC.EQ.0) EAV=ECAS
      IF(NCRVEC.gt.0) then
* Core shift has been used
        Call mma_deallocate(CRVEC)
        Call mma_deallocate(CRPROJ)
      END IF
      Call Timing(dum1,dum2,time3(2),dum3)
      TimeCIOpt = TimeCIOpt + time3(2) - time3(1)
      TimeWfn = TimeWfn + time3(2) - time1(1)
      time1(1) = time3(2)
*
* Calculation of natural orbitals. These orbitals are stored on
* JOBIPH in IADR15(12), followed by the occupation numbers.
* Usage: Only for one-electron properties
* Note: in an average calculation natural orbitals are obtained
* for all roots. Stored sequentially at IADR15(12) on JOBPIH.
*PAM2009 This is done in NATORB. Note that it places the resulting
* natural orbitals and occ nos into JOBIPH. Thus the arguments nr
* 5 and 6 are simply those for the last root that NATORB processed,
* and they should not be reused for anything after return from NATORB.
* NATORB args: Arg1 is current CMO coeffs, used in CI;
*  all the rest should be regarded as scratch.
*
      Call StatusLine('RASSCF: ','Compute natural orbitals')
      IPR=0
      IF(IPRLOC(6).EQ.4) IPR=5
      IF(IPRLOC(6).EQ.5) IPR=10
      Call mma_allocate(CMON,NTOT2,Label='CMON')
      Call mma_allocate(OCCX,NTOT,Label='OCCX')
      Call mma_allocate(scr1,MAX(NACPAR,NO2M),Label='scr1')
      Call mma_allocate(scr2,NO2M,Label='scr2')
      Call mma_allocate(SMAT,NTOT1,Label='SMAT')
*PAM2009 NATORB args: Arg1 is current CMO coeffs, used in CI;
*  all the rest should be regarded as scratch.
      CALL NATORB_RASSCF(CMO,scr1,scr2,SMAT,CMON,OCCX)
      Call mma_deallocate(scr1)
      Call mma_deallocate(scr2)
*PAM2009 Deallocate CMON, OCCX.
      Call mma_deallocate(CMON)
      Call mma_deallocate(OCCX)

*
* Compute transition density matrices
      If (KeyTDM) Then
#ifdef _HDF5_
         Call mma_allocate(Tmp,NConf,Label='Tmp')
         Call mma_allocate(VecL,NConf,Label='VecL')
         Call mma_allocate(VecR,NConf,Label='VecR')
         Call mma_allocate(kcnf,NACTEL,Label='kcnf')
         Call mma_allocate(Dtmp,NAC*NAC,Label='Dtmp')
         Call mma_allocate(DStmp,NAC*NAC,Label='DStmp')
         jDisk=IADR15(4)
         Call DDafile(JOBIPH,2,Tmp,nConf,jDisk)
         Do jRoot=2,lRoots
*           Read and reorder the left CI vector
            Call DDafile(JOBIPH,2,Tmp,nConf,jDisk)
            Call Reord2(NAC,NACTEL,STSYM,1,
     &                  CONF,CFTP,Tmp,VecL,kcnf)
            kDisk=IADR15(4)
            Do kRoot=1,jRoot-1
*              Read and reorder the right CI vector
               Call DDafile(JOBIPH,2,Tmp,nConf,kDisk)
               Call Reord2(NAC,NACTEL,STSYM,1,
     &                     CONF,CFTP,Tmp,VecR,kcnf)
*              Compute TDM and store in h5 file
               Call Lucia_Util('Densi',
     &                         CI_Vector=VecL(:),
     &                         RVec=VecR(:))
               idx=(jRoot-2)*(jRoot-1)/2+kRoot
               Call mh5_put_dset(wfn_transdens,Dtmp(1:NAC*NAC),
     &              [NAC,NAC,1], [0,0,idx-1])
               If (iSpin.gt.1)
     &         Call mh5_put_dset(wfn_transsdens,DStmp(1:NAC**2),
     &              [NAC,NAC,1], [0,0,idx-1])
            End Do
         End Do
         Call mma_deallocate(TMP)
         Call mma_deallocate(VecL)
         Call mma_deallocate(VecR)
         Call mma_deallocate(kcnf)
         Call mma_deallocate(Dtmp)
         Call mma_deallocate(DStmp)
#else
         Call WarningMessage(1,'HDF5 support disabled, '//
     &                         'TDM keyword ignored.')
#endif
      End If

      if (KeySSCR) then
        call spin_correlation_driver(orb_range_p, orb_range_q, iroot)
        call mma_deallocate(orb_range_p)
        call mma_deallocate(orb_range_q)
      end if

*
*****************************************************************
* Export all information relevant to geometry optimizations.
* Save also the reaction field operator.
      Call Timing(dum1,dum2,time2(1),dum3)
      If (iRlxRoot.eq.0) iRlxRoot=iRoot(1)
*
* Replace average occ Fock with occ Fock for state iRlxRoot
* and densities with the densities for state iRLXRoot
c      write(6,*) 'I am in RASSCF before call to PutRlx!'
      If ( ITERM.ne.99 ) Then
         Call mma_allocate(Dens,nTot1,Label='Dens')
         Call PutRlx(DMAT,DSPN,PMAT,Dens,CMO)
         Call Export1(IFINAL,CMO,DMAT,PMAT,Dens,FockOcc)
         Call mma_deallocate(Dens)
      End If
      Call Timing(dum1,dum2,time2(2),dum3)
      TimeRelax = time2(2) - time2(1)
*****************************************************************
*
      EMY=EMY+CASDFT_Funct
*                                                                      *
************************************************************************
*                                                                      *
*     Output section
*                                                                      *
************************************************************************
*                                                                      *

      Call StatusLine('RASSCF: ','Printing results')
      IF (IPRLEV.GE.USUAL .AND..NOT.lOPTO) THEN
        Write(LF,*)
        Write(LF,'(6X,A)') repeat('*',120)
        Write(LF,'(6X,A,118X,A)') '*','*'
        Write(LF,'(6X,A,52X,A,53X,A)')
     &        '*','Final results','*'
        Write(LF,'(6X,A,118X,A)') '*','*'
        Write(LF,'(6X,A)') repeat('*',120)
        Write(LF,*)
      END IF
#ifdef _FDE_
      ! Embedding
      if (embpot) then
       Eemb=embPotEneMODensities(D1I, D1A, embInt,
     &                nBas, nTot2, nSym)
       Write(LF,*) "Final energy from embedding potential: ", Eemb
       ! Write out ESP on grid if requested
       if (embWriteEsp) then
        Call Get_iScalar('Unique atoms',nNuc)
        Call embPotOutputMODensities(nNuc,nSym,D1I,D1A,nBas,nTot2)
       end if
      end if
#endif

      If ( ITERM.ne.99 ) THEN
       If (.not.DoSplitCAS) then
        CALL OUTCTL(CMO,OCCN,SMAT,lOPTO)
       else
        CALL OUTCTLSplit(CMO,OCCN,SMAT,lOPTO)
       end if
      End If

      Call mma_deallocate(SMAT)

*
* Write information for MOLDEN
*
c  i_root=0 gives state- and spin-averaged natural orbitals
c  i_root>0 gives natural spin orbitals for that root
      mroots=nroots
      if ((iSpin.eq.1).and.(nroots.eq.1)) mroots=0
      Do i_root=0,mroots
        Call Interf(i_root,FDIAG,1,0)
      End Do

* Create output orbital files:
      Call OrbFiles(JOBIPH,IPRLEV)

************************************************************************
******************           Closing up RASSCF       *******************
************************************************************************

2010   continue

      If (DoCholesky.and.ALGO.eq.2) Then
         Call mma_deallocate(Qmat)
      EndIf

*  Release  some memory allocations
      Call mma_deallocate(DIAF)
      Call mma_deallocate(FockOcc)
      Call mma_deallocate(FI)
      Call mma_deallocate(FA)
      Call mma_deallocate(D1I)
      Call mma_deallocate(D1A)
      Call mma_deallocate(OccN)
      Call mma_deallocate(CMO)
#ifdef _DMRG_
* Free RDMs for the reaction field reference root in QCMaquis calculations
      if (doDMRG.and.PCM_On()) then
        Call mma_deallocate(RF1)
        if (twordm_qcm) then
          Call mma_deallocate(RF2)
        end if
      end if
#endif
      if (lRF) call DWSol_final()
*     call DWSCF_final()

c deallocating TUVX memory...
      Call mma_deallocate(TUVX)
      Call mma_deallocate(DSPN)
      Call mma_deallocate(DMAT)
      Call mma_deallocate(PMAT)
      Call mma_deallocate(PA)
!Leon: The velociraptor comes! xkcd.com/292/
9989  Continue
*
* release SEWARD
      Call ClsSew
* ClsSew is needed for releasing memory used by integral_util, rys... which is allocated when MC-PDFT run is performed.

*---  Finalize Cholesky information if initialized
      if (DoCholesky)then
         Call Cho_X_Final(irc)
         if (irc.ne.0) then
            Write(LF,*)'RASSCF: Cho_X_Final fails with return code ',irc
            Write(LF,*)' Try to recover. Calculation continues.'
         endif
         If (Do_OFemb) Then
            Call mma_deallocate(FMaux)
            Call OFE_print(EAV)
         EndIf
      endif

!      do i=1,NTOT2
!        write(*,*) "A,I",D1A(i),D1I(i)
!      end do
      Call mma_deallocate(CleanMask,safe='*')

*
* Skip Lucia stuff if NECI or BLOCK-DMRG is on
      If (.not. any([allocated(CI_solver), DumpOnly,
     &              doDMRG, doBlockDMRG])) then
        Call Lucia_Util('CLOSE')
      end if


      Call StatusLine('RASSCF: ','Finished.')
      If (IPRLEV.GE.2) Write(LF,*)
      if(ifvb.eq.1) call make_close_rvb
cvv call to grid is moved up, in order to call clssew safely..
c      If (iCIonly.eq.0) Then
c        Call Grid_driver(-1,'RASSCF','RASORB',iR)
c      End If

      Call Timing(dum1,dum2,time1(2),dum3)
      TimeTotal = time1(2)
      TimeOutput = TimeOutput + time1(2) - time1(1)
      IF (IPRLEV.GE.3) THEN
       Call PrtTim
       Call FastIO('STATUS')
      END IF
      Call ClsFls_RASSCF()

*
c Rc_RAS  =  0 : The RASSCF wave function is converged
c         = 16 : The RASSCF wave function is not(!) converged
c         = 99 : The RASSCF energy is divergent or
c                the CI and SX energies differ
      Rc_RAS = ITERM
      Rc_RAS = Max(RC_RAS,Rc_CI)
      Rc_RAS = Max(RC_RAS,Rc_SX)
      If (Rc_Ras.eq.0) then
         ireturn=_RC_ALL_IS_WELL_
      Else If (Rc_Ras.eq.16) then
         ireturn=_RC_NOT_CONVERGED_
      Else
         Call WarningMessage(2,'Something is wrong: Did CI fail?')
         ireturn=_RC_GENERAL_ERROR_
      End If
*
      If (Do_OFemb) Then
         If (isStructure().eq.1) Then
            If (iReturn.ne._RC_ALL_IS_WELL_) Then
               Call WarningMessage(1,'RASSCF: non-zero return code.')
            EndIf
            iReturn=_RC_CONTINUE_LOOP_
            Call Check_FThaw(iReturn)
         EndIf
      EndIf

      if (.not. (iDoGas .or. doDMRG .or. doBlockDMRG
     &          .or. allocated(CI_solver) .or. DumpOnly)) then
        Call MKGUGA_FREE(SGS,CIS,EXS)
      end if

      if (DoFaro) then
         call faroald_free()
         call citrans_free()
      end if

      if (allocated(CI_solver)) then
          call CI_solver%cleanup()
          deallocate(CI_solver)
      end if

* DMRG: Save results for other use
! ==========================================================
      if(doDMRG)then
#ifdef _DMRG_
        !Leon: Generate 4-RDM evaluation templates for NEVPT2

        !In the new interface the EvRDM keyword will be ignored.
        !Instead, NEVPT2Prep will always generate the template.
        !RDM evaluation will now happen in the NEVPT2 module
        !where NEVPT2 either can attempt to compute it directly
        !with the new interface or read from QCMaquis HDF5 result
        !file.
        if (DoNEVPT2Prep) then
          if (MAXIT.eq.0) then
            write(6,*) " --- DMRG-SCF iterations are skipped, only "//
     &        "QCMaquis input for higher-order RDMs will be generated."
          endif
          if (NACTEL.gt.3) then ! Ignore 4-RDM if we have <4 electrons
          do i=1,NROOTS
              Write (6,'(a)') 'Writing 4-RDM QCMaquis template'//
     &   ' for state '//str(i)
              call qcmaquis_interface_prepare_hirdm_template(
     &        filename="meas-4rdm."//str(i-1)//".in",
     &        state=i-1,
     &        tpl=TEMPLATE_4RDM)
              call qcmaquis_mpssi_transform(
     &             trim(qcmaquis_param%workdir)//'/'//
     &             trim(qcmaquis_param%project_name), i)
          end do
          else
            write(6,*) "Skipping 4-RDM QCMaquis template generation "//
     &        "since we have less than 4 electrons."
          end if
          ! Generate 3-TDM templates
          if(NACTEL.gt.2) then ! but only if we have more than 3 el.
          do i=1,NROOTS
            do j=i+1,NROOTS
              Write (6,'(a)') 'Writing 3-TDM QCMaquis template'//
     &   ' for states '//str(i)//" and "//str(j)
              call qcmaquis_interface_prepare_hirdm_template(
     &        filename="meas-3tdm."//str(i-1)//"."//str(j-1)//".in",
     &        state=i-1,
     &        state_j=j-1,
     &        tpl=TEMPLATE_TRANSITION_3RDM)
            end do
          end do
          else
            write(6,*) "Skipping 3-RDM QCMaquis template generation "//
     &        "since we have less than 3 electrons."
          end if
        end if
        ! is it really needed in the times of Fortran 2008?
        call qcmaquis_interface_deinit
#endif
      end if
! ==========================================================
! Exit

 9990 Continue
C Close the one-electron integral file:
      iRC=-1
      iOpt=0
      Call ClsOne(iRC,iOpt)
      If (iRC.ne.0) Then
         Write (6,*) 'Error when trying to close the one-electron'
         Write (6,*) 'integral file.'
         Call Quit(_RC_INTERNAL_ERROR_)
      End If

      If (IfVB.ne.2) Then
        DO I=10,99
          INQUIRE(UNIT=I,OPENED=IfOpened)
          IF (IfOpened.and.I.ne.19) CLOSE (I)
        END DO
        Close(LUInput)
      End If

      return

      end subroutine rasscf
