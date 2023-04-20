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
*               2016, Andrew M. Sand                                   *
************************************************************************
      SUBROUTINE MCPDFT(IRETURN)
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
*                 A program for MC-PDFT calculations                   *
*                 Called after RASSCF is called.                       *
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
*     March '98.                                                       *
*                                                                      *
*     Modified AMS Feb 2016 - separate MCPDFT from RASSCF              *
************************************************************************

      use csfbas, only: CONF, KCFTP
      use hybridpdft, only: do_hybrid
      use Fock_util_global, only: ALGO, DoActive, DoCholesky
      use OFembed, only: Do_OFemb, FMaux
      use UnixInfo, only: ProgName
      use stdalloc, only : mma_allocate, mma_deallocate
      use write_pdft_job, only: iwjob, writejob
      use sxci_pdft, only: idxsx
      use mspdft, only: dogradmspd, mspdftmethod, do_rotate, iF1MS,
     &                  iF2MS, iFxyMS, iFocMS, iDIDA, IP2MOt, D1AOMS,
     &                  D1SAOMS, doNACMSPD, cmsNACstates, doMECIMSPD
      Implicit Real*8 (A-H,O-Z)

#include "WrkSpc.fh"
#include "wadr.fh"
#include "rasdim.fh"
#include "warnings.h"
#include "input_ras_mcpdft.fh"
#include "rasscf.fh"
#include "general.fh"
#include "gas.fh"
#include "bk_approx.fh"
#include "output_ras.fh"
#include "rctfld.fh"
#include "timers.fh"
#include "casvb.fh"
#include "rasscf_lucia.fh"
#include "lucia_ini.fh"
#include "gugx.fh"
#include "pamint.fh"
#include "ciinfo.fh"
      Integer LRState,NRState         ! storing info in Do_Rotate.txt
      Integer LHrot,NHrot             ! storing info in H0_Rotate.txt
      Real*8  MSPDFTShift
      Logical lshiftdiag
      CHARACTER(Len=18)::MatInfo
      Integer LXScratch,NXScratch
      INTEGER LUMS,IsFreeUnit
      Dimension WGRONK(2)
      External IsFreeUnit

      Logical DSCF
      Character*80 Line
      Logical IfOpened
      Logical Found
      Character(len=9),DIMENSION(:),Allocatable::VecStat
      CHARACTER(Len=9)::StatVec
      CHARACTER(Len=30)::mspdftfmt
      Logical RefBas

* --------- Cholesky stuff:
#include "chopar.fh"
#include "chotime.fh"
* --------- End Cholesky stuff
      Character*8 EMILOOP

      External RasScf_Init_m
      External Scan_Inp_m
      integer iRef_E,IAD19
      integer IADR19(1:15)
      integer NMAYBE,KROOT
      real*8 EAV
!
      real*8, allocatable :: PLWO(:)
      integer ivkcnf
      Dimension Dummy(1)
* Set status line for monitor:
      Call StatusLine('MCPDFT:',' Just started.')
* Set the return code(s)
      IRETURN=_RC_ALL_IS_WELL_

* Local print level in this routine:
      IPRLEV=IPRLOC(1)

* Set some Cholesky stuff
      DoActive=.true.

* Set variable IfVB to check if this is a VB job.
      IfVB=0
      If (ProgName(1:5).eq.'casvb') IfVB=2
* Default option switches and values, and initial data.
      EAV = 0.0d0
      Call RasScf_Init_m()
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


      If (IfVB /= 2) then

* Make a copy, upper-cased, left-adjusted, of the input between and including
* the '&MCPDFT' and the 'End of input' markers, skipping all lines beginning
* with '*' or '!' or ' '  when left-adjusted, and replacing any rightmost
* substring beginning with '!' with blanks.
* That copy will be in file 'CleanInput', and its unit number is returned
* as LUInput in common (included file input_ras_mcpdft.fh) by the following call:
        Call cpinp_(LUInput,iRc)
* If something wrong with input file:
        If (iRc.ne._RC_ALL_IS_WELL_) Then
         Call WarningMessage(2,'Input file is unusable.')
         Write(6,*)' MCPDFT Error: Could not make a clean copy of'
         Write(6,*)' the input file. This is an unexpected bug.'
         IRETURN=_RC_INTERNAL_ERROR_
         GOTO 9990
        End If


! Scan the input file for keywords:
        Call Scan_Inp_m(iRc)
! If something wrong with input file:
        If (iRc.ne._RC_ALL_IS_WELL_) Then
         If (IPRLOC(1).GE.TERSE) Then
          Call WarningMessage(2,'Scanning input file failed.')
! Calling again, now with iRc indicating an error, will echo the keywords:
          Call Scan_Inp_m(iRc)
         End If
         IRETURN=_RC_INPUT_ERROR_
         GOTO 9990
        End If

! Local print level in this routine:
        IPRLEV=IPRLOC(1)
*
      end if
* Open files
      Call OpnFls_RASSCF_m(DSCF,DoCholesky)

* Some preliminary input data:
      Call Rd1Int_m()
      If ( .not.DSCF ) Call Rd2Int_MCPDFT

* Process the input:
      Call Proc_InpX(DSCF,iRc)
* If something goes wrong in proc_inp:
      If (iRc.ne._RC_ALL_IS_WELL_) Then
        If (IPRLEV.ge.TERSE) Then
          Call WarningMessage(2,'Input processing failed.')
          Write(6,*)' RASSCF Error: Proc_Inp failed unexpectedly.'
          Write(6,*)' Here is a printing of the input file that'
          Write(6,*)' was processed:'
          Rewind(LUInput)
  15      Continue
          Read(LuInput,'(A80)',End=16,Err=16) Line
          Write(6,*) Line
          Go To 15
  16      Continue
        End If
        IRETURN=iRc
        GOTO 9990
      End If


* Local print level may have changed:
      IPRLEV=IPRLOC(1)


      Call InpPri_m()

*--------------------------------------------------------
*
* Allocate various matrices
*
      Call GetMem('FI','Allo','Real',LFI,NTOT1)
      Call GetMem('FA','Allo','Real',LFA,NTOT1)
      Call GetMem('D1I','Allo','Real',LD1I,NTOT2)
      Call GetMem('D1A','Allo','Real',LD1A,NTOT2)
      Call GetMem('D1tot','Allo','Real',LD1tot,NTOT1)
      Call GetMem('OCCN','Allo','Real',LOCCN,NTOT)
      Call GetMem('LCMO','Allo','Real',LCMO,NTOT2)
      allocate(PLWO(1:NACPAR))
      PLWO(:) = 0
!
*
      LTUVX=1
      LDMAT=1
      LDSPN=1
      LPMAT=1
      LPA  =1
      If ( NAC.GT.0 ) then

        Call GetMem('TUVX','Allo','Real',LTUVX,NACPR2)
        Call FZero(Work(LTUVX),NACPR2)
        ltuvx_cvb=ltuvx
        Call GetMem('DMAT','Allo','Real',LDMAT,NACPAR)
        Call GetMem('DSPN','Allo','Real',LDSPN,NACPAR)
        Call GetMem('PMAT','Allo','Real',LPMAT,NACPR2)
        Call GetMem('P2AS','Allo','Real',LPA,NACPR2)
        call dcopy_(NACPAR,[0.0d0],0,Work(LDMAT),1)
        call dcopy_(NACPAR,[0.0d0],0,Work(LDSPN),1)
      Else
        LTUVX = ip_Dummy
        ltuvx_cvb=ltuvx
        LDMAT = ip_Dummy
        LDSPN = ip_Dummy
        LPMAT = ip_Dummy
        LPA   = ip_Dummy
      End If
*
* Get start orbitals

* Initialize OCCN array, to prevent false alarms later from
* automated detection of using uninitialized variables:
      call dcopy_(NTot,[0.0D0],0,Work(lOCCN),1)

* PAM03: Note that removal of linear dependence may change the nr
* of secondary/deleted orbitals, affecting some of the global
* variables: NSSH(),NDEL(),NORB(),NTOT3, etc etc
      Call ReadVc_m(Work(LCMO),Work(lOCCN),
     &             WORK(LDMAT),WORK(LDSPN),WORK(LPMAT),WORK(LPA))
* Only now are such variables finally known.
      If (IPRLOC(1).GE.DEBUG) Then
        CALL TRIPRT('Averaged one-body density matrix, D, in RASSCF',
     &              ' ',Work(LDMAT),NAC)
        CALL TRIPRT('Averaged one-body spin density matrix DS, RASSCF',
     &              ' ',Work(LDSPN),NAC)
        CALL TRIPRT('Averaged two-body density matrix, P',
     &              ' ',WORK(LPMAT),NACPAR)
        CALL TRIPRT('Averaged antisym 2-body density matrix PA RASSCF',
     &              ' ',WORK(LPA),NACPAR)
      END IF
*
* Allocate core space for dynamic storage of data
*
      CALL ALLOC_m()

      Call Timing(Swatch,Swatch,Ebel_1,Swatch)

      ECAS   = 0.0d0
      Call GetMem('FOcc','ALLO','REAL',ipFocc,nTot1)

      if (KSDFT(1:5).eq.'T:'.or. KSDFT(1:3).eq.'FT:') Then
        KSDFT_TEMP=KSDFT
        KSDFT='SCF'
        ExFac=1.0D0
      else
        KSDFT_TEMP=KSDFT
      end if

      If (iCIRST.eq.1) Then
        Call GetMem('TmpDMAT','Allo','Real',ipTmpDMAT,NACPAR)
        call dcopy_(NACPAR,Work(LDMAT),1,Work(ipTmpDMAT),1)
        If (NASH(1).ne.NAC) Call DBLOCK_m(Work(ipTmpDMAT))
        Call Get_D1A_RASSCF_m(Work(LCMO),Work(ipTmpDMAT),WORK(LD1A))
        Call GetMem('TmpDMAT','Free','Real',ipTmpDMAT,NACPAR)
        DoActive = .true.
      Else
        lRf = .false.
        IF(KSDFT_TEMP(1:2).ne.'T:') Then
          KSDFT='SCF'
          ExFac=1.0D0
        end IF
        Call dcopy_(NTOT2,[0.0D0],0,WORK(LD1A),1)
        DoActive = .false.
      End If

! 413 Continue

!AMS start-
! - Read in the CASSCF Energy from JOBIPH file.  These values are not
! used in calculations, but are merely reprinted as the reference energy
! for each calculated MC-PDFT energy.
!

      iJOB=0
      Call GetMem('REF_E','ALLO','REAL',iRef_E,lroots)
      Call Fzero(Work(iRef_E),lroots)
        Call f_Inquire('JOBOLD',Found)
        if (.not.found) then
          Call f_Inquire('JOBIPH',Found)
          if(Found) JOBOLD=JOBIPH
        end if
        If (Found) iJOB=1
        If (iJOB.eq.1) Then
           if(JOBOLD.le.0) Then
             JOBOLD=20
             Call DaName(JOBOLD,'JOBOLD')
           end if
        end if
       IADR19(:)=0
       IAD19=0
      Call IDaFile(JOBOLD,2,IADR19,15,IAD19)
      jdisk = IADR19(6)
!I must read from the 'old' JOBIPH file.
      Call GetMem('ELIST','ALLO','REAL',iEList,MXROOT*MXITER)
      Call DDaFile(JOBOLD,2,Work(iEList),MXROOT*MXITER,jdisk)
      NMAYBE=0
      DO IT=1,MXITER
        AEMAX=0.0D0
        DO I=1,MXROOT
          E=WORK(iEList+MXROOT*(IT-1)+(I-1))
          AEMAX=MAX(AEMAX,ABS(E))
        END DO
        IF(ABS(AEMAX).LE.1.0D-12) GOTO 11
        NMAYBE=IT
      END DO
  11  CONTINUE
      Do_Rotate=.false.
      IF(iMSPDFT==1) Then
       call f_inquire('ROT_HAM',Do_Rotate)
       If(.not.Do_Rotate) Then
        write(6,'(6X,A,A)')'keyword "MSPD" is used but ',
     &  'the file of rotated Hamiltonian is not found.'
        write(6,'(6X,2a)')'Performing regular (state-',
     &   'specific) MC-PDFT calculation'
       End If
      End IF
      IF(Do_Rotate) Then
        write(6,'(6X,80A)') ('=',i=1,80)
        write(6,*)
        write(6,'(6X,A,A)')'keyword "MSPD" is used and ',
     &  'file recording rotated hamiltonian is found. '
        write(6,*)
        write(6,'(6X,A,A)')
     &  'Switching calculation to Multi-State Pair-Density ',
     &  'Functional Theory (MS-PDFT) '
        write(6,'(6X,A)')'calculation.'
        write(6,*)
        NHRot=lroots**2
        CALL GETMEM('HRot','ALLO','REAL',LHRot,NHRot)
        LUMS=12
        LUMS=IsFreeUnit(LUMS)
        CALL Molcas_Open(LUMS,'ROT_HAM')
        Do Jroot=1,lroots
          read(LUMS,*) (Work(LHRot+Jroot-1+(Kroot-1)*lroots)
     &                 ,kroot=1,lroots)
        End Do
        Read(LUMS,'(A18)') MatInfo
        MSPDFTMethod=' MS-PDFT'
        IF(trim(adjustl(MatInfo)).eq.'an unknown method') THEN
         write(6,'(6X,A,A)')'The MS-PDFT calculation is ',
     & 'based on a user-supplied rotation matrix.'
        ELSE
         write(6,'(6X,A,A,A)')'The MS-PDFT method is ',
     &   trim(adjustl(MatInfo)),'.'
        If(trim(adjustl(MatInfo)).eq.'XMS-PDFT') MSPDFTMethod='XMS-PDFT'
        If(trim(adjustl(MatInfo)).eq.'CMS-PDFT') MSPDFTMethod='CMS-PDFT'
        If(trim(adjustl(MatInfo)).eq.'VMS-PDFT') MSPDFTMethod='VMS-PDFT'
        If(trim(adjustl(MatInfo)).eq.'FMS-PDFT') MSPDFTMethod='FMS-PDFT'
        ENDIF
        write(6,*)
        write(6,'(6X,80A)') ('=',i=1,80)
        write(6,*)
        Close(LUMS)
        do KROOT=1,lROOTS
          ENER(IROOT(KROOT),1)=Work((LHRot+(Kroot-1)*lroots+
     &                                     (KROOT-1)))
           EAV = EAV + ENER(IROOT(KROOT),ITER) * WEIGHT(KROOT)
           Work(iRef_E + KROOT-1) = ENER(IROOT(KROOT),1)
        end do
      Else
        do KROOT=1,lROOTS
          ENER(IROOT(KROOT),1)=Work(iEList+MXROOT*(NMAYBE-1) +
     &                                     (KROOT-1))
           EAV = EAV + ENER(IROOT(KROOT),ITER) * WEIGHT(KROOT)
           Work(iRef_E + KROOT-1) = ENER(IROOT(KROOT),1)
        end do
      End IF!End IF for Do_Rotate=.true.

      IF(doNACMSPD) Then
        write(6,'(6X,80A)') ('=',i=1,80)
        write(6,*)
        write(6,'(6X,A,I3,I3)')'keyword NAC is used for states:',
     & cmsNACstates(1), cmsNACstates(2)
        write(6,*)
        write(6,'(6X,80A)') ('=',i=1,80)
        call Put_lScalar('isCMSNAC        ', doNACMSPD)
        call Put_iArray('cmsNACstates    ', cmsNACstates, 2)
      ELSE
        cmsNACstates(1) = iRlxRoot
        cmsNACstates(2) = 0
        call Put_lScalar('isCMSNAC        ', doNACMSPD)
        call Put_iArray('cmsNACstates    ', cmsNACstates, 2)
      End IF!End IF for doNACMSPD=.true.

      IF(doMECIMSPD) Then
        write(6,'(6X,80A)') ('=',i=1,80)
        write(6,*)
        write(6,'(6X,A,I3,I3)')'keyword MECI is used for states:'
        write(6,*)
        write(6,'(6X,80A)') ('=',i=1,80)
        call Put_lScalar('isMECIMSPD      ', doMECIMSPD)
      ELSE
        call Put_lScalar('isMECIMSPD      ', doMECIMSPD)
      End IF

      Call GetMem('ELIST','FREE','REAL',iEList,MXROOT*MXITER)
      If(JOBOLD.gt.0.and.JOBOLD.ne.JOBIPH) Then
        Call DaClos(JOBOLD)
        JOBOLD=-1
      End if

*
      IF (KSDFT_TEMP(1:2).eq.'T:'.or.KSDFT_TEMP(1:3).eq.'FT:') Then
        KSDFT=KSDFT_TEMP
        ExFac=0.0d0
      end IF
*
* Transform two-electron integrals and compute at the same time
* the Fock matrices FI and FA
*
      Call Timing(Swatch,Swatch,Fortis_1,Swatch)
      If (.not.DoCholesky .or. ALGO.eq.1) Then
         Call GetMem('PUVX','Allo','Real',LPUVX,NFINT)
         Call FZero(Work(LPUVX),NFINT)
      Else
         LPUVX=ip_Dummy
      EndIf
      Call Get_D1I_RASSCF_m(Work(LCMO),Work(lD1I))

      DoActive = .true.

      IPR=0
      IF(IPRLOC(2).EQ.debug) IPR=5
      IF(IPRLOC(2).EQ.insane) IPR=10

      CALL TRACTL2(WORK(LCMO),WORK(LPUVX),WORK(LTUVX),WORK(LD1I),
     &              WORK(LFI),WORK(LD1A),WORK(LFA),IPR,lSquare,ExFac)

      Call Put_dArray('Last orbitals',Work(LCMO),ntot2)

      if (doGSOR) then
        Call f_Inquire('JOBOLD',Found)
        if (.not.found) then
          Call f_Inquire('JOBIPH',Found)
          if(Found) JOBOLD=JOBIPH
        end if
        If (Found) iJOB=1
        If (iJOB.eq.1) Then
           if(JOBOLD.le.0) Then
             JOBOLD=20
             Call DaName(JOBOLD,'JOBOLD')
           end if
        end if
        IADR19(:)=0
        IAD19=0
        LUCT=87
        LUCT=IsFreeUnit(LUCT)
        CALL Molcas_Open(LUCT,'CI_THETA')

        Call IDaFile(JOBOLD,2,IADR19,15,IAD19)
        CALL GETMEM('CIVEC','ALLO','REAL',LW4,NCONF)
        CALL GETMEM('Dtmp ','ALLO','REAL',LW6,NACPAR)
        CALL GETMEM('DStmp','ALLO','REAL',LW7,NACPAR)
        CALL GETMEM('Ptmp ','ALLO','REAL',LW8,NACPR2)
        CALL GETMEM('PAtmp','ALLO','REAL',LW9,NACPR2)
        CALL GETMEM('Pscr','ALLO','REAL',LW10,NACPR2)

        call dcopy_(NACPAR,[0.0D0],0,WORK(LW6),1)
        call dcopy_(NACPAR,[0.0D0],0,WORK(LW7),1)
        call dcopy_(NACPR2,[0.0D0],0,WORK(LW8),1)
        call dcopy_(NCONF,[0.0D0],0,WORK(LW4),1)
        iDisk = IADR19(4)
        jDisk = IADR19(3)

        Call GetMem('CIVtmp','Allo','Real',LW11,nConf)
        DO jRoot=1,lroots
          do i=1,nconf
            read(LUCT,*) Work(LW4-1+i)
          end do
          Call DDafile(JOBOLD,1,Work(LW4),nConf,iDisk)
          call getmem('kcnf','allo','inte',ivkcnf,nactel)
          Call Reord2(NAC,NACTEL,STSYM,1,
     &                CONF,iWork(KCFTP),
     &                Work(LW4),Work(LW11),iWork(ivkcnf))
          Call dcopy_(nconf,Work(LW11),1,Work(LW4),1)
          call getmem('kcnf','free','inte',ivkcnf,nactel)
          C_Pointer = Lw4
          CALL GetMem('Lucia','Allo','Real',Lucia_Base, 1)
!Andrew - changed here
          CALL Lucia_Util('Densi',ip_Dummy,iDummy,Dummy)
          If (IFCAS.GT.2 .OR. iDoGAS) Then
            Call CISX_m(IDXSX,Work(LW6),Work(LW7),Work(LW8),
     &              Work(LW9),Work(LW10))
          End If
          CALL GetMem('Lucia','Free','Real',Lucia_Base, 1)

          Call DDafile(JOBOLD,1,Work(LW6),NACPAR,jDisk)
          Call DDafile(JOBOLD,1,Work(LW7),NACPAR,jDisk)
          Call DDafile(JOBOLD,1,Work(LW8),NACPR2,jDisk)
          Call DDafile(JOBOLD,1,Work(LW9),NACPR2,jDisk)
        end do
        Close(LUCT)
        Call fCopy('JOBIPH','JOBGS',ierr)

       end if!DoGSOR



      If (.not.DoCholesky .or. ALGO.eq.1) Then
         if(dogradmspd) then
           CALL Put_dArray('TwoEIntegral    ',Work(LPUVX),nFINT)
         end if
         Call GetMem('PUVX','Free','Real',LPUVX,NFINT)
      EndIf

      Call Timing(Swatch,Swatch,Fortis_2,Swatch)
      Fortis_2 = Fortis_2 - Fortis_1
      Fortis_3 = Fortis_3 + Fortis_2

      IF(KSDFT_TEMP(1:2).eq.'T:'.or. KSDFT_TEMP(1:3).eq.'FT:') Then
        IF(DoGradMSPD) THEN
          Call GetMem('F1MS' ,'Allo','Real',iF1MS ,nTot1*nRoots)
          Call GetMem('FocMS','Allo','Real',iFocMS,nTot1*nRoots)
          Call GetMem('FxyMS','Allo','Real',iFxyMS,nTot4*nRoots)
          Call GetMem('F2MS' ,'Allo','Real',iF2MS ,nACPR2*nRoots)
          Call GetMem('P2MO' ,'Allo','Real',iP2MOt,nACPR2*nRoots)
          Call GetMem('DIDA' ,'Allo','Real',iDIDA ,nTot1*(nRoots+1))
          Call GetMem('D1AOMS' ,'Allo','Real',D1AOMS,nTot1*nRoots)
          if (ispin.ne.1) then
            Call GetMem('D1SAOMS' ,'Allo','Real',D1SAOMS,nTot1*nRoots)
          end if
          Call FZero(Work(iP2MOt),lRoots*NACPR2)
        END IF
        CALL GETMEM('CASDFT_Fock','ALLO','REAL',LFOCK,NACPAR)

        Call MSCtl(Work(LCMO),Work(LFOCK),Work(LFI),Work(LFA),
     &       Work(iRef_E))

        If(IWJOB==1.and.(.not.Do_Rotate)) Call writejob(iadr19)

        If (Do_Rotate) Then
          NHRot=lroots**2
          Do Jroot=1,lroots
            Work(LHRot+Jroot-1+(Jroot-1)*lroots)=Work(iRef_E-1+Jroot)
          End DO
          Write(6,'(6X,80a)') ('*',i=1,80)
          Write(6,*)
          Write(6,'(34X,2A)')MSPDFTMethod,' FINAL RESULTS'
          Write(6,*)
          Write(6,'(6X,80a)') ('*',i=1,80)
          Write(6,*)

          lshiftdiag=.false.
          CALL shiftdiag(WORK(LHRot),MSPDFTShift,lshiftdiag,lRoots,10)
          if(.not.do_hybrid) then
            write(6,'(6X,2A)') MSPDFTMethod,' Effective Hamiltonian'
          else
            write(6,'(6X,3A)')
     &         'Hybrid ',MSPDFTMethod,' Effective Hamiltonian'
          end if
          if(lshiftdiag) then
            write(6,'(6X,A,F9.2,A)')
     &          '(diagonal values increased by',-MSPDFTShift,' hartree)'
            Do JRoot=1,lRoots
              Work(LHRot+Jroot-1+(Jroot-1)*lroots)=
     &                  Work(LHRot+Jroot-1+(Jroot-1)*lroots)-MSPDFTShift
            End Do
          end if
          Call RecPrt(' ','(7X,10(F9.6,1X))',Work(LHRot),lroots,lroots)
          write (6,*)
*MS-PDFT    To diagonalize the final MS-PDFT effective H matrix.
*MS-PDFT    Eigenvectors will be stored in LRState. This notation for the
*MS-PDFT    address here is the same for the rotated space in XMS-CASPT2.
          NRState=NHRot
          CALL GETMEM('RotStat','ALLO','REAL',LRState,NRState)
          Call FZero(Work(LRState),NRState)
          Call Dsyev_('V','U',lroots,Work(LHRot),lroots,Work(LRState),
     &               WGRONK,-1,INFO)
          NXScratch=Int(WGRONK(1))
          Call GetMem('XScratch','Allo','Real',LXScratch,NXScratch)
          Call Dsyev_('V','U',lroots,Work(LHRot),lroots,Work(LRState),
     &               Work(LXScratch),NXScratch,INFO)

          if(lshiftdiag) then
            Do Jroot=1,lRoots
              Work(LRState+Jroot-1)=Work(LRState+Jroot-1)+MSPDFTShift
            End Do
          end if

          if(.not.do_hybrid) then
            write(6,'(6X,2A)')MSPDFTMethod,' Energies:'
            Do Jroot=1,lroots
              write(6,'(6X,3A,1X,I4,3X,A13,F18.8)')
     &            '::    ',MSPDFTMethod,' Root',
     &            Jroot,'Total energy:',Work(LRState+Jroot-1)
            End Do
          else
            write(6,'(6X,3A)')'Hybrid ',MSPDFTMethod,' Energies:'
            Do Jroot=1,lroots
              write(6,'(6X,4A,1X,I4,3X,A13,F18.8)')
     &             '::    ','Hybrid ',MSPDFTMethod,' Root',
     &              Jroot,'Total energy:',Work(LRState+Jroot-1)
            End Do
          end if
          Call Put_iScalar('Number of roots',nroots)
          Call Put_dArray('Last energies',WORK(LRState),nroots)
          Call Put_dScalar('Last energy',WORK(LRState+iRlxRoot-1))
          Write(6,*)
          CALL mma_allocate(VecStat,lRoots)
          Do Jroot=1,lRoots
            write(StatVec,'(A5,I4)')'Root ',JRoot
            VecStat(JRoot)=StatVec
          End Do
          if(.not.do_hybrid) then
            write(6,'(6X,2A)')MSPDFTMethod,' Eigenvectors:'
          else
            write(6,'(6X,3A)')'Hybrid ',MSPDFTMethod,' Eigenvectors:'
          end if
          write(6,'(7X,A)')'Intermediate-state Basis'
          write(mspdftfmt,'(A4,I5,A9)')
     &           '(6X,',lRoots,'(A10,5X))'
          write(6,mspdftfmt)((VecStat(JRoot)),JRoot=1,lroots)
*Added by Chen to write energies and states of MS-PDFT into JOBIPH
          If(IWJOB==1) Call writejob(iadr19,LREnergy=LRState,LRot=LHRot)
          Call RecPrt(' ','(7X,10(F9.6,6X))',
     &               Work(LHRot),lroots,lroots)
         if(DoGradMSPD) then
           if(DoNACMSPD) then
              Call MSPDFTNAC_Misc(LHRot)
           ELSE
              Call MSPDFTGrad_Misc(LHRot)
           end if
           Call GetMem('F1MS' ,'Free','Real',iF1MS , nTot1*nRoots)
           Call GetMem('F2MS' ,'Free','Real',iF2MS ,nACPR2*nRoots)
           Call GetMem('FxyMS','Free','Real',iFxyMS, nTot4*nRoots)
           Call GetMem('P2MO' ,'Free','Real',iP2MOt,nACPR2*nRoots)
           Call GetMem('FocMS','Free','Real',iFocMS, nTot1*nRoots)
           Call GetMem('DIDA' ,'Free','Real',iDIDA ,nTot1*(nRoots+1))
           Call GetMem('D1AOMS' ,'Free','Real',D1AOMS,nTot1*nRoots)
           if (ispin.ne.1)
     &     Call GetMem('D1SAOMS' ,'Free','Real',D1SAOMS,nTot1*nRoots)
          end if

          refbas=.false.
          call f_inquire('ROT_VEC',RefBas)
          Call GetMem('XScratch','FREE','Real',LXScratch,NXScratch)
*print MS-PDFT final states in basis of reference states
*re-use RotStat, XScratch and LRState
          if(RefBas) then
            NXScratch=NHRot
            Call GetMem('XScratch','ALLO','Real',LXScratch,NXScratch)
            Call FZero(Work(LXScratch),NXScratch)
            Call FZero(Work(LRState)  ,NXScratch)
            CALL ReadMat2('ROT_VEC',MatInfo,WORK(LRState),
     &                     lRoots,lRoots,7,18,'T')
            CALL DGEMM_('n','n',lRoots,lRoots,lRoots,1.0d0,
     &          Work(LRState),lRoots,Work(LHRot),lRoots,0.0d0,
     &          Work(LXScratch),lRoots)
            write(6,'(7X,A)')'Reference-state Basis'
            write(6,mspdftfmt)((VecStat(JRoot)),JRoot=1,lroots)
            Call RecPrt(' ','(7X,10(F9.6,6X))',
     &                Work(LXScratch),lroots,lroots)
            CALL PrintMat2('FIN_VEC',MatInfo,WORK(LXScratch),
     &                      lRoots,lRoots,7,18,'T')
            Call GetMem('XScratch','FREE','Real',LXScratch,NXScratch)
          end if
*        Gradient part
          if(DoGradMSPD) then
            Call Put_iScalar('Number of roots',nroots)
            Call Put_cArray('Relax Method','MSPDFT  ',8)
            Call Put_cArray('MCLR Root','****************',16)
            Call Put_iScalar('Relax CASSCF root',irlxroot)
          end if
          Write(6,'(6X,80a)') ('*',i=1,80)
          CALL GETMEM('HRot','FREE','REAL',LHRot,NHRot)
          CALL GETMEM('RotStat','FREE','REAL',LRState,NRState)
          CALL mma_deallocate(VecStat)
        End If
        CALL GETMEM('CASDFT_Fock','FREE','REAL',LFOCK,NACPAR)
      END IF


      ICICH=0
*****************************************************************************************
***************************           Closing up MC-PDFT      ***************************
*****************************************************************************************

************************************************************************
*^follow closing up MC-PDFT
*
* release SEWARD
*
      Call ClsSew()
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

*  Release  some memory allocations
      Call GetMem('FOCC','FREE','REAL',ipFocc,NTOT1)
      Call GetMem('FI','Free','Real',LFI,NTOT1)
      Call GetMem('FA','Free','Real',LFA,NTOT1)
      Call GetMem('D1I','Free','Real',LD1I,NTOT2)
      Call GetMem('D1A','Free','Real',LD1A,NTOT2)
      Call GetMem('D1tot','Free','Real',lD1tot,NTOT1)
      Call GetMem('LCMO','Free','Real',LCMO,NTOT2)
      Call GetMem('REF_E','Free','REAL',iRef_E,lroots)
      Call GetMem('OCCN','Free','Real',LOCCN,NTOT)

      If ( NAC.GT.0 ) then
        Call GetMem('DMAT','free','Real',LDMAT,NACPAR)
        Call GetMem('DSPN','free','Real',LDSPN,NACPAR)
        Call GetMem('PMAT','free','Real',LPMAT,NACPR2)
        Call GetMem('P2AS','free','Real',LPA,NACPR2)
        Call GetMem('TUVX','free','Real',LTUVX,NACPR2)
      End if

      If (iClean.eq.1) Call Free_iWork(ipCleanMask)
*
*
       if (doGSOR) then
          CALL GETMEM('CIVEC','FREE','REAL',LW4,NCONF)
          CALL GETMEM('Dtmp ','FREE','REAL',LW6,NACPAR)
          CALL GETMEM('DStmp','FREE','REAL',LW7,NACPAR)
          CALL GETMEM('Ptmp ','FREE','REAL',LW8,NACPR2)
          CALL GETMEM('PAtmp','FREE','REAL',LW9,NACPR2)
          CALL GETMEM('Pscr','FREE','REAL',LW10,NACPR2)
          Call GetMem('CIVtmp','FREE','Real',LW11,nConf)
          Call Lucia_Util('CLOSE',iDummy,iDummy,Dummy)
          Call MKGUGA_FREE_m
       end if

*
      Call StatusLine('MCPDFT:','Finished.')
      If (IPRLEV.GE.2) Write(LF,*)


      Call Timing(Swatch,Swatch,Ebel_3,Swatch)
      IF (IPRLEV.GE.3) THEN
       Call PrtTim_m
       Call FastIO('STATUS')
      END IF
      Call ClsFls_RASSCF_m()

      If (Do_OFemb) Then
         Call GetEnvF('EMIL_InLoop',EMILOOP)
         If (EMILOOP.eq.' ') EMILOOP='0'
         If (EMILOOP(1:1).ne.'0') Then
            If (iReturn.ne._RC_ALL_IS_WELL_) Then
               Call WarningMessage(1,'RASSCF: non-zero return code.')
            EndIf
            iReturn=_RC_CONTINUE_LOOP_
            Call Check_FThaw(iReturn)
         EndIf
      EndIf

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
      End

