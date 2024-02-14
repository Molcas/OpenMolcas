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

      use csfbas, only: CONF
      use glbbas, only: CFTP
      use Fock_util_global, only: DoCholesky
      use write_pdft_job, only: iwjob, writejob
      use sxci, only: idxsx
      use mspdft_grad, only: dogradmspd
      use mspdft, only: mspdftmethod, do_rotate, iF1MS,
     &                  iF2MS, iFxyMS, iFocMS, iDIDA, IP2MOt, D1AOMS,
     &                  D1SAOMS, doNACMSPD, cmsNACstates, doMECIMSPD,
     &                  mspdft_finalize
      use mcpdft_output, only: terse, debug, insane, lf, iPrLoc
      use mspdft_util, only: replace_diag
      use rctfld_module
      use rasscf_lucia, only: PAtmp, Pscr, CIVEC, Ptmp, DStmp, Dtmp
      use stdalloc, only: mma_allocate, mma_deallocate
      use lucia_interface, only: lucia_util
      use wadr, only: DMAT, PMAT, PA, FockOcc, TUVX, FI, FA, DSPN,
     &                D1I, D1A, OccN, CMO
      use gugx, only: IFCAS

      Implicit Real*8 (A-H,O-Z)

#include "WrkSpc.fh"
#include "rasdim.fh"
#include "warnings.h"
#include "input_ras_mcpdft.fh"
#include "rasscf.fh"
#include "general.fh"
#include "gas.fh"
#include "timers.fh"
#include "lucia_ini.fh"
#include "pamint.fh"
#include "ciinfo.fh"
      Integer LHrot,NHrot             ! storing info in H0_Rotate.txt

      CHARACTER(Len=18)::MatInfo
      INTEGER LUMS,IsFreeUnit
      External IsFreeUnit

      Logical DSCF
      Character*80 Line
      Logical IfOpened
      Logical Found

      External RasScf_Init_m
      External Scan_Inp_m
      integer iRef_E,IAD19
      integer IADR19(1:15)
      integer NMAYBE,KROOT
      real*8 EAV
!
      real*8, allocatable :: PLWO(:), CIV(:)
      integer ivkcnf
* Set status line for monitor:
      Call StatusLine('MCPDFT:',' Just started.')
* Set the return code(s)
      IRETURN=_RC_ALL_IS_WELL_

* Local print level in this routine:
      IPRLEV=IPRLOC(1)

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
       write(lf,*)' MCPDFT Error: Could not make a clean copy of'
       write(lf,*)' the input file. This is an unexpected bug.'
       IRETURN=_RC_INTERNAL_ERROR_
       GOTO 9990
      End If


! Scan the input file for keywords:
      Call Scan_Inp_m(iRc)
! If something wrong with input file:
      If (iRc.ne._RC_ALL_IS_WELL_) Then
       If (IPRLOC(1).GE.TERSE) Then
        Call WarningMessage(2,'Scanning input file failed.')
! Calling again, now with iRc indicating an error,
! will echo the keywords:
        Call Scan_Inp_m(iRc)
       End If
       IRETURN=_RC_INPUT_ERROR_
       GOTO 9990
      End If

! Local print level in this routine:
      IPRLEV=IPRLOC(1)

* Open files
      Call OpnFls_RASSCF_m(DSCF)

* Some preliminary input data:
      Call Rd1Int()
      If ( .not.DSCF ) Call Rd2Int_RASSCF()

* Process the input:
      Call Proc_InpX(DSCF,iRc)
* If something goes wrong in proc_inp:
      If (iRc.ne._RC_ALL_IS_WELL_) Then
        If (IPRLEV.ge.TERSE) Then
          Call WarningMessage(2,'Input processing failed.')
          write(lf,*)' MC-PDFT Error: Proc_Inp failed unexpectedly.'
          write(lf,*)' Here is a printing of the input file that'
          write(lf,*)' was processed:'
          Rewind(LUInput)
  15      Continue
          Read(LuInput,'(A80)',End=16,Err=16) Line
          write(lf,*) Line
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
      Call mma_allocate(FI,NTOT1,Label='FI')
      Call mma_allocate(FA,NTOT1,Label='FA')
      Call mma_allocate(D1I,NTOT2,Label='D1I')
      Call mma_allocate(D1A,NTOT2,Label='D1A')
      Call mma_allocate(OCCN,NTOT,Label='OccN')
      Call mma_allocate(CMO,NTOT2,Label='CMO')
      allocate(PLWO(1:NACPAR))
      PLWO(:) = 0
!
*
      Call mma_allocate(TUVX,NACPR2,Label='TUVX')
      TUVX(:)=0.0D0
      Call mma_allocate(DSPN,NACPAR,Label='DSPN')
      Call mma_allocate(DMAT,NACPAR,Label='DMAT')
      DMAT(:)=0.0D0
      Call mma_allocate(PMAT,NACPR2,Label='PMAT')
      Call mma_allocate(PA,NACPR2,Label='PA')
*
* Get start orbitals

* Initialize OCCN array, to prevent false alarms later from
* automated detection of using uninitialized variables:
      OccN(:)=0.0D0

* PAM03: Note that removal of linear dependence may change the nr
* of secondary/deleted orbitals, affecting some of the global
* variables: NSSH(),NDEL(),NORB(),NTOT3, etc etc
      Call ReadVc_m(CMO,OCCN,DMAT,DSPN,PMAT,PA)
* Only now are such variables finally known.
      If (IPRLOC(1).GE.DEBUG) Then
        CALL TRIPRT('Averaged one-body density matrix, D, in RASSCF',
     &              ' ',DMAT,NAC)
        CALL TRIPRT('Averaged one-body spin density matrix DS, RASSCF',
     &              ' ',DSPN,NAC)
        CALL TRIPRT('Averaged two-body density matrix, P',
     &              ' ',PMAT,NACPAR)
        CALL TRIPRT('Averaged antisym 2-body density matrix PA RASSCF',
     &              ' ',PA,NACPAR)
      END IF
*
* Allocate core space for dynamic storage of data
*
      CALL ALLOC()

      Call Timing(dum1,dum2,Ebel_1,dum3)

      ECAS   = 0.0d0
      Call mma_allocate(FockOcc,nTot1,Label='FockOcc')

      ! I guess we spoof for the 2-electron part? Im not sure..
        KSDFT_TEMP=KSDFT
        KSDFT='SCF'
        ExFac=1.0D0

      Call GetMem('TmpDMAT','Allo','Real',ipTmpDMAT,NACPAR)
      call dcopy_(NACPAR,DMAT,1,Work(ipTmpDMAT),1)
      If (NASH(1).ne.NAC) then
        Call DBLOCK(Work(ipTmpDMAT))
      end if
      Call Get_D1A_RASSCF(CMO,Work(ipTmpDMAT),D1A)
      Call GetMem('TmpDMAT','Free','Real',ipTmpDMAT,NACPAR)

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

      IF(iMSPDFT==1) Then
       call f_inquire('ROT_HAM',Do_Rotate)
       If(.not.Do_Rotate) Then
        write(lf,'(6X,A,A)')'keyword "MSPD" is used but ',
     &  'the file of rotated Hamiltonian is not found.'
        write(lf,'(6X,2a)')'Performing regular (state-',
     &   'specific) MC-PDFT calculation'
       End If
      End IF
      IF(Do_Rotate) Then
        write(lf,'(6X,80A)') ('=',i=1,80)
        write(lf,*)
        write(lf,'(6X,A,A)')'keyword "MSPD" is used and ',
     &  'file recording rotated hamiltonian is found. '
        write(lf,*)
        write(lf,'(6X,A,A)')
     &  'Switching calculation to Multi-State Pair-Density ',
     &  'Functional Theory (MS-PDFT) '
        write(lf,'(6X,A)')'calculation.'
        write(lf,*)
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
         write(lf,'(6X,A,A)')'The MS-PDFT calculation is ',
     & 'based on a user-supplied rotation matrix.'
        ELSE
         write(lf,'(6X,A,A,A)')'The MS-PDFT method is ',
     &   trim(adjustl(MatInfo)),'.'
        If(trim(adjustl(MatInfo)).eq.'XMS-PDFT') MSPDFTMethod='XMS-PDFT'
        If(trim(adjustl(MatInfo)).eq.'CMS-PDFT') MSPDFTMethod='CMS-PDFT'
        If(trim(adjustl(MatInfo)).eq.'VMS-PDFT') MSPDFTMethod='VMS-PDFT'
        If(trim(adjustl(MatInfo)).eq.'FMS-PDFT') MSPDFTMethod='FMS-PDFT'
        ENDIF
        write(lf,*)
        write(lf,'(6X,80A)') ('=',i=1,80)
        write(lf,*)
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

      ! now we reset..
        KSDFT=KSDFT_TEMP
        ExFac=0.0d0
*
* Transform two-electron integrals and compute at the same time
* the Fock matrices FI and FA
*
      Call Timing(dum1,dum2,Fortis_1,dum3)
      Call GetMem('PUVX','Allo','Real',LPUVX,NFINT)
      Call FZero(Work(LPUVX),NFINT)
      Call Get_D1I_RASSCF(CMO,D1I)

      IPR=0
      IF(IPRLOC(2).EQ.debug) IPR=5
      IF(IPRLOC(2).EQ.insane) IPR=10

      CALL TRACTL2(CMO,WORK(LPUVX),TUVX,D1I,
     &             FI,D1A,FA,IPR,lSquare,ExFac)

      Call Put_dArray('Last orbitals',CMO,ntot2)

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
        CALL mma_allocate(CIVEC,NCONF,Label='CIVEC')
        CALL mma_allocate(Dtmp,NAC**2,Label='Dtmp')
        CALL mma_allocate(DStmp,NAC**2,Label='DStmp')
        CALL mma_allocate(Ptmp,NACPR2,Label='Ptmp')
        CALL mma_allocate(PAtmp,NACPR2,Label='PAtmp')
        CALL mma_allocate(Pscr,NACPR2,Label='Pscr')

        Dtmp(:)=0.0D0
        DStmp(:)=0.0D0
        Ptmp(:)=0.0D0
        CIVEC(:)=0.0D0
        iDisk = IADR19(4)
        jDisk = IADR19(3)

        Call mma_allocate(CIV,nConf,Label='CIV')
        DO jRoot=1,lroots
          do i=1,nconf
            read(LUCT,*) CIVEC(i)
          end do
          Call DDafile(JOBOLD,1,CIVEC,nConf,iDisk)
          call getmem('kcnf','allo','inte',ivkcnf,nactel)
          Call Reord2(NAC,NACTEL,STSYM,1,
     &                CONF,CFTP,CIVEC,CIV,iWork(ivkcnf))
          Call dcopy_(nconf,CIV,1,CIVEC,1)
          call getmem('kcnf','free','inte',ivkcnf,nactel)
          CALL Lucia_Util('Densi',
     &                    CI_Vector=CIVEC(:))
          If (IFCAS > 2) Then
            Call CISX(IDXSX,Dtmp,DStmp,Ptmp,PAtmp,PScr)
          End If

          Call DDafile(JOBOLD,1,Dtmp,NACPAR,jDisk)
          Call DDafile(JOBOLD,1,DStmp,NACPAR,jDisk)
          Call DDafile(JOBOLD,1,Ptmp,NACPR2,jDisk)
          Call DDafile(JOBOLD,1,PAtmp,NACPR2,jDisk)
        end do
        Close(LUCT)
        Call fCopy('JOBIPH','JOBGS',ierr)

       end if!DoGSOR

      if(dogradmspd) then
        CALL Put_dArray('TwoEIntegral    ',Work(LPUVX),nFINT)
      end if
      Call GetMem('PUVX','Free','Real',LPUVX,NFINT)

      Call Timing(dum1,dum2,Fortis_2,dum3)
      Fortis_2 = Fortis_2 - Fortis_1
      Fortis_3 = Fortis_3 + Fortis_2

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

      ! This is where MC-PDFT actually computes the PDFT energy for
      ! each state
      ! only after 500 lines of nothing above...
      Call MSCtl(CMO,FI,FA,Work(iRef_E))

      ! I guess iRef_E now holds the MC-PDFT energy for each state??

        If(IWJOB==1.and.(.not.Do_Rotate)) Call writejob(iadr19)

        If (Do_Rotate) Then
          call replace_diag(work(lhrot), work(iref_e), lroots)
          call mspdft_finalize(work(lhrot), lroots, irlxroot, iadr19)
        End If

      ! Free up some space
      CALL GETMEM('CASDFT_Fock','FREE','REAL',LFOCK,NACPAR)

      if (do_rotate) then
        CALL GETMEM('HRot','FREE','REAL',LHRot,NHRot)
        if(DoGradMSPD) then
          Call GetMem('F1MS' ,'Free','Real',iF1MS , nTot1*nRoots)
          Call GetMem('F2MS' ,'Free','Real',iF2MS ,nACPR2*nRoots)
          Call GetMem('FxyMS','Free','Real',iFxyMS, nTot4*nRoots)
          Call GetMem('P2MO' ,'Free','Real',iP2MOt,nACPR2*nRoots)
          Call GetMem('FocMS','Free','Real',iFocMS, nTot1*nRoots)
          Call GetMem('DIDA' ,'Free','Real',iDIDA ,nTot1*(nRoots+1))
          Call GetMem('D1AOMS' ,'Free','Real',D1AOMS,nTot1*nRoots)
          if (ispin.ne.1)
     &       Call GetMem('D1SAOMS' ,'Free','Real',D1SAOMS,nTot1*nRoots)
        end if
      end if

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
           Write(LF,*)'MC-PDFT: Cho_X_Final fails with return code ',irc
           Write(LF,*)' Try to recover. Calculation continues.'
         endif
      endif

*  Release  some memory allocations
      Call mma_deallocate(FockOcc)
      Call mma_deallocate(FI)
      Call mma_deallocate(FA)
      Call mma_deallocate(D1I)
      Call mma_deallocate(D1A)
      Call mma_deallocate(OccN)
      Call mma_deallocate(CMO)
      Call GetMem('REF_E','Free','REAL',iRef_E,lroots)

      Call mma_deallocate(DMAT)
      Call mma_deallocate(DSPN)
      Call mma_deallocate(PMAT)
      Call mma_deallocate(PA)
      Call mma_deallocate(TUVX)

*
*
       if (doGSOR) then
          CALL mma_deallocate(CIVEC)
          CALL mma_deallocate(Dtmp)
          CALL mma_deallocate(DStmp)
          CALL mma_deallocate(Ptmp)
          CALL mma_deallocate(PAtmp)
          CALL mma_deallocate(Pscr)
          CALL mma_deallocate(CIV)
          Call Lucia_Util('CLOSE')
          Call MKGUGA_FREE()
       end if

*
      Call StatusLine('MCPDFT:','Finished.')
      If (IPRLEV.GE.2) Write(LF,*)

      Call Timing(dum1,dum2,Ebel_3,dum3)
      IF (IPRLEV.GE.3) THEN
       Call PrtTim()
       Call FastIO('STATUS')
      END IF
      Call ClsFls_RASSCF_m()

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

      DO I=10,99
        INQUIRE(UNIT=I,OPENED=IfOpened)
        IF (IfOpened.and.I.ne.19) CLOSE (I)
      END DO
      Close(LUInput)

      End

