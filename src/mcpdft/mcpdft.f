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
      use Fock_util_global, only: DoCholesky
      use mcpdft_input, only: mcpdft_options, parse_input
      use write_pdft_job, only: writejob
      use mspdft, only: mspdft_finalize, heff, mspdft_init
      use printlevel, only: terse, debug, insane, usual
      use mcpdft_output, only: lf, iPrLoc
      use mspdft_util, only: replace_diag
      use stdalloc, only: mma_allocate, mma_deallocate
      use wadr, only: DMAT, PMAT, PA, FockOcc, TUVX, FI, FA, DSPN,
     &                D1I, D1A, OccN, CMO
      use rasscf_global, only: ECAS, ExFac, IPR, ITER,
     &                         lRoots, lSquare, NAC, NACPAR, NACPR2,
     &                         nFint, iRoot, Weight,
     &                         Ener

      Implicit None

      Integer iReturn
#include "rasdim.fh"
#include "warnings.h"
#include "general.fh"
#include "timers.fh"

      Logical IfOpened
      Logical Found

      integer IAD19
      integer IADR19(1:15)
      integer NMAYBE,KROOT
      real*8 EAV
!
      Real*8, allocatable :: TmpDMat(:), Ref_E(:), EList(:,:),
     &                       PUVX(:)
      Logical DSCF
      Real*8 AEMAX, dum1, dum2, dum3, E
      Integer I, iJOB, iPrLev, iRC, IT, jDisk

      Call StatusLine('MCPDFT:',' Just started.')
      IRETURN=_RC_ALL_IS_WELL_

! Local print level in this routine:
      IPRLEV=IPRLOC(1)

! Default option switches and values, and initial data.
      EAV = 0.0d0
      call mcpdft_init()

      call parse_input()


! Local print level in this routine:
      IPRLEV=IPRLOC(1)

      Call open_files_mcpdft(DSCF)

! Some preliminary input data:
      Call Rd1Int()
      If ( .not. DSCF ) then
        Call Rd2Int_RASSCF()
      end if

! Process the input:
      Call Proc_InpX(DSCF,iRc)

! If something goes wrong in proc_inp:
      If (iRc.ne._RC_ALL_IS_WELL_) Then
        If (IPRLEV.ge.TERSE) Then
          Call WarningMessage(2,'Input processing failed.')
          write(lf,*)' MC-PDFT Error: Proc_Inp failed unexpectedly.'
        End If
        IRETURN=iRc
        call close_files()
        Return
      End If


! Local print level may have changed:
      IPRLEV=IPRLOC(1)

      if(mcpdft_options%mspdft) then
        call mspdft_init()
      endif

      Call InpPri_m()

! Allocate various matrices
      Call mma_allocate(FI,NTOT1,Label='FI')
      Call mma_allocate(FA,NTOT1,Label='FA')
      Call mma_allocate(D1I,NTOT2,Label='D1I')
      Call mma_allocate(D1A,NTOT2,Label='D1A')
      Call mma_allocate(OCCN,NTOT,Label='OccN')
      Call mma_allocate(CMO,NTOT2,Label='CMO')
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

      Call mma_allocate(TmpDMAT,NACPAR,Label='TmpDMat')
      call dcopy_(NACPAR,DMAT,1,TmpDMAT,1)
      If (NASH(1).ne.NAC) then
        Call DBLOCK(TmpDMAT)
      end if
      Call Get_D1A_RASSCF(CMO,TmpDMAT,D1A)
      Call mma_deallocate(TmpDMAT)

!AMS start-
! - Read in the CASSCF Energy from JOBIPH file.  These values are not
! used in calculations, but are merely reprinted as the reference energy
! for each calculated MC-PDFT energy.
      iJOB=0
      Call mma_allocate(Ref_E,lroots,Label='Ref_E')
      Ref_E(:)=0.0D0
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
      Call mma_allocate(EList,MXROOT,MXITER,Label='EList')
      Call DDaFile(JOBOLD,2,EList,MXROOT*MXITER,jdisk)
      NMAYBE=0
      DO IT=1,MXITER
        AEMAX=0.0D0
        DO I=1,MXROOT
          E=EList(I,IT)
          AEMAX=MAX(AEMAX,ABS(E))
        END DO
        IF(ABS(AEMAX).LE.1.0D-12) GOTO 11
        NMAYBE=IT
      END DO
  11  CONTINUE

      IF(mcpdft_options%mspdft) Then
        do KROOT=1,lROOTS
           ener(iroot(kroot),1)=heff(kroot,kroot)
           EAV = EAV + ENER(IROOT(KROOT),ITER) * WEIGHT(KROOT)
           Ref_E(KROOT) = ENER(IROOT(KROOT),1)
        end do
      Else
        do KROOT=1,lROOTS
          ENER(IROOT(KROOT),1)=EList(KRoot,NMAYBE)
           EAV = EAV + ENER(IROOT(KROOT),ITER) * WEIGHT(KROOT)
           Ref_E(KROOT) = ENER(IROOT(KROOT),1)
        end do
      End IF

      call mma_deallocate(elist)
      If(JOBOLD.gt.0.and.JOBOLD.ne.JOBIPH) Then
        Call DaClos(JOBOLD)
        JOBOLD=-1
      End if

      IF(mcpdft_options%nac) Then
        IF(IPRLEV.ge.USUAL) THEN
        write(6,'(6X,A)') repeat('=',80)
        write(6,*)
        write(6,'(6X,A,I3,I3)')'keyword NAC is used for states:',
     &         mcpdft_options%nac_states(1),
     &         mcpdft_options%nac_states(2)
        write(6,*)
        write(6,'(6X,A)') repeat('=',80)
        END IF
      ELSE
        mcpdft_options%nac_states(1) = mcpdft_options%rlxroot
        mcpdft_options%nac_states(2) = 0
      End IF
      call Put_lScalar('isCMSNAC        ', mcpdft_options%nac)
      call Put_iArray('cmsNACstates    ', mcpdft_options%nac_states, 2)

      IF(mcpdft_options%meci .and. iprlev .ge. usual) Then
        write(lf,'(6X,A)') repeat('=',80)
        write(lf,*)
        write(lf,'(6X,A,I3,I3)')'keyword MECI is used for states:'
        write(lf,*)
        write(lf,'(6X,A)') repeat('=',80)
      End IF
      call Put_lScalar('isMECIMSPD      ', mcpdft_options%meci)


*
* Transform two-electron integrals and compute at the same time
* the Fock matrices FI and FA
*
      Call Timing(dum1,dum2,Fortis_1,dum3)
      Call mma_allocate(PUVX,NFINT,Label='PUVX')
      PUVX(:)=0.0D0
      Call Get_D1I_RASSCF(CMO,D1I)

      IPR=0
      IF(IPRLOC(2).EQ.debug) IPR=5
      IF(IPRLOC(2).EQ.insane) IPR=10

      CALL TRACTL2(CMO,PUVX,TUVX,D1I,FI,D1A,FA,IPR,lSquare,ExFac)

      Call Put_dArray('Last orbitals',CMO,ntot2)

      if(mcpdft_options%grad .and. mcpdft_options%mspdft) then
        CALL Put_dArray('TwoEIntegral    ',PUVX,nFINT)
      end if
      Call mma_deallocate(PUVX)

      Call Timing(dum1,dum2,Fortis_2,dum3)
      Fortis_2 = Fortis_2 - Fortis_1
      Fortis_3 = Fortis_3 + Fortis_2

      ! This is where MC-PDFT actually computes the PDFT energy for
      ! each state
      ! only after 500 lines of nothing above...
      Call MSCtl(CMO,FI,Ref_E)

      ! I guess Ref_E now holds the MC-PDFT energy for each state??

      If(mcpdft_options%wjob .and.(.not.mcpdft_options%mspdft)) then
        Call writejob(iadr19)
      end if

        If (mcpdft_options%mspdft) Then
          call replace_diag(heff, ref_e, lroots)
          call mspdft_finalize(lroots, iadr19)
        End If

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
      Call mma_deallocate(REF_E)

      Call mma_deallocate(DMAT)
      Call mma_deallocate(DSPN)
      Call mma_deallocate(PMAT)
      Call mma_deallocate(PA)
      Call mma_deallocate(TUVX)

      Call StatusLine('MCPDFT:','Finished.')
      If (IPRLEV.GE.2) Write(LF,*)

      Call Timing(dum1,dum2,Ebel_3,dum3)
      IF (IPRLEV.GE.3) THEN
       Call PrtTim()
       Call FastIO('STATUS')
      END IF

      Call close_files()

      Contains
      subroutine close_files()
      Integer I
      call close_files_mcpdft()
      DO I=10,99
        INQUIRE(UNIT=I,OPENED=IfOpened)
        IF (IfOpened.and.I.ne.19) CLOSE (I)
      END DO
      End subroutine close_files

      End SUBROUTINE MCPDFT

