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
      use definitions,only:iwp,wp,u6
      use constants,only:zero
      use Fock_util_global, only: DoCholesky
      use mcpdft_input, only: mcpdft_options, parse_input
      use write_pdft_job, only: writejob
      use mspdft, only: mspdft_finalize, heff, mspdft_init
      use printlevel, only: terse, debug, insane
      use mcpdft_output, only: iPrLoc
      use mspdft_util, only: replace_diag
      use stdalloc, only: mma_allocate, mma_deallocate
      use wadr, only: FockOcc
      use general_data,only:ntot1,ntot2
      use rasscf_global, only: ExFac, IPR,
     &                         lRoots, lSquare, NACPR2,
     &                         nFint

      Implicit None

      integer(kind=iwp), intent(out) :: iReturn
#include "warnings.h"
#include "timers.fh"

      logical(kind=iwp) :: dscf
      real(kind=wp), allocatable :: ref_e(:), PUVX(:), D1I(:), D1A(:)
      real(kind=wp), allocatable :: cmo(:), FI(:), FA(:), tuvx(:)
      Real(kind=wp) :: dum1, dum2, dum3
      integer(kind=iwp) :: iPrLev, iRC, state

      Call StatusLine('MCPDFT: ','Just started.')
      IRETURN=_RC_ALL_IS_WELL_

! Local print level in this routine:
      IPRLEV=IPRLOC(1)

! Default option switches and values, and initial data.
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
          write(u6,*)' MC-PDFT Error: Proc_Inp failed unexpectedly.'
        End If
        IRETURN=iRc
        Return
      End If


! Local print level may have changed:
      IPRLEV=IPRLOC(1)

      if(mcpdft_options%mspdft) then
        call mspdft_init()
      endif

      Call InpPri_m()

! Get start orbitals
      Call mma_allocate(CMO,NTOT2,Label='CMO')
      Call ReadVc_m(CMO)

! Allocate core space for dynamic storage of data
! Really just determine size of 2-e integrals
      CALL ALLOC()

      Call Timing(dum1,dum2,Ebel_1,dum3)

! This needs to be put somehwere else...
      Call mma_allocate(FockOcc,nTot1,Label='FockOcc')


! - Read in the CASSCF Energy from JOBIPH file.
      Call mma_allocate(ref_e,lroots,Label='Ref_E')
      ref_e(:)=0.0D0
      call ref_energy(ref_e,lroots)

! Transform two-electron integrals and compute at the same time
! the Fock matrices FI and FA
      Call Timing(dum1,dum2,Fortis_1,dum3)

      IPR=0
      IF(IPRLOC(2).EQ.debug) IPR=5
      IF(IPRLOC(2).EQ.insane) IPR=10



      ! so this does 2 things, first it puts the integrals
      ! into the file LUINTM, the second is that we write
      ! the integrals to the runfile anyways for mspdft grad

      if(mcpdft_options%grad .and. mcpdft_options%mspdft) then
        Call mma_allocate(D1I,NTOT2,Label='D1I')
        Call mma_allocate(D1A,NTOT2,Label='D1A')
        Call mma_allocate(FI,NTOT1,Label='FI')
        Call mma_allocate(FA,NTOT1,Label='FA')
        call mma_allocate(puvx,nfint,label='PUVX')
        Call mma_allocate(TUVX,NACPR2,Label='TUVX')
        D1I(:) = zero
        D1A(:) = zero
      CALL TRACTL2(CMO,PUVX,TUVX,D1I,
     &             FI,D1A,FA,IPR,lSquare,ExFac)
        CALL Put_dArray('TwoEIntegral    ',PUVX,nFINT)
        call mma_deallocate(PUVX)
        Call mma_deallocate(FI)
        Call mma_deallocate(FA)
        Call mma_deallocate(TUVX)
        Call mma_deallocate(D1I)
        Call mma_deallocate(D1A)
      end if

      Call Timing(dum1,dum2,Fortis_2,dum3)
      Fortis_2 = Fortis_2 - Fortis_1
      Fortis_3 = Fortis_3 + Fortis_2

      ! This is where MC-PDFT actually computes the PDFT energy for
      ! each state
      ! only after 500 lines of nothing above...
      Call MSCtl(CMO,ref_e)
      Call mma_deallocate(CMO)

      ! I guess ref_e now holds the MC-PDFT energy for each state??

      If(mcpdft_options%wjob .and.(.not.mcpdft_options%mspdft)) then
        Call writejob(ref_e,lroots)
      end if

      If (mcpdft_options%mspdft) Then
        Call Put_cArray('Relax Method','MSPDFT  ',8)
        call replace_diag(heff, ref_e, lroots)
        call mspdft_finalize(lroots)
      else
        Call Put_cArray('Relax Method','MCPDFT  ',8)
        if(iprlev >= terse) then
        do state=1,lroots
          call PrintResult(u6,'(6X,A,I3,A,F16.8)',
     &          'MCPDFT root number',state,
     &          ' Total energy:',[ref_e(state)],1)
        enddo
      endif
      End If

!***********************************************************************
!**************************           Closing up MC-PDFT      **********
!***********************************************************************

! release SEWARD
      Call ClsSew()

! Finalize Cholesky information if initialized
      if (DoCholesky)then
         Call Cho_X_Final(irc)
         if (irc.ne.0) then
           Write(u6,*)'MC-PDFT: Cho_X_Final fails with return code ',irc
           Write(u6,*)' Try to recover. Calculation continues.'
         endif
      endif

! Release  some memory allocations
      Call mma_deallocate(FockOcc)
      call mma_deallocate(ref_e)


      Call StatusLine('MCPDFT: ','Finished.')
      If (IPRLEV.GE.2) Write(u6,*)

      Call Timing(dum1,dum2,Ebel_3,dum3)
      IF (IPRLEV.GE.3) THEN
       Call PrtTim()
       Call FastIO('STATUS')
      END IF

      call close_files_mcpdft()

      End SUBROUTINE MCPDFT

