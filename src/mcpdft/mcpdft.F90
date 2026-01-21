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
!               2016, Andrew M. Sand                                   *
!               2024, Matthew R. Hennefarth                            *
!***********************************************************************

subroutine mcpdft(ireturn)

use Fock_util_global, only: DoCholesky
use timers, only: TimeInput, TimeOutput, TimeTotal, TimeTrans
use write_pdft_job, only: writejob
use PrintLevel, only: DEBUG, INSANE, TERSE
use mcpdft_output, only: iPrLoc
use mcpdft_input, only: mcpdft_options, parse_input
use mspdft_util, only: replace_diag
use mspdft, only: heff, mspdft_init, mspdft_finalize
use wadr, only: FockOcc
use general_data, only: nelec3, nhole1, ntot1, ntot2
use rasscf_global, only: ExFac, IPR, lRoots, lSquare, NACPR2, nFint
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: iReturn
#include "warnings.h"
integer(kind=iwp) :: iPrLev, iRC, state
logical(kind=iwp) :: dscf
real(kind=wp) :: dum1, dum2, dum3, time1(2), time2(2)
character(len=8) :: Method
real(kind=wp), allocatable :: cmo(:), D1A(:), D1I(:), e_mcscf(:), e_states(:), FA(:), FI(:), PUVX(:), tuvx(:)

call StatusLine('MCPDFT:',' Just started.')
IRETURN = _RC_ALL_IS_WELL_

! Local print level in this routine:
IPRLEV = IPRLOC(1)

! Default option switches and values, and initial data.
call mcpdft_init()

call parse_input()

! Local print level in this routine:
IPRLEV = IPRLOC(1)

call open_files_mcpdft(DSCF)

! Some preliminary input data:
call Rd1Int()
if (.not. DSCF) call Rd2Int_RASSCF()

! Process the input:
call Proc_InpX(DSCF,iRc)

! Delayed checks
! (This should be in mcpdft_init, but nhole1, nelec3 are set in proc_inpx)
if (mcpdft_options%grad) then
  if ((nhole1 > 0) .or. (nelec3 > 1)) then
    call WarningMessage(2,'MC-PDFT gradients with a RASSCF wave function not supported')
    write(u6,*) ' ************* ERROR **************'
    write(u6,*) ' MC-PDFT gradients are not         '
    write(u6,*) ' implemented with RASSCF           '
    write(u6,*) ' **********************************'
    call Quit_OnUserError()
  end if
end if

! If something goes wrong in proc_inp:
if (iRc /= _RC_ALL_IS_WELL_) then
  if (IPRLEV >= TERSE) then
    call WarningMessage(2,'Input processing failed.')
    write(u6,*) ' MC-PDFT Error: Proc_Inp failed unexpectedly.'
  end if
  IRETURN = iRc
  return
end if

! Local print level may have changed:
IPRLEV = IPRLOC(1)

if (mcpdft_options%mspdft) call mspdft_init()

call InpPri_m()

! Get start orbitals
call mma_allocate(CMO,NTOT2,Label='CMO')
call ReadVc_m(CMO)

! Allocate core space for dynamic storage of data
! Really just determine size of 2-e integrals
call ALLOC()

call Timing(dum1,dum2,time1(1),dum3)
TimeInput = time1(1)

! This needs to be put somehwere else...
call mma_allocate(FockOcc,nTot1,Label='FockOcc')

! - Read in the CASSCF Energy from JOBIPH file.
call mma_allocate(e_mcscf,lroots,Label='e_mcscf')
e_mcscf(:) = Zero
call ref_energy(e_mcscf,lroots)

! Transform two-electron integrals and compute at the same time
! the Fock matrices FI and FA
call Timing(dum1,dum2,time2(1),dum3)

! so this does 2 things, first it puts the integrals
! into the file LUINTM, the second is that we write
! the integrals to the runfile anyways for mspdft grad

if (mcpdft_options%grad .and. mcpdft_options%mspdft) then
  call mma_allocate(D1I,NTOT2,Label='D1I')
  call mma_allocate(D1A,NTOT2,Label='D1A')
  call mma_allocate(FI,NTOT1,Label='FI')
  call mma_allocate(FA,NTOT1,Label='FA')
  call mma_allocate(puvx,nfint,label='PUVX')
  call mma_allocate(TUVX,NACPR2,Label='TUVX')
  D1I(:) = Zero
  D1A(:) = Zero
  IPR = 0
  if (IPRLOC(2) == DEBUG) IPR = 5
  if (IPRLOC(2) == INSANE) IPR = 10
  call TRACTL2(CMO,PUVX,TUVX,D1I,FI,D1A,FA,IPR,lSquare,ExFac)
  call Put_dArray('TwoEIntegral',PUVX,nFINT)
  call mma_deallocate(PUVX)
  call mma_deallocate(FI)
  call mma_deallocate(FA)
  call mma_deallocate(TUVX)
  call mma_deallocate(D1I)
  call mma_deallocate(D1A)
end if

call Timing(dum1,dum2,time2(2),dum3)
TimeTrans = TimeTrans+time2(2)-time2(1)

! This is where MC-PDFT actually computes the PDFT energy for each state
! only after 500 lines of nothing above...
call mma_allocate(e_states,lroots,label='e_states')
call compute_mcpdft_energy(CMO,e_mcscf,e_states)
call mma_deallocate(CMO)

call Timing(dum1,dum2,time1(1),dum3)

! pdft energy now stored in e_states

if (mcpdft_options%wjob .and. (.not. mcpdft_options%mspdft)) call writejob(e_states,lroots)

if (mcpdft_options%mspdft) then
  Method = 'MSPDFT'
  call replace_diag(heff,e_states,lroots)
  call mspdft_finalize(lroots)
else
  Method = 'MCPDFT'
  ! The following put_* calls should be combined with the MS-PDFT
  ! ones if a gradient call is initiated...
  call put_darray('Last energies',e_states,lroots)
  if (mcpdft_options%grad) then
    call put_dscalar('Last energy',e_states(mcpdft_options%rlxroot))
    call put_iscalar('Relax CASSCF root',mcpdft_options%rlxroot)
    call put_carray('MCLR Root','****************',16)
  end if
  if (iprlev >= TERSE) then
    do state=1,lroots
      call PrintResult(u6,'(6X,A,I3,A,F16.8)','MCPDFT root number',state,' Total energy:',[e_states(state)],1)
    end do
  end if
end if
if ((nhole1 > 0) .or. (nelec3 > 0)) Method(7:8) = '-R'
call Put_cArray('Relax Method',Method,8)

!***********************************************************************
!*****************         Closing up MC-PDFT        *******************
!***********************************************************************

! release SEWARD
call ClsSew()

! Finalize Cholesky information if initialized
if (DoCholesky) then
  call Cho_X_Final(irc)
  if (irc /= 0) then
    write(u6,*) 'MC-PDFT: Cho_X_Final fails with return code ',irc
    write(u6,*) ' Try to recover. Calculation continues.'
  end if
end if

! Release  some memory allocations
call mma_deallocate(FockOcc)
call mma_deallocate(e_mcscf)
call mma_deallocate(e_states)

call StatusLine('MCPDFT:','Finished.')
if (IPRLEV >= 2) write(u6,*)

call Timing(dum1,dum2,time1(2),dum3)
TimeTotal = time1(2)
TimeOutput = TimeOutput+time1(2)-time1(1)
if (IPRLEV >= 3) then
  call PrtTim()
  call FastIO('STATUS')
end if

call close_files_mcpdft()

end subroutine mcpdft
