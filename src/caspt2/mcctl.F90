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
! Copyright (C) 1998, Per Ake Malmqvist                                *
!               2019, Stefano Battaglia                                *
!***********************************************************************
!--------------------------------------------*
! 1998  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*

subroutine MCCTL(HEFF,NSTATE,JSTATE)

use caspt2_global, only: iPrGlb
use PrintLevel, only: VERBOSE
use caspt2_module, only: E2Corr, mState, NLYRoot
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: NSTATE, JSTATE
real(kind=wp), intent(inout) :: HEFF(NSTATE,NSTATE)
integer(kind=iwp) :: ISTATE
real(kind=wp) :: DVALUE, TOTCPU1, TOTCPU2, TOTWALL1, TOTWALL2
character(len=160) :: string
real(kind=wp), allocatable :: cpu_timing(:), wall_timing(:)

call mma_allocate(cpu_timing,nstate,'timing in mcctl')
call mma_allocate(wall_timing,nstate,'timing in mcctl')
cpu_timing(:) = Zero
wall_timing(:) = Zero

! The ket state is JSTATE.
! Loop over the bra states

do ISTATE=1,NSTATE

  write(string,'(A,I0,A,I0,A,I0)') 'Multistate coupling between state ',ISTATE,' and ',JSTATE,' out of ',NSTATE
  call StatusLine('CASPT2: MCCTL: ',string)
  TOTCPU1 = Zero
  TOTWALL1 = Zero
  TOTCPU2 = Zero
  TOTWALL2 = Zero
  call CWTIME(TOTCPU1,TOTWALL1) !start clock for total time

  if (ISTATE == JSTATE) then
    HEFF(ISTATE,JSTATE) = HEFF(ISTATE,JSTATE)+E2CORR
  else
    ! Compute the effective Hamiltonian:
    call HEFVAL(ISTATE,JSTATE,DVALUE)
    HEFF(ISTATE,JSTATE) = HEFF(ISTATE,JSTATE)+DVALUE
  end if

  call CWTIME(TOTCPU2,TOTWALL2)
  cpu_timing(istate) = TOTCPU2-TOTCPU1
  wall_timing(istate) = TOTWALL2-TOTWALL1

end do

if ((IPRGLB >= VERBOSE) .or. (NLYROOT /= 0)) then
  write(u6,*)
  write(u6,*) 'Hamiltonian Effective Couplings'
  write(u6,*) '-------------------------------'
  write(u6,*)
  write(u6,'(10X,6X,A3,I4,A3)') ' | ',MSTATE(JSTATE),' > '
  do ISTATE=1,NSTATE
    write(u6,'(A3,I4,A3,ES22.14)') ' < ',MSTATE(ISTATE),' | ',HEFF(ISTATE,JSTATE)
  end do
end if

if ((IPRGLB >= VERBOSE) .or. (NLYROOT /= 0)) then
  write(u6,*)
  write(u6,'(A,I4,A)') 'Time spent for multi-state couplings for root ',MSTATE(JSTATE),':'
  write(u6,*) '----------------- CPU TIME  -------- WALL TIME'
  do ISTATE=1,NSTATE
    write(u6,'(A3,I4,A,F18.3,2x,F18.3)') ' < ',MSTATE(ISTATE),' |',cpu_timing(istate),wall_timing(istate)
  end do
end if
call mma_deallocate(cpu_timing)
call mma_deallocate(wall_timing)

end subroutine MCCTL
