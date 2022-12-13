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
! Copyright (C) 2018, Denis Jelovina                                   *
!***********************************************************************

! if act="C" checks if array x(N) is identical across processes
!            returns stat=.true. if data are identical
!            oterwise returns stat=.false.
!            value of stat is rank-independent
! if act="S" copy data from master to all processes

#include "compiler_features.h"
#ifdef _MOLCAS_MPP_

subroutine check_parallel_data(x,n,stat,act)

use Para_Info, only: MyRank, nProcs
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: n
real(kind=wp), intent(inout) :: x(n)
logical(kind=iwp), intent(out) :: stat
character, intent(in) :: act
integer(kind=iwp) :: irank, k, this
real(kind=wp), allocatable :: x_prll(:,:)

stat = .true.
if (nProcs == 1) return

call mma_allocate(x_prll,n,nProcs,label='x_prll')
x_prll(:,:) = Zero

this = MyRank+1
x_prll(:,this) = x
call GADsum(x_prll,n*nProcs)

if (act == 'C') then
  do irank=1,nProcs
    ! comparing MyRank vs rank irank
    if (irank == this) cycle
    do k=1,n
      if (x_prll(k,irank) /= x_prll(k,this)) then
        stat = .false.
        exit
      end if
    end do
  end do

else if (act == 'S') then
  ! copy data from master to MyRank
  if (this /= 1) x(:) = x_prll(:,1)

else
  write(u6,*) 'check_parallel_data(), illegal value:'
  write(u6,*) 'act=',act
  write(u6,*) 'correct function call!!'
  call abort()
end if

call mma_deallocate(x_prll)

return

end subroutine check_parallel_data

#elif !defined (EMPTY_FILES)

! Some compilers do not like empty files
#include "macros.fh"
dummy_empty_procedure(check_parallel_data)

#endif
