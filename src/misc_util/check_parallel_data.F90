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
!            returns stat=.true. if data are idential
!            oterwise returns stat=.false.
!            value of stat is rank-independent
! if act="S" copy data from masterto all processes

#include "compiler_features.h"
#ifdef _MOLCAS_MPP_

subroutine check_parallel_data(x,n,stat,act)

use Para_Info, only: MyRank, nProcs
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: n
real(kind=wp) :: x(n)
logical(kind=iwp) :: stat
character :: act
#include "WrkSpc.fh"
integer(kind=iwp) :: i, irank, itemp, j, k, ns

stat = .true.
if (nProcs == 1) return

ns = n*nProcs
call getmem('x_prll','ALLO','REAL',itemp,ns)
do i=1,ns
  Work(itemp+i-1) = 0
end do

j = MyRank*n
do i=1,n
  Work(j+itemp+i-1) = x(i)
end do
call GADsum(Work(itemp),ns)

if (act == 'C') then
  do irank=0,nProcs-1
    ! comparing MyRank vs rank j
    if (irank == MyRank) cycle
    i = MyRank*n
    j = irank*n
    do k=1,n
      if (Work(i+itemp+k-1) /= Work(j+itemp+k-1)) then
        stat = .false.
        exit
      end if
    end do
  end do

else if (act == 'S') then
  ! copy data from master to MyRank
  if (MyRank /= 0) then
    j = 0
    do i=1,n
      x(i) = Work(j+itemp+i-1)
    end do
  end if

else
  write(u6,*) 'check_parallel_data(), illegal value:'
  write(u6,*) 'act=',act
  write(u6,*) 'correct function call!!'
  call abort()
end if

call getmem('x_prll','FREE','REAL',itemp,ns)

return

end subroutine check_parallel_data

#elif !defined (EMPTY_FILES)

! Some compilers do not like empty files
#include "macros.fh"
dummy_empty_procedure(check_parallel_data)

#endif
