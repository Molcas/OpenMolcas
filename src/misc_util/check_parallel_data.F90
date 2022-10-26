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

! if action="C" checks if array x(N) is identical across processes
!               returns status=.true. if data are idential
!               oterwise returns status=.false.
!                value of status is rank-independent
! if action="S" copy data from masterto all processes

#include "compiler_features.h"
#ifdef _MOLCAS_MPP_

subroutine check_parallel_data(x,n,status,action)

use Para_Info, only: MyRank, nProcs

implicit none
integer :: n
real*8 :: x(n)
logical :: status
integer :: itemp, ns, j, i, k, irank
character :: action
#include "WrkSpc.fh"

status = .true.
if (nProcs == 1) return

ns = n*nProcs
call getmem('x_prll','ALLO','REAL',itemp,ns)
do i=1,ns
  Work(itemp+i-1) = 0
end do
!
j = MyRank*n
do i=1,n
  Work(j+itemp+i-1) = x(i)
end do
call GADsum(Work(itemp),ns)
!
if (action == 'C') then
  do irank=0,nProcs-1
!     comparing MyRank vs rank j
    if (irank == MyRank) cycle
    i = MyRank*n
    j = irank*n
    do k=1,n
      if (Work(i+itemp+k-1) /= Work(j+itemp+k-1)) then
        status = .false.
        exit
      end if
    end do
  end do

else if (action == 'S') then
!     copy data from master to MyRank
  if (MyRank /= 0) then
    j = 0
    do i=1,n
      x(i) = Work(j+itemp+i-1)
    end do
  end if
!
else
  write(6,*) 'check_parallel_data(), illegal value:'
  write(6,*) 'action=',action
  write(6,*) 'correct function call!!'
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
