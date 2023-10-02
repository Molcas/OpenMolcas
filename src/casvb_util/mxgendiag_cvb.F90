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
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

subroutine mxgendiag_cvb(a,s,eigval,n)

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: n
real(kind=wp) :: a(n,n), s(n,n), eigval(n)
integer(kind=iwp) :: info, lwrk
real(kind=wp) :: wrk(1)
real(kind=wp), allocatable :: wrk2(:)

info = 0
lwrk = -1
call dsygv_(1,'V','U',n,a,n,s,n,eigval,wrk,lwrk,info)
lwrk = nint(wrk(1))
call mma_allocate(wrk2,lwrk,label='wrk2')
call dsygv_(1,'V','U',n,a,n,s,n,eigval,wrk2,lwrk,info)
call mma_deallocate(wrk2)
if (info /= 0) then
  write(u6,*) ' Error in generalized diagonalization!'
  write(u6,*) ' Dsygv exited with code:',info
  call abend_cvb()
end if

return

end subroutine mxgendiag_cvb
