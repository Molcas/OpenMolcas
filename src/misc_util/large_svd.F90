!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine large_svd(m,n,amat,umat,vmat,svals)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: m, n
real(kind=wp) :: amat(m,*), umat(m,*), vmat(n,*), svals(*)
#include "WrkSpc.fh"
integer(kind=iwp) :: info, ipwork, lwork, nm
real(kind=wp) :: wrk1_lapack(1)

! Note that dgesvd returns V**T, not V.
nm = min(n,m)
!write(u6,*) ' In large_svd. Calling dgesvd:'
call dgesvd_('S','S',m,n,amat,m,svals,umat,m,vmat,nm,wrk1_lapack,-1,info)
!write(u6,*) ' large_svd back from dgesvd'
lwork = int(wrk1_lapack(1))
!write(u6,*) ' lwork:',lwork
call getmem('lapckwrk','allo','real',ipwork,lwork)
!write(u6,*) ' Calling dgesvd again:'
call dgesvd_('S','S',m,n,amat,m,svals,umat,m,vmat,nm,work(ipwork),lwork,info)
!write(u6,*) ' large_svd back from dgesvd'
call getmem('lapckwrk','free','real',ipwork,lwork)

return

end subroutine large_svd
