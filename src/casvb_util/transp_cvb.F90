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

subroutine transp_cvb(a,b,n1,n2)
! Transposes matrix A; A and B may share memory.

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: n1, n2
real(kind=wp) :: a(n1,n2), b(n2,n1)
#include "WrkSpc.fh"
integer(kind=iwp) :: i, i1, iskip, j
integer(kind=iwp), external :: mstackr_cvb

i1 = mstackr_cvb(n2*n1)
iskip = -n2+i1-1
do i=1,n1
  iskip = iskip+n2
  do j=1,n2
    work(j+iskip) = a(i,j)
  end do
end do
call fmove_cvb(work(i1),b,n2*n1)
call mfreer_cvb(i1)

return

end subroutine transp_cvb
