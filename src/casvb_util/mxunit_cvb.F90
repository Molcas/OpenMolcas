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

!IFG trivial
subroutine mxunit_cvb(a,n)

use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: n
real(kind=wp) :: a(n,n)
integer(kind=iwp) :: i

call fzero(a,n*n)
do i=1,n
  a(i,i) = One
end do

return

end subroutine mxunit_cvb
