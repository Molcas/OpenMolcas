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

subroutine noperm(wrk,wrksize,a,b,post)
! realize mapping without permutation
! define %d, %i

use ccsd_global, only: Map_Type, nsym
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: wrksize
real(kind=wp), intent(inout) :: wrk(wrksize)
type(Map_Type), intent(in) :: a
type(Map_Type), intent(inout) :: b
integer(kind=iwp), intent(out) :: post
integer(kind=iwp) :: ib

! def b%i

b%i(1:nsym,1:nsym,1:nsym) = a%i(1:nsym,1:nsym,1:nsym)

! def initial values

b%d(0,:) = a%d(0,:)

post = b%pos0
do ib=1,a%d(0,5)
  b%d(ib,2:6) = a%d(ib,2:6)
  b%d(ib,1) = post
  post = post+b%d(ib,2)

  call map11(wrk(a%d(ib,1)),wrk(b%d(ib,1)),a%d(ib,2),1)

end do

return

end subroutine noperm
