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

subroutine noperm(wrk,wrksize,a,b,posst)
! realize mapping without permutation
! define %d, %i

use ccsd_global, only: Map_Type, nsym
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: wrksize, posst
real(kind=wp) :: wrk(wrksize)
type(Map_Type) :: a, b
integer(kind=iwp) :: i, ib, j, k, nhelp

! def b%i

do k=1,nsym
  do j=1,nsym
    do i=1,nsym
      b%i(i,j,k) = a%i(i,j,k)
    end do
  end do
end do

! def initial values

do nhelp=1,6
  b%d(0,nhelp) = a%d(0,nhelp)
end do

posst = b%pos0
do ib=1,a%d(0,5)
  do nhelp=2,6
    b%d(ib,nhelp) = a%d(ib,nhelp)
  end do
  b%d(ib,1) = posst
  posst = posst+b%d(ib,2)

  call map11(wrk(a%d(ib,1)),wrk(b%d(ib,1)),a%d(ib,2),1)

end do

return

end subroutine noperm
