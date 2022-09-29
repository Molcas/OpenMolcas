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

subroutine grow_l2(A,B,nc,nv,dima,dimb,lasta,lastb)
! this routine does:
!
! grow A(A,B,m) from the blocked cholesky vectors B(a',b',m)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nc, nv, dima, dimb, lasta, lastb
real(kind=wp), intent(inout) :: A(nv,nv,nc)
real(kind=wp), intent(in) :: B(dima,dimb,nc)
integer(kind=iwp) :: i2, i3

do i3=1,nc
  do i2=1,dimb
    A(lasta+1:lasta+dima,lastb+i2,i3) = B(:,i2,i3)
    A(lastb+i2,lasta+1:lasta+dima,i3) = B(:,i2,i3)
  end do
end do

return

end subroutine grow_l2
