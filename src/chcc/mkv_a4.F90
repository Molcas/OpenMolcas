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

subroutine MkV_A4(Vp,V,dimb,dima,no,dimij)
! this routine does:
! Vp(a,b,ij) <- (ai|bj) from V(b,j,a,i)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dimb, dima, no, dimij
real(kind=wp), intent(out) :: Vp(dima,dimb,dimij)
real(kind=wp), intent(in) :: V(dimb,no,dima,no)
integer(kind=iwp) :: b, i, ij, j

ij = 0
do i=1,no
  do j=1,i
    ij = ij+1
    do b=1,dimb
      Vp(:,b,ij) = V(b,j,:,i)
    end do
  end do
end do

return

end subroutine MkV_A4
