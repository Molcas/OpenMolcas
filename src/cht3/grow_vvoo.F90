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

subroutine grow_vvoo(A,B,no,nv,dima,dimb,lasta,lastb)
! this routine does:
!
! grow A(1324)/(vvoo) from the blocked cholesky vectors B(1234)/(vovo)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: no, nv, dima, dimb, lasta, lastb
real(kind=wp), intent(inout) :: A(nv,nv,no,no)
real(kind=wp), intent(in) :: B(dima,no,dimb,no)
integer(kind=iwp) :: i2, i3, i4

!!write(u6,'(A,4(i3,x))') 'AA lasta, lastb, dima, dimb ',lasta,lastb,dima,dimb

do i4=1,no
  do i3=1,no
    do i2=1,dimb
      A(lasta+1:lasta+dima,lastb+i2,i3,i4) = B(:,i3,i2,i4)
    end do
  end do
end do

return

end subroutine grow_vvoo
