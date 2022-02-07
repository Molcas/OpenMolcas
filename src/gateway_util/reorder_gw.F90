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

subroutine Reorder_GW(A,B,k,l,n,m)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: k, l, n, m
real(kind=wp), intent(in) :: A(k,l,n,m)
real(kind=wp), intent(out) :: B(k,n,l,m)
integer(kind=iwp) :: l_, n_, m_

do m_=1,m
  do n_=1,n
    do l_=1,l
      B(:,n_,l_,m_) = A(:,l_,n_,m_)
    end do
  end do
end do

return

end subroutine Reorder_GW
