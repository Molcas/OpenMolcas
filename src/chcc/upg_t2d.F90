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

subroutine UpG_T2d(T2,dima,adda)
! upgrade T2 (diagonal)

use Index_Functions, only: nTri_Elem
use chcc_global, only: no, T2c
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dima, adda
real(kind=wp), intent(in) :: T2(nTri_Elem(dima),no,no)
integer(kind=iwp) :: a, ab, j

do j=1,no
  ab = 0
  do a=1,dima
    T2c(adda+a,adda+1:adda+a-1,:,j) = T2(ab+1:ab+a-1,:,j)
    T2c(adda+1:adda+a,adda+a,j,:) = T2(ab+1:ab+a,:,j)
    ab = ab+a
  end do
end do

return

end subroutine UpG_T2d
