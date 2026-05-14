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

subroutine GetTauHlp2(Tau,T1,dima,adda,no,nv)
! Make Tau for aGrp == bGrp

use Index_Functions, only: nTri_Elem
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dima, adda, no, nv
real(kind=wp), intent(inout) :: Tau(nTri_Elem(dima),no,no)
real(kind=wp), intent(in) :: T1(nv,no)
integer(kind=iwp) :: a, ab, i

do i=1,no
  ab = 0
  do a=1,dima
    Tau(ab+1:ab+a,i,:) = Tau(ab+1:ab+a,i,:)+T1(adda+a,i)*T1(adda+1:adda+a,:)
    ab = ab+a
  end do
end do

return

end subroutine GetTauHlp2
