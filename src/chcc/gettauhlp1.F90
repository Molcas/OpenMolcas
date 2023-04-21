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

subroutine GetTauHlp1(Tau,T1,dima,dimb,adda,addb,no,nv)
! Make Tau for aGrp /= bGrp

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dima, dimb, adda, addb, no, nv
real(kind=wp), intent(inout) :: Tau(dima,dimb,no,no)
real(kind=wp), intent(in) :: T1(nv,no)
integer(kind=iwp) :: b, j

do j=1,no
  do b=1,dimb
    Tau(:,b,:,j) = Tau(:,b,:,j)+T1(addb+b,j)*T1(adda+1:adda+dima,:)
  end do
end do

return

end subroutine GetTauHlp1
