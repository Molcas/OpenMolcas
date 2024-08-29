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

subroutine ResortD(D_Old,D_New,iBas,iCmp,jBas,jCmp)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iBas, iCmp, jBas, jCmp
real(kind=wp), intent(in) :: D_Old(iBas,jBas,iCmp,jCmp)
real(kind=wp), intent(out) :: D_New(iBas,iCmp,jBas,jCmp)
integer(kind=iwp) :: iC, jB, jC

do jC=1,jCmp
  do jB=1,jBas
    do iC=1,iCmp
      D_New(:,iC,jB,jC) = D_Old(:,jB,iC,jC)
    end do
  end do
end do

return

end subroutine ResortD
