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

subroutine Fix_Coeff(nPrim,nCntrc,Coeff_c,Coeff_p,Mode)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nPrim, nCntrc
real(kind=wp), intent(inout) :: Coeff_c(nPrim,nCntrc)
real(kind=wp), intent(in) :: Coeff_p(nPrim,nPrim)
character, intent(in) :: Mode
integer(kind=iwp) :: iP

! Put in the normalization constant for the product basis function.

if (Mode == 'F') then
  do iP=1,nPrim
    Coeff_c(iP,:) = Coeff_c(iP,:)/Coeff_p(iP,iP)
  end do
else
  do iP=1,nPrim
    Coeff_c(iP,:) = Coeff_c(iP,:)*Coeff_p(iP,iP)
  end do
end if

return

end subroutine Fix_Coeff
