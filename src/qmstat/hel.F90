!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) Anders Ohrn                                            *
!***********************************************************************
!  Hel
!
!> @brief
!>   Couple the electrostatic part of the solvent with the QM-region.
!>   Only include the static part, no polarization at this moment
!> @author A. Ohrn
!>
!> @details
!> (2) The electrostatics.
!>
!> @param[in]  Eint   The static field from the solvent on the QM molecule centers
!> @param[in]  itri   Number of elements in triangular \f$ H \f$-matrix
!> @param[in]  ici    Number of MME-centers
!> @param[in]  ql     MME-charges, obtained from the MME
!> @param[in]  dil    MME-dipoles
!> @param[in]  qqxxyy MME-quadrupoles.
!> @param[out] vmat   The electrostatic part of the solute-solvent interaction matrix
!***********************************************************************

subroutine Hel(Eint,itri,ici,ql,dil,qqxxyy,vmat)

use Constants, only: Zero, Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: itri, ici
real(kind=wp), intent(in) :: Eint(ici,10), Ql(itri,ici), Dil(itri,3,ici), QQxxyy(itri,6,ici)
real(kind=wp), intent(out) :: Vmat(itri)
integer(kind=iwp) :: i, j

! Zeros
Vmat(:) = Zero

! The electrostatic perturbation: <psi_i|V_el|psi_j>
do i=1,itri
  do j=1,ici
    Vmat(i) = Vmat(i)+Eint(j,1)*Ql(i,j)+ &
              Eint(j,2)*Dil(i,1,j)+Eint(j,3)*Dil(i,2,j)+Eint(j,4)*Dil(i,3,j)+ &
              Eint(j,5)*QQxxyy(i,1,j)+Eint(j,7)*QQxxyy(i,3,j)+Eint(j,10)*QQxxyy(i,6,j)+ &
              Eint(j,6)*QQxxyy(i,2,j)*Two+Eint(j,8)*QQxxyy(i,4,j)*Two+Eint(j,9)*QQxxyy(i,5,j)*Two
  end do
end do

return

end subroutine Hel
