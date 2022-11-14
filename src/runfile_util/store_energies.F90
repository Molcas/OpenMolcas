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

subroutine Store_Energies(nEnergies,Energies,iEnergy)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nEnergies, iEnergy
real(kind=wp), intent(in) :: Energies(nEnergies)

call Put_iScalar('Number of roots',nEnergies)
call Put_dArray('Last energies',Energies,nEnergies)
call Put_dScalar('Last energy',Energies(iEnergy))

return

end subroutine Store_Energies
