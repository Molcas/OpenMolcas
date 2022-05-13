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

function nBas_Eff(NrExp,NrBas,Cff,nExp_Eff)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: nBas_Eff
integer(kind=iwp), intent(in) :: NrExp, NrBas, nExp_Eff
real(kind=wp), intent(in) :: Cff(NrExp,NrBas)
integer(kind=iwp) :: iBas, iExp

nBas_Eff = NrBas

do iBas=1,NrBas

  do iExp=1,nExp_Eff

    if (Cff(iExp,iBas) /= Zero) then
      nBas_Eff = NrBas-iBas+1
      return
    end if

  end do

end do

return

end function nBas_Eff
