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

subroutine Do_Index(Indx,NrBas,NrBas_Eff,iCmp)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: NrBas_Eff, iCmp, NrBas
integer(kind=iwp), intent(out) :: Indx(NrBas_Eff,iCmp)
integer(kind=iwp) :: iAdd, iB, iB_Eff, iC

iAdd = NrBas-NrBas_Eff
do iB_Eff=1,NrBas_Eff
  iB = iB_Eff+iAdd
  do iC=1,iCmp
    Indx(iB_Eff,iC) = (iC-1)*NrBas+iB
  end do
end do

return

end subroutine Do_Index
