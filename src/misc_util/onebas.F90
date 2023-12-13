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

subroutine OneBas(Label)
!***********************************************************************
!                                                                      *
!     Change nBas in OneDat                                            *
!                                                                      *
!***********************************************************************

use OneDat, only: nBas, nSym
use Definitions, only: iwp, u6

implicit none
character(len=*), intent(in) :: Label
integer(kind=iwp) :: IntBas(8)

if (Label == 'CONT') then
  call Get_iArray('nBas',IntBas,nSym)
else if (Label == 'PRIM') then
  call Get_iArray('nBas_Prim',IntBas,nSym)
else
  write(u6,*) 'OneBas: Illegal Label value!'
  write(u6,*) 'Value: ',Label
  call Abend()
end if
nBas(1:nSym) = IntBas(1:nSym)

return

end subroutine OneBas
