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

subroutine Cho_P_OpenR(iOpt)
!
! Purpose: open (iOpt=1) or close (iOpt=2) files for storing global
!          reduced set indices.

use Cholesky, only: LuRed_G
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: iOpt
character(len=5) :: FNRed
character(len=*), parameter :: SecNam = 'Cho_P_OpenR'

if (iOpt == 1) then
  LuRed_G = 7
  FNRed = 'CHRED'
  call DAName_MF_WA(LuRed_G,FNRed)
else if (iOpt == 2) then
  if (LuRed_G > 0) call DAClos(LuRed_G)
else
  call Cho_Quit('iOpt error in '//SecNam,104)
end if

end subroutine Cho_P_OpenR
