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

subroutine LDF_CheckThrs()

use Constants, only: Zero

implicit none
#include "localdf.fh"

if (Thr_Accuracy < Zero) then
  call WarningMessage(2,'LDF: Thr_Accuracy<0')
  call Quit_OnUserError()
end if
if (Thr_LDFPrescreen < Zero) then
  call WarningMessage(2,'LDF: Thr_LDFPrescreen<0')
  call Quit_OnUserError()
end if
Thr_LDFPrescreen = min(Thr_LDFPrescreen,Thr_Accuracy)

return

end subroutine LDF_CheckThrs
