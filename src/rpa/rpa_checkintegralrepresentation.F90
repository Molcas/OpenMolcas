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

subroutine RPA_CheckIntegralRepresentation()

use RPA_globals, only: doCD, doDF, doLDF

implicit none

if (.not. (doCD .or. doDF .or. doLDF)) then
  call RPA_Warn(2,'RPA requires CD, DF, or LDF. Conventional integrals not implemented.')
end if

end subroutine RPA_CheckIntegralRepresentation
