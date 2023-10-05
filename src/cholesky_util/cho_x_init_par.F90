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

subroutine Cho_X_Init_Par(irc,isDF)
!
! Purpose: setup for parallel Cholesky/DF.

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: irc
logical(kind=iwp), intent(in) :: isDF

if (isDF) then
  call Cho_X_Init_Par_DF(irc)
else
  call Cho_X_Init_Par_GenBak()
  call Cho_X_Init_Par_Cho(irc)
end if

end subroutine Cho_X_Init_Par
