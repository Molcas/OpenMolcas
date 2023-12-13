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
! Copyright (C) 1986, Per E. M. Siegbahn                               *
!***********************************************************************

subroutine JTIME(ITIME)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(out) :: ITIME
real(kind=wp) :: DTIM, DUM1, DUM2, DUM3

call TIMING(DTIM,DUM1,DUM2,DUM3)
ITIME = int(DTIM)

return

end subroutine JTIME
