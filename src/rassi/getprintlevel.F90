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

subroutine GETPRINTLEVEL()

use rassi_aux, only: ipglob
use Definitions, only: iwp

implicit none
integer(kind=iwp), external :: iPrintLevel
logical(kind=iwp), external :: Reduce_Prt

! Global print level
IPGLOB = iPrintLevel(-1)
! Change print level when in an optimization loop
if (Reduce_Prt()) IPGLOB = max(IPGLOB-2,0)

end subroutine GETPRINTLEVEL
