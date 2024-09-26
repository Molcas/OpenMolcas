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

module Save_Stuff

use Definitions, only: wp

implicit none
private

real(kind=wp) :: DltNTh_old, DThr_Old, EThr_old, FThr_Old, ThrInt_old

public :: DltNTh_old, DThr_Old, EThr_old, FThr_Old, ThrInt_old

end module Save_Stuff
