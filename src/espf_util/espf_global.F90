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

module espf_global

use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp) :: MMIterMax
real(kind=wp) :: ConvF
integer(kind=iwp), parameter :: MxExtPotComp = 10, &
                                QM = 0, MMI = 1, MMO = 2
character(len=*), parameter :: TPRDefName = 'topol.tpr'

public :: ConvF, MMI, MMIterMax, MMO, MxExtPotComp, QM, TPRDefName

end module espf_global
