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

subroutine Init_GTList()

use TList_Mod, only: GT_Status, iTCnSt
#ifdef _MOLCAS_MPP_
use TList_Mod, only: nTasks, igaTsk
#endif
use Para_Info, only: nProcs, Is_Real_Par

implicit none

if (GT_Status) return
GT_Status = .true.

iTCnSt = 1
if ((.not. Is_Real_Par()) .or. (nProcs == 1)) return
#ifdef _MOLCAS_MPP_
! create global tasklist...
call GATskL(.true.,nTasks,igaTsk)
#endif

end subroutine Init_GTList
