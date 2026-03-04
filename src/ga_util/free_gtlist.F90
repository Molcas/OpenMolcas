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

subroutine Free_GTList()

use Para_Info, only: Is_Real_Par, nProcs
use TList_Mod, only: GT_Status, iTCnSt
#ifdef _MOLCAS_MPP_
use TList_Mod, only: igaTsk, nTasks
#endif

implicit none

if (.not. GT_Status) return
GT_Status = .false.

iTCnSt = 1
if ((.not. Is_Real_Par()) .or. (nProcs == 1)) return
#ifdef _MOLCAS_MPP_
! create global tasklist...
call GATskL(.false.,nTasks,igaTsk)
#endif

end subroutine Free_GTList
