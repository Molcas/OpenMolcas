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

subroutine ReInit_GTList()

use TList_Mod, only: GT_Status, iTCnSt
#ifdef _MOLCAS_MPP_
use TList_Mod, only: iGATsk
#endif
use Para_Info, only: Is_Real_Par, nProcs
use Definitions, only: u6

implicit none

if (.not. GT_Status) then
  write(u6,*) 'ReInit_GTList: List not active!'
  call Abend()
end if
iTCnSt = 1
if ((.not. Is_Real_Par()) .or. (nProcs == 1)) return
#ifdef _MOLCAS_MPP_
! initialize global tasklist...
call GATskL_Zero(igaTsk)
#endif

end subroutine ReInit_GTList
