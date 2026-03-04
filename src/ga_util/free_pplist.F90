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

subroutine Free_PPList()

use TList_Mod, only: PP_Status, TskL
use Para_Info, only: Is_Real_Par, nProcs
use stdalloc, only: mma_deallocate

implicit none

if (.not. allocated(TskL)) return
PP_Status = .false.

if ((.not. Is_Real_Par()) .or. (nProcs == 1)) return
call mma_deallocate(TskL)

end subroutine Free_PPList
