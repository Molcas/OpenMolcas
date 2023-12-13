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

subroutine NQGrid_Init()

use Grid_On_Disk, only: Final_Grid, G_S, iDisk_Grid, iDisk_Set, Intermediate, Lu_Grid, Old_Functional_Type, Regenerate
use nq_Info, only: Other_Type
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: iDisk, iDum(1)

!                                                                      *
!***********************************************************************
!                                                                      *
! Make the grid file dirty

! Open the file.
Lu_Grid = 77
call DaName_MF_WA(Lu_Grid,'NQGRID')

! Write the status flag and disk addresses fo the sets.

iDisk_Set(Final_Grid) = -1
iDisk_Set(Intermediate) = -1
G_S(Final_Grid) = Regenerate
G_S(Intermediate) = Regenerate
Old_Functional_Type = Other_Type

iDisk_Grid = 0
call iDaFile(Lu_Grid,1,G_S,2,iDisk_Grid)
iDisk = iDisk_Grid
call iDaFile(Lu_Grid,1,iDisk_Set,2,iDisk_Grid)
iDum(1) = Old_Functional_Type
call iDaFile(Lu_Grid,1,iDum,1,iDisk_Grid)

iDisk_Set(Final_Grid) = iDisk_Grid
iDisk_Set(Intermediate) = iDisk_Grid

call iDaFile(Lu_Grid,1,iDisk_Set,2,iDisk)

call DaClos(Lu_Grid)
!                                                                      *
!***********************************************************************
!                                                                      *

return

end subroutine NQGrid_Init
