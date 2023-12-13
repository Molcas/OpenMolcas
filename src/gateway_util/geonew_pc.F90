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
! Copyright (C) 1991,2003, Roland Lindh                                *
!***********************************************************************

subroutine GeoNew_PC()
!***********************************************************************
!                                                                      *
! Object: to pick up the geometry from a special file. This will only  *
!         make any difference of there exist a file otherwise SEWARD   *
!         will use the geometry as specified by the standard input     *
!         file.                                                        *
!                                                                      *
!     Author: Roland Lindh, Dep. of Theoretical Chemistry,             *
!             University of Lund, SWEDEN                               *
!             March '91                                                *
!                                                                      *
!     Modified to work with point charges. RL 20030507                 *
!***********************************************************************

use RunFile_procedures, only: Get_PC_Coord_New
use External_Centers, only: nData_XF, XF
use stdalloc, only: mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: lBuf, nAtoms
real(kind=wp), allocatable :: CN(:)

! Check if there is a data field called 'GeoNewPC'

call Get_PC_Coord_New(CN,lBuf)
nAtoms = lbuf/nData_XF

! Quit if the datafield 'NewGeom' is not available

if (lBuf == 0) then
  return
end if

! Replace coodinates read in subroutine input

call dcopy_(nAtoms*nData_XF,CN,1,XF,1)
write(u6,*)
write(u6,'(A)') '    Point Charge data read from RUNFILE'
write(u6,*)

! Epilogue, end

call mma_deallocate(CN)

return

end subroutine GeoNew_PC
