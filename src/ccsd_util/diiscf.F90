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

subroutine diiscf(diispoint,ndiis)
! this routine closes 1-ndiis Temp files for DIIS procedure
! lun's are stored in stack (diispoint)
!
! diispoint - stack of lun numbers (I)
! ndiis     - size of diis procedure (I)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: diispoint(4), ndiis
integer(kind=iwp) :: lun, rc

if (ndiis > 0) then
  lun = diispoint(1)
  call filemanager(3,lun,rc)
end if

if (ndiis > 1) then
  lun = diispoint(2)
  call filemanager(3,lun,rc)
end if

if (ndiis > 2) then
  lun = diispoint(3)
  call filemanager(3,lun,rc)
end if

if (ndiis > 3) then
  lun = diispoint(4)
  call filemanager(3,lun,rc)
end if

return

end subroutine diiscf
