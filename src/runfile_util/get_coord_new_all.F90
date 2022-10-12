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
! Copyright (C) Roland Lindh                                           *
!***********************************************************************
!  Get_Coord_New_All
!
!> @brief
!>   Get the new coordinates from RUNFILE
!> @author R. Lindh
!>
!> @details
!> Place Cartesian coordinates (in a.u.) into array \p Coord_All(3,*).
!>
!> @param[out] Coord_All  Array of coordinates
!> @param[in]  nAtoms_All Number of atoms
!***********************************************************************

subroutine Get_Coord_New_All(Coord_All,nAtoms_All)

use RunFile_procedures, only: Get_Coord_New
use stdalloc, only: mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nAtoms_All
real(kind=wp), intent(out) :: Coord_All(3,nAtoms_All)
integer(kind=iwp) :: nAtoms, nAtoms_Allx
real(kind=wp), allocatable :: CU(:,:)

call Get_nAtoms_All(nAtoms_Allx)
if (nAtoms_All /= nAtoms_Allx) then
  write(u6,*) 'Get_Coord_New_All: nAtoms_All /= nAtoms_Allx'
  write(u6,*) 'nAtoms_All=',nAtoms_All
  write(u6,*) 'nAtoms_Allx=',nAtoms_Allx
  call Abend()
end if
call Get_Coord_New(CU,nAtoms)
call Get_Coord_All_(CU,nAtoms,Coord_All,nAtoms_All)
call mma_deallocate(CU)

return

end subroutine Get_Coord_New_All
