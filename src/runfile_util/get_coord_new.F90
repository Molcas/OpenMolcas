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
!  Get_Coord_New
!
!> @brief
!>   Get the updated/new symmetry-unique Cartesian coordinates of the basis set centers
!> @author R. Lindh
!>
!> @details
!> The utility will read the updated/new symmetry-unique Cartesian coordinates of the basis set centers from the runfile.
!>
!> @param[out] CN     Array of the symmetry-unique Cartesian coordinates of the basis set centers
!> @param[out] nAtoms Number of symmetry-unique basis set centers
!***********************************************************************

#ifdef _IN_MODULE_

subroutine Get_Coord_New(CN,nAtoms)

use stdalloc, only: mma_allocate
use Definitions, only: wp, iwp

implicit none
real(kind=wp), allocatable, intent(out) :: CN(:,:)
integer(kind=iwp), intent(out) :: nAtoms
integer(kind=iwp) :: nAtoms3
logical(kind=iwp) :: Found
character(len=*), parameter :: Label = 'GeoNew'

call qpg_dArray(Label,Found,nAtoms3)
nAtoms = nAtoms3/3
if ((.not. Found) .or. (nAtoms3 == 0)) return
call mma_allocate(CN,3,nAtoms)
call Get_dArray(Label,CN,nAtoms3)

return

end subroutine Get_Coord_New

#endif
