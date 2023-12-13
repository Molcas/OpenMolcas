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
! Copyright (C) 2013, Ignacio Fdez. Galvan                             *
!***********************************************************************
!  Fix_Symmetry
!
!> @brief
!>   Fix the symmetry of a structure by removing off-zero values
!> @author Ignacio Fdez. Galv&aacute;n
!>
!> @details
!> Fix the symmetry of a structure, by making sure that the coordinates that
!> should be zero remain there, removing numerical inaccuracies.
!>
!> @param[in,out] Coord Cartesian coordinates of the structure to fix (unique atoms)
!> @param[in]     nAt   Number of atoms in the structure
!> @param[in]     Stab  Stabilizers for each atom
!***********************************************************************

subroutine Fix_Symmetry(Coord,nAt,Stab)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nAt, Stab(nAt)
real(kind=wp), intent(inout) :: Coord(3,nAt)
integer(kind=iwp) :: iAt, j
real(kind=wp), parameter :: thr = 1.0e-12_wp

do iAt=1,nAt
  do j=1,3
    if (btest(Stab(iAt),j-1)) then
      if (abs(Coord(j,iAt)) > thr) call WarningMessage(1,'Significant deviation from symmetry axis.')
      Coord(j,iAt) = Zero
    end if
  end do
end do

end subroutine Fix_Symmetry
