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

subroutine GetCnt(NATOMS)
! Purpose: Read data from ONEINT

use Cntrl, only: Coor
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: NATOMS

! Read NATOMS, nr of symmetry-unique atoms
call Get_iScalar('Unique atoms',NATOMS)

! Read COOR, cartesian coordinates of each atom
call Get_dArray('Unique Coordinates',COOR,3*nAtoms)

end subroutine GetCnt
