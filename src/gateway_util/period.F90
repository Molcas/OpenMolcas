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

module Period
! mAdCell, lthCell - information on atoms in the unit cell
! If iAt = 1..lthCell, then
!    iWork(mAdCell+(iAt-1)) - total sequence # of the cell's atom in the full list
! VCell   - unit cell vectors
! ispread - how many steps to spread the unit cell along each unit cell vector

integer mAdCell, lthCell, nAtC
equivalence(lthCell,nAtC)

integer ispread(3)
real*8 VCell(3,3)
logical Cell_l
integer, allocatable :: AdCell(:)

end module Period
