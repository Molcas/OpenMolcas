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
! Copyright (C) 2020, Ignacio Fdez. Galvan                             *
!***********************************************************************

subroutine Write_Input()

use stdalloc, only: mma_Allocate, mma_Deallocate
use Definitions, only: wp, iwp
use Constants, only: Angstrom

implicit none
integer(kind=iwp) :: LU, nAtoms, i
real(kind=wp), allocatable :: Coord(:,:)
character(len=2), allocatable :: Symbol(:)
integer(kind=iwp), external :: IsFreeUnit

! read data from RUNFILE
call Get_nAtoms_All(nAtoms)
call mma_Allocate(Coord,3,nAtoms,label='Coord')
call mma_Allocate(Symbol,nAtoms,label='Symbol')
call Get_Coord_All(Coord,nAtoms)
call Get_Name_All(Symbol)

! write interface input file
LU = IsFreeUnit(11)
call Molcas_Open(LU,'INPUT')
write(LU,100) '[XYZ]'
write(LU,101) nAtoms
write(LU,100) 'angstrom'
do i=1,nAtoms
  write(LU,102) Symbol(i),Angstrom*Coord(:,i)
end do
close(LU)

! clean up
call mma_Deallocate(Coord)
call mma_Deallocate(Symbol)

return

100 format(a)
101 format(i6)
102 format(a2,1x,3f20.12)

end subroutine Write_Input
