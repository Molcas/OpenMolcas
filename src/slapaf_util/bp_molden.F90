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
! Copyright (C) 2026, Ignacio Fdez. Galvan                             *
!***********************************************************************

subroutine BP_Molden(nAtom,x_vector,y_vector,Asymmetry)

use Symmetry_Info, only: nIrrep
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nAtom
real(kind=wp), intent(in) :: x_vector(3,nAtom), y_vector(3,nAtom), Asymmetry
integer(kind=iwp) :: iAtom, LuBP
real(kind=wp), allocatable :: Coord(:,:)
character(len=2), allocatable :: Element(:)
integer(kind=iwp), external :: isFreeUnit

! Only without symmetry
if (nIrrep /= 1) return

!                                                                      *
!***********************************************************************
!                                                                      *
! Open input file for MOLDEN

LuBP = isFreeUnit(9)
call molcas_open(LuBP,'MD_BP')
!                                                                      *
!***********************************************************************
!                                                                      *
write(LuBP,*) '[Molden Format]'
!                                                                      *
!***********************************************************************
!                                                                      *
! Write fake frequencies to Molden input file.
! The frequencies are proportional to the size of the NAC vector along that direction

write(LuBP,*) '[N_FREQ]'
write(LuBP,*) 2
write(LuBP,*) '[FREQ]'
write(LuBP,*) sqrt((One-Asymmetry**2)/(One-Asymmetry))
write(LuBP,*) sqrt((One-Asymmetry**2)/(One+Asymmetry))
!                                                                      *
!***********************************************************************
!                                                                      *
! Write coordinates of all centers

call mma_allocate(Coord,3,nAtom,label='Coord')
call mma_allocate(Element,nAtom,label='Element')
call Get_Coord_All(Coord,nAtom)
call Get_Name_All(Element)
write(LuBP,*) '[NATOM]'
write(LuBP,*) nAtom

write(LuBP,*) '[FR-COORD]'
do iAtom=1,nAtom
  write(LuBP,'(A,3E16.8)') Element(iAtom),Coord(:,iAtom)
end do
call mma_deallocate(Coord)
call mma_deallocate(Element)
!                                                                      *
!***********************************************************************
!                                                                      *
! Write the branching plane vectors

write(LuBP,*) '[FR-NORM-COORD]'
write(LuBP,*) 'vibration ',1
do iAtom=1,nAtom
  write(LuBP,*) x_vector(:,iAtom)
end do
write(LuBP,*) 'vibration ',2
do iAtom=1,nAtom
  write(LuBP,*) y_vector(:,iAtom)
end do

close(LuBP)

end subroutine BP_Molden
