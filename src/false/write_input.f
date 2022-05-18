************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2020, Ignacio Fdez. Galvan                             *
************************************************************************

      subroutine Write_Input()
      use stdalloc, only: mma_Allocate, mma_Deallocate
      implicit none
      integer :: LU,nAtoms,i
      real*8, allocatable :: Coord(:,:)
      character(len=2), allocatable :: Symbol(:)
      integer, external :: IsFreeUnit
#include "constants2.fh"

      ! read data from RUNFILE
      call Get_nAtoms_All(nAtoms)
      call mma_Allocate(Coord,3,nAtoms,label='Coord')
      call mma_Allocate(Symbol,nAtoms,label='Symbol')
      call Get_Coord_All(Coord,nAtoms)
      call Get_Name_All(Symbol)

      ! write interface input file
      LU=IsFreeUnit(11)
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

100   format(A)
101   format(I6)
102   format(A2,1X,3F20.12)
      end subroutine Write_Input
