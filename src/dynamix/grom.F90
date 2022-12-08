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
! Copyright (C) 2008, Igor Schapiro                                    *
!***********************************************************************
!
! *********************************************************************
! *                                                                   *
! *  Writes out the forces and energies for Gromacs                   *
! *                                                                   *
! * 18/01/2008                                                        *
! * Igor Schapiro                                                     *
! *                                                                   *
! *********************************************************************

subroutine GROM(irc)

use Dynamix_Globals, only: iPrint, INSANE
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: irc
integer(kind=iwp) :: natom, i, j, filenum
character(len=80) :: filename
real(kind=wp), allocatable :: xyz(:), force(:)
character(len=2), allocatable :: atom(:)
integer(kind=iwp), external :: IsFreeUnit
#include "warnings.h"

if (IPRINT == INSANE) write(u6,*) ' Entering GROM'

write(u6,*) '**** Writes out Forces and Energies for Gromacs ****'

call DxRdNAtomStnd(natom)
call mma_allocate(atom,natom)
call mma_allocate(xyz,natom*3)
call mma_allocate(force,natom*3)

! Read atom, their coordinates and forces

call DxRdStnd(natom,atom,xyz,force)

! Write the energies and forces to file

filenum = IsFreeUnit(81)
filename = 'MOL2GROM'
call Molcas_Open(filenum,filename)
write(filenum,*) natom
do i=1,natom
  write(filenum,'(3es20.10)') (force((i-1)*3+j),j=1,3)
end do
close(filenum)

call mma_deallocate(atom)
call mma_deallocate(xyz)
call mma_deallocate(force)

irc = _RC_ALL_IS_WELL_

return

end subroutine GROM
