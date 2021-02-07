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

#include "warnings.fh"
#include "Molcas.fh"
#include "prgm.fh"
#include "stdalloc.fh"
parameter(ROUTINE='GROM')
#include "MD.fh"
#include "WrkSpc.fh"
external IsFreeUnit
integer natom, i, j, irc, file, IsFreeUnit
character filname*80
real*8, allocatable :: xyz(:), force(:)
character, allocatable :: atom(:)*2

if (IPRINT == INSANE) write(6,*) ' Entering ',ROUTINE

write(6,*) '**** Writes out Forces and Energies for Gromacs ****'

call DxRdNAtomStnd(natom)
call mma_allocate(atom,natom)
call mma_allocate(xyz,natom*3)
call mma_allocate(force,natom*3)

! Read atom, their coordinates and forces

call DxRdStnd(natom,atom,xyz,force)

! Write the energies and forces to file

file = IsFreeUnit(81)
filname = 'MOL2GROM'
call Molcas_Open(file,filname)
write(file,*) natom
do i=1,natom
  write(file,'(3D20.10)') (force((i-1)*3+j),j=1,3)
end do
close(file)

call mma_deallocate(atom)
call mma_deallocate(xyz)
call mma_deallocate(force)

irc = _RC_ALL_IS_WELL_

return

end subroutine GROM
