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
!  Get_nAtoms_All
!
!> @author R. Lindh
!>
!> @details
!> Get number of all atoms (not only symmetry unique) from RUNFILE.
!>
!> @param[out] nAtoms_All Number of all atoms in the molecule
!***********************************************************************

subroutine Get_nAtoms_All(nAtoms_All)

use stdalloc, only: mma_allocate, mma_deallocate
use Symmetry_Info, only: iOper, nIrrep, Symmetry_Info_Get
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(out) :: nAtoms_All
integer(kind=iwp) :: Active = 0, iAll_Atom, iChAtom, iCoSet(0:7), iGen(3), iUnique_Atom, nAtoms, nCoSet, nGen
real(kind=wp), allocatable :: Coord(:,:)
integer(kind=iwp), external :: iChxyz

call Get_iScalar('Unique atoms',nAtoms)
call mma_allocate(Coord,3,nAtoms,label='Coord')
call Get_dArray('Unique Coordinates',Coord,3*nAtoms)

!write(u6,*) 'Enter Get_nAtoms_All'
!                                                                      *
!***********************************************************************
!                                                                      *
if (Active == 0) then
  call Symmetry_Info_Get()
  Active = 1
end if
!write(u6,*) 'Get_nAtoms_All: nIrrep',nIrrep
!write(u6,*) 'Get_nAtoms_All: iOper',(iOper(i),i=0,nIrrep-1)
!                                                                      *
!***********************************************************************
!                                                                      *
nGen = 0
if (nIrrep == 2) nGen = 1
if (nIrrep == 4) nGen = 2
if (nIrrep == 8) nGen = 3
if (nGen >= 1) iGen(1) = iOper(1)
if (nGen >= 2) iGen(2) = iOper(2)
if (nGen == 3) iGen(3) = iOper(4)
!write(u6,*) 'nGen=',nGen
!write(u6,*) 'iGen=',iGen
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute total number of centers.

iAll_Atom = 0
do iUnique_Atom=1,nAtoms
  !write(u6,*) 'iUnique_Atom=',iUnique_Atom

  iChAtom = iChxyz(Coord(:,iUnique_Atom),iGen,nGen)
  !write(u6,*) 'iChAtom=',iChAtom
  call CoSet(iCoSet,nCoSet,iChAtom)
  !write(u6,*) 'nCoSet=',nCoSet
  !write(u6,*) 'iCoSet=',(iCoSet(i),i=0,nCoset-1)

  iAll_Atom = iAll_Atom+nCoSet

end do

nAtoms_all = iAll_Atom
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_deallocate(Coord)

!write(u6,*) 'nAtoms_All=',nAtoms_All
!write(u6,*) 'Exit Get_nAtoms_All'
return

end subroutine Get_nAtoms_All
