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

subroutine Get_nAtoms_All_(Coord_Unique_Atoms,nUnique_Atoms,nAll_Atoms)

use Symmetry_Info, only: nIrrep, iOper, Symmetry_Info_Get

implicit real*8(a-h,o-z)
#include "real.fh"
integer iGen(3), iCoSet(0:7)
real*8 Coord_Unique_Atoms(3,nUnique_Atoms)
integer, save :: Active = 0

!write(6,*) 'Enter Get_nAtoms_All_'
!                                                                      *
!***********************************************************************
!                                                                      *
if (Active == 0) then
  call Symmetry_Info_Get()
  Active = 1
end if
!write(6,*) 'Get_nAtoms_All_: nIrrep',nIrrep
!write(6,*) 'Get_nAtoms_All_: iOper',(iOper(i),i=0,nIrrep-1)
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
!write(6,*) 'nGen=',nGen
!write(6,*) 'iGen=',iGen
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute total number of centers.

iAll_Atom = 0
do iUnique_Atom=1,nUnique_Atoms
  !write(6,*) 'iUnique_Atom=',iUnique_Atom

  iChAtom = iChxyz(Coord_Unique_Atoms(1,iUnique_Atom),iGen,nGen)
  !write(6,*) 'iChAtom=',iChAtom
  call CoSet(iCoSet,nCoSet,iChAtom)
  !write(6,*) 'nCoSet=',nCoSet
  !write(6,*) 'iCoSet=',(iCoSet(i),i=0,nCoset-1)

  iAll_Atom = iAll_Atom+nCoSet

end do

nAll_Atoms = iAll_Atom
!                                                                      *
!***********************************************************************
!                                                                      *
!write(6,*) 'Exit Get_nAtoms_All_'

return

end subroutine Get_nAtoms_All_
