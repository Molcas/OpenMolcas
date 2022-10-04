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

subroutine Get_Name_All_(Coord_Unique_Atoms,nUnique_Atoms,nAll_Atoms,Element_Unique,Element)

use Symmetry_Info, only: iOper, nIrrep, Symmetry_Info_Get
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: nUnique_Atoms, nAll_Atoms
real(kind=wp) :: Coord_Unique_Atoms(3,nUnique_Atoms)
character(len=*) :: Element_Unique(nUnique_Atoms), Element(*)
integer(kind=iwp) :: Active = 0, i, iAll_Atom, iChAtom, iCoSet(0:7), iGen(3), iUnique_Atom, nCoSet, nGen
integer(kind=iwp), external :: iChxyz

!                                                                      *
!***********************************************************************
!                                                                      *
if (Active == 0) then
  call Symmetry_Info_Get()
  Active = 1
end if
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
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute total number of centers.

iAll_Atom = 0
do iUnique_Atom=1,nUnique_Atoms

  iChAtom = iChxyz(Coord_Unique_Atoms(1,iUnique_Atom),iGen,nGen)
  call CoSet(iCoSet,nCoSet,iChAtom)

  do i=1,nCoSet
    iAll_Atom = iAll_Atom+1
    Element(iAll_Atom) = Element_Unique(iUnique_Atom)
  end do

end do

nAll_Atoms = iAll_Atom
!                                                                      *
!***********************************************************************
!                                                                      *
!write(u6,*) 'Exit Get_nAtoms_All_'

return

end subroutine Get_Name_All_
