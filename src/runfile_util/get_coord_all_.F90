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

subroutine Get_Coord_All_(Coord_Unique,nUnique_Atoms,Coord_All,nAll_Atoms)

use Symmetry_Info, only: iOper, nIrrep, Symmetry_Info_Get
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nUnique_Atoms, nAll_Atoms
real(kind=wp), intent(in) :: Coord_Unique(3,nUnique_Atoms)
real(kind=wp), intent(out) :: Coord_All(3,nAll_Atoms)
integer(kind=iwp) :: Active = 0, iAll_Atom, iChAtom, iCo, iCoSet(0:7,0:7), iGen(3), iStab(0:7), iUnique_Atom, MaxDCR, nCoSet, &
                     nGen, nStab
integer(kind=iwp), external :: iChxyz

!write(u6,*) 'Enter Get_Coord_All_'
!                                                                      *
!***********************************************************************
!                                                                      *
if (Active == 0) then
  call Symmetry_Info_Get()
  Active = 1
end if
!write(u6,*) 'Get_Coord_All_: nIrrep=',nIrrep
!write(u6,*) 'Get_Coord_All_: iOper=',(iOper(i),i=0,nIrrep-1)
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
! Generate list of all coordinates, index arrays, etc.

iAll_Atom = 0
MaxDCR = 0
do iUnique_Atom=1,nUnique_Atoms
  !write(u6,*) 'iUnique_Atom=',iUnique_Atom

  iChAtom = iChxyz(Coord_Unique(:,iUnique_Atom),iGen,nGen)
  !write(u6,*) 'iChAtom=',iChAtom
  call Stblz(iChAtom,nStab,iStab,MaxDCR,iCoSet)
  nCoSet = nIrrep/nStab
  !write(u6,*) 'In Get_Coord_All'
  !write(u6,*) 'nCoset=',nCoset
  !write(u6,*) 'iCoset=',(iCoset(i,0),i=0,nCoset-1)

  do iCo=0,nCoSet-1
    !write(u6,*) 'In Get_Coord_All'
    !write(u6,*) 'iCo,iCoSet(iCo,0)=',iCo,iCoSet(iCo,0)
    iAll_Atom = iAll_Atom+1
    call OA(iCoSet(iCo,0),Coord_Unique(:,iUnique_Atom),Coord_All(:,iAll_Atom))
  end do

end do

!call RecPrt('Coord_Unique',' ',Coord_Unique,3,nUnique_Atoms)
!call RecPrt('Coord_All',' ',Coord_All,3,nAll_Atoms)
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Get_Coord_All_
