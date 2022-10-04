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

use Symmetry_Info, only: nIrrep, iOper, Symmetry_Info_Get
implicit real*8(a-h,o-z)
#include "real.fh"
integer iGen(3), iCoSet(0:7,0:7), iStab(0:7)
real*8 Coord_Unique(3,nUnique_Atoms), Coord_All(3,nAll_Atoms)
integer, save :: Active = 0

!write(6,*) 'Enter Get_Coord_All_'
!                                                                      *
!***********************************************************************
!                                                                      *
if (Active == 0) then
  call Symmetry_Info_Get()
  Active = 1
end if
!write(6,*) 'Get_Coord_All_: nIrrep=',nIrrep
!write(6,*) 'Get_Coord_All_: iOper=',(iOper(i),i=0,nIrrep-1)
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
  !write(6,*) 'iUnique_Atom=',iUnique_Atom

  iChAtom = iChxyz(Coord_Unique(1,iUnique_Atom),iGen,nGen)
  !write(6,*) 'iChAtom=',iChAtom
  call Stblz(iChAtom,nStab,iStab,MaxDCR,iCoSet)
  nCoSet = nIrrep/nStab
  !write(6,*) 'In Get_Coord_All'
  !write(6,*) 'nCoset=',nCoset
  !write(6,*) 'iCoset=',(iCoset(i,0),i=0,nCoset-1)

  do iCo=0,nCoSet-1
    !write(6,*) 'In Get_Coord_All'
    !write(6,*) 'iCo,iCoSet(iCo,0)=',iCo,iCoSet(iCo,0)
    iAll_Atom = iAll_Atom+1
    call OA(iCoSet(iCo,0),Coord_Unique(1:3,iUnique_Atom),Coord_All(1:3,iAll_Atom))
  end do

end do

!call RecPrt('Coord_Unique',' ',Coord_Unique,3,nUnique_Atoms)
!call RecPrt('Coord_All',' ',Coord_All,3,nAll_Atoms)
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Get_Coord_All_
