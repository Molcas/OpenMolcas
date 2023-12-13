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

subroutine Get_Name_All(LblCnt)

use stdalloc, only: mma_allocate, mma_deallocate
use Symmetry_Info, only: iOper, nIrrep, Symmetry_Info_Get
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
character(len=*), intent(_OUT_) :: LblCnt(*)
#include "Molcas.fh"
integer(kind=iwp) :: Active = 0, i, iAll_Atom, iChAtom, iCoSet(0:7), iGen(3), iUnique_Atom, nAtoms, nCoSet, nGen
character(len=len(LblCnt)), allocatable :: LblCnt_Unique(:)
real(kind=wp), allocatable :: Coord(:,:)
integer(kind=iwp), external :: iChxyz

call Get_iScalar('Unique atoms',nAtoms)
call mma_allocate(Coord,3,nAtoms,label='Coord')
call mma_allocate(LblCnt_Unique,nAtoms,label='LblCnt_Unique')
call Get_dArray('Unique Coordinates',Coord,3*nAtoms)
select case (len(LblCnt))
  case (2)
    call Get_Name(LblCnt_Unique)
  case (LenIn)
    call Get_cArray('Unique Atom Names',LblCnt_Unique,LenIn*nAtoms)
  case default
    call SysAbendMsg('Get_Name_All','Wrong character length','Aborting')
end select
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
do iUnique_Atom=1,nAtoms

  iChAtom = iChxyz(Coord(:,iUnique_Atom),iGen,nGen)
  call CoSet(iCoSet,nCoSet,iChAtom)

  do i=1,nCoSet
    iAll_Atom = iAll_Atom+1
    LblCnt(iAll_Atom) = LblCnt_Unique(iUnique_Atom)
  end do

end do

!***********************************************************************
!                                                                      *
!write(u6,*) 'Exit Get_Name_All'

call mma_deallocate(Coord)
call mma_deallocate(LblCnt_Unique)

return

end subroutine Get_Name_All
