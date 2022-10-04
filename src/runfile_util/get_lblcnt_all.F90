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

subroutine Get_LblCnt_All(xLblCnt)

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
#include "Molcas.fh"
character(len=LenIn) :: xLblCnt(*)
integer(kind=iwp) :: nAtoms, nAtoms_all
character(len=LenIn) ::  xLblCnt_Unique(MxAtom)
real(kind=wp), allocatable :: Coord(:,:)

call Get_iScalar('Unique atoms',nAtoms)
call mma_allocate(Coord,3,nAtoms,Label='Coord')
call Get_dArray('Unique Coordinates',Coord,3*nAtoms)
call Get_Name(xLblCnt_Unique)
call Get_cArray('Unique Atom Names',xLblCnt_Unique,LENIN*nAtoms)
call Get_Name_All_(Coord,nAtoms,nAtoms_all,xLblCnt_Unique,xLblCnt)
call mma_deallocate(Coord)

return

end subroutine Get_LblCnt_All
