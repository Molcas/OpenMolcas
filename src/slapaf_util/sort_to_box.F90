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

subroutine Sort_to_Box(Coor,nAtoms,iTab,nMax,nx,ny,nz,iBox,iANr,xmin,ymin,zmin,Box_Size)

use Definitions, only: wp, iwp
!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: nAtoms, nMax, nx, ny, nz, iANr(nAtoms)
real(kind=wp), intent(in) :: Coor(3,nAtoms), xmin, ymin, zmin, Box_Size
integer(kind=iwp), intent(out) :: iTab(0:nMax,nx,ny,nz), iBox(3,nAtoms)
integer(kind=iwp) :: iAtom, ix, iy, iz, Nr
real(kind=wp) :: x, y, z
integer(kind=iwp), external :: iTabRow

!                                                                      *
!***********************************************************************
!                                                                      *

iTab(:,:,:,:) = 0

do iAtom=1,nAtoms

  if (iTabRow(iAnr(iAtom)) == 0) cycle

  x = Coor(1,iAtom)-xmin
  ix = int(x/Box_Size)+1
  y = Coor(2,iAtom)-ymin
  iy = int(y/Box_Size)+1
  z = Coor(3,iAtom)-zmin
  iz = int(z/Box_Size)+1
  iBox(1,iAtom) = ix
  iBox(2,iAtom) = iy
  iBox(3,iAtom) = iz

  Nr = iTab(0,ix,iy,iz)+1
  if (Nr > nMax) then
    call WarningMessage(2,'Sort_to_Box: Nr > nMax')
    call Abend()
  end if
# ifdef _DEBUGPRINT_
  write(u6,*) 'Sort_to_Box: ix,iy,iz,Nr,iAtom=',ix,iy,iz,Nr,iAtom
# endif

  iTab(0,ix,iy,iz) = Nr
  iTab(Nr,ix,iy,iz) = iAtom

end do

return

end subroutine Sort_to_Box
