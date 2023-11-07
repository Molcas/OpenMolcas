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

! This subroutine should be in a module, to avoid explicit interfaces
#ifndef _IN_MODULE_
#error "This file must be compiled inside a module"
#endif

subroutine Box(Coor,mTtAtm,iANr,TabB,TabA,nBonds,nMax)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Two, Eight
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: mTtAtm, iANr(mTtAtm)
real(kind=wp), intent(in) :: Coor(3,mTtAtm)
integer(kind=iwp), intent(out) :: nBonds, nMax
integer(kind=iwp), allocatable, intent(out) :: TabB(:,:), TabA(:,:,:)
integer(kind=iwp) :: iAtom, nBondMax, nx, ny, nz
real(kind=wp) :: adjust, Box_Size, ThrB, xmax, xmin, ymax, ymin, zmax, zmin
integer(kind=iwp), allocatable :: iBox(:,:), Tab(:,:,:,:)

!                                                                      *
!***********************************************************************
!                                                                      *
!#define _DEBUGPRINT_
!                                                                      *
!***********************************************************************
!                                                                      *
if (mTtAtm < 2) then
  write(u6,*) 'Too few atoms to relax: mTtAtm=',mTtAtm
  call WarningMessage(2,'mTtAtm < 2')
  call Abend()
end if

ThrB = 0.4_wp
#ifdef _DEBUGPRINT_
call RecPrt('Box: Coor',' ',Coor,3,mTtAtm)
write(u6,*) 'Box: ThrB=',ThrB
#endif

xmin = 1.0e10_wp
ymin = 1.0e10_wp
zmin = 1.0e10_wp
xmax = -1.0e10_wp
ymax = -1.0e10_wp
zmax = -1.0e10_wp

! Establish boundaries

do iAtom=1,mTtAtm
  xmin = min(xmin,Coor(1,iAtom))
  xmax = max(xmax,Coor(1,iAtom))
  ymin = min(ymin,Coor(2,iAtom))
  ymax = max(ymax,Coor(2,iAtom))
  zmin = min(zmin,Coor(3,iAtom))
  zmax = max(zmax,Coor(3,iAtom))
end do
xmin = xmin-1.0e-2_wp
xmax = xmax+1.0e-2_wp
ymin = ymin-1.0e-2_wp
ymax = ymax+1.0e-2_wp
zmin = zmin-1.0e-2_wp
zmax = zmax+1.0e-2_wp

Box_Size = Eight  ! a.u.
nx = max(1,int((xmax-xmin)/Box_Size)+1)
adjust = (nx*Box_size-(xmax-xmin))/Two
xmin = xmin-adjust
xmax = xmax+adjust
ny = max(1,int((ymax-ymin)/Box_Size)+1)
adjust = (ny*Box_size-(ymax-ymin))/Two
ymin = ymin-adjust
ymax = ymax+adjust
nz = max(1,int((zmax-zmin)/Box_Size)+1)
adjust = (nz*Box_size-(zmax-zmin))/Two
zmin = zmin-adjust
zmax = zmax+adjust
#ifdef _DEBUGPRINT_
write(u6,*) 'nx,ny,nz=',nx,ny,nz
#endif

nMax = 100
!nf nMax = 40
! AOM Fixed this size to account for double VdW counting
nBondMax = mTtAtm*(mTtAtm+1)
! AOM
call mma_allocate(TabB,3,nBondMax,Label='TabB')
call mma_allocate(TabA,[1,2],[0,nMax],[1,mTtAtm],Label='TabA')
call mma_allocate(Tab,[0,nMax],[1,nx],[1,ny],[1,nz],Label='Tab')
call mma_allocate(iBox,3,mTtAtm,Label='iBox')

call Sort_to_Box(Coor,mTtAtm,Tab,nMax,nx,ny,nz,iBox,iANr,xmin,ymin,zmin,Box_Size)

call Find_Bonds(Coor,mTtAtm,Tab,nMax,nx,ny,nz,iBox,iANr,TabB,nBonds,nBondMax,TabA,ThrB)

call mma_deallocate(iBox)
call mma_deallocate(Tab)

return

end subroutine Box
