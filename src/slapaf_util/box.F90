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

subroutine Box(Coor,mTtAtm,iANr,TabB,TabA,nBonds,nMax)

implicit real*8(a-h,o-z)
#include "real.fh"
#include "stdalloc.fh"
real*8 Coor(3,mTtAtm)
integer iANr(mTtAtm)
integer, allocatable :: TabB(:,:), TabA(:,:,:), iBox(:,:), Tab(:,:,:,:)

!                                                                      *
!***********************************************************************
!                                                                      *
!define _DEBUGPRINT_
!                                                                      *
!***********************************************************************
!                                                                      *
if (mTtAtm < 2) then
  write(6,*) 'Too few atoms to relax: mTtAtm=',mTtAtm
  call WarningMessage(2,'mTtAtm < 2')
  call Abend()
end if

ThrB = 0.40d0
#ifdef _DEBUGPRINT_
call RecPrt('Box: Coor',' ',Coor,3,mTtAtm)
write(6,*) 'Box: ThrB=',ThrB
#endif

xmin = 1.0D+10
ymin = 1.0D+10
zmin = 1.0D+10
xmax = -1.0D+10
ymax = -1.0D+10
zmax = -1.0D+10

! Establish boundaries

do iAtom=1,mTtAtm
  xmin = min(xmin,Coor(1,iAtom))
  xmax = max(xmax,Coor(1,iAtom))
  ymin = min(ymin,Coor(2,iAtom))
  ymax = max(ymax,Coor(2,iAtom))
  zmin = min(zmin,Coor(3,iAtom))
  zmax = max(zmax,Coor(3,iAtom))
end do
xmin = xmin-1.0D-2
xmax = xmax+1.0D-2
ymin = ymin-1.0D-2
ymax = ymax+1.0D-2
zmin = zmin-1.0D-2
zmax = zmax+1.0D-2

Box_Size = 8.0d0  ! a.u.
nx = max(1,int((xmax-xmin)/Box_Size)+1)
adjust = (dble(nx)*Box_size-(xmax-xmin))/Two
xmin = xmin-adjust
xmax = xmax+adjust
ny = max(1,int((ymax-ymin)/Box_Size)+1)
adjust = (dble(ny)*Box_size-(ymax-ymin))/Two
ymin = ymin-adjust
ymax = ymax+adjust
nz = max(1,int((zmax-zmin)/Box_Size)+1)
adjust = (dble(nz)*Box_size-(zmax-zmin))/Two
zmin = zmin-adjust
zmax = zmax+adjust
#ifdef _DEBUGPRINT_
write(6,*) 'nx,ny,nz=',nx,ny,nz
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
