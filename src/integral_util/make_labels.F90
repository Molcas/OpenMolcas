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

!#define _DEBUGPRINT_
subroutine Make_Labels(LblCbs,LblSbs,MxFnc,lMax)

use define_af, only: AngTp
use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: MxFnc, lMax
character(len=8), intent(inout) :: LblCBs(MxFnc), LblSBs(MxFnc)
integer(kind=iwp) :: i, ix, ixyz, iy, iyMax, iz, l, lxyz, m, n
character(len=3) :: Sgn

! Generate cartesian labels

lxyz = 0
do ixyz=0,lMax
  do ix=ixyz,0,-1
    iyMax = ixyz-ix
    do iy=iyMax,0,-1
      lxyz = lxyz+1
      iz = ixyz-ix-iy
      ! Form labels for cartesian basis functions
      write(LblCBs(lxyz),'(A,3I2.2)') AngTp(ixyz),ix,iy,iz
    end do
  end do
end do
if (lMax >= 0) LblCBs(1) = '01s     '
if (lMax >= 1) then
  LblCBs(2) = '02px    '
  LblCBs(3) = '02py    '
  LblCBs(4) = '02pz    '
end if

! Do the same for the spherical gaussians.

i = 0
do n=0,lMax
  do l=n,0,-2
    do m=-l,l
      if (m < 0) then
        Sgn = '-  '
      else if (m > 0) then
        Sgn = '+  '
      else
        Sgn = '   '
      end if
      i = i+1
      write(LblSbs(i),'(I2.2,A,I2.2,A)') n+1,AngTp(l),abs(m),Sgn
    end do
  end do
end do
#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) 'lMax,MxFnc=',lMax,MxFnc
write(u6,*)
write(u6,*) 'LblCBs:'
do ixyz=1,lxyz
  write(u6,'(A)') LblCBs(ixyz)
end do
write(u6,*) 'LblSBs:'
do ixyz=1,i
  write(u6,'(A)') LblCBs(ixyz)
end do
write(u6,*)
#endif

return

end subroutine Make_Labels
