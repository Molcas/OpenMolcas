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
! Copyright (C) 2009, Roland Lindh                                     *
!***********************************************************************

subroutine DstChk(xyz,Lbls,mCentr)
!***********************************************************************
!                                                                      *
! Object: to check that the structure is realistic.                    *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!***********************************************************************

use Constants, only: Zero, Angstrom
use Definitions, only: wp, iwp, u6

implicit none
#include "LenIn.fh"
integer(kind=iwp), intent(in) :: mCentr
real(kind=wp), intent(in) :: xyz(3,mCentr)
character(len=LenIn), intent(in) :: Lbls(mCentr)
integer(kind=iwp) :: icc, iLarge, jcc
real(kind=wp) :: R, RMax, RMin, x1, x2, y1, y2, z1, z2

if (mCentr < 5) return

iLarge = 0
do icc=1,mCentr
  if (index('1234567890',Lbls(icc)(2:2)) == 0) iLarge = 1
end do
if (iLarge == 1) return

RMax = Zero
RMin = huge(RMin)
do icc=1,mCentr
  x1 = xyz(1,icc)
  y1 = xyz(2,icc)
  z1 = xyz(3,icc)
  do jcc=1,icc-1
    x2 = xyz(1,jcc)
    y2 = xyz(2,jcc)
    z2 = xyz(3,jcc)
    R = sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
    RMin = min(RMin,R)
    RMax = max(RMax,r)
  end do
end do

if (RMax*Angstrom < 0.7_wp) then
  write(u6,*) 'All bonds shorter than 0.7 angstrom, this is probably wrong!'
  write(u6,*) 'The program will stop execution. To proceed, correct the '
  write(u6,*) 'input or use the "Expert" keyword to force execution'
  call AbEnd()
end if
if (RMin*Angstrom > 2.8_wp) then
  write(u6,*) 'All bonds longer than 2.8 angstrom, this is probably wrong!'
  write(u6,*) 'The program will stop execution. To proceed, correct the '
  write(u6,*) 'input or use the "Expert" keyword to force execution'
  call AbEnd()
end if

return

end subroutine DstChk
