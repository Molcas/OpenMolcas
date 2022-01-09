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

implicit real*8(A-H,O-Z)

#include "print.fh"
#include "real.fh"
#include "angstr.fh"
#include "Molcas.fh"
real*8 xyz(3,mCentr)
character*(LENIN) Lbls(mCentr)

lu = 6

if (mCentr < 5) Go To 99

iLarge = 0
do icc=1,mCentr
  if (index('1234567890',Lbls(icc)(2:2)) == 0) iLarge = 1
end do
if (iLarge == 1) goto 99

RMax = 0.0d0
RMin = 1.0d10
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

if (RMax*Angstr < 0.7d0) then
  write(lu,*) 'All bonds shorter than 0.7 angstrom, this is probably wrong!'
  write(lu,*) 'The program will stop execution. To proceed, correct the '
  write(lu,*) 'input or use the "Expert" keyword to force execution'
  call AbEnd()
end if
if (RMin*Angstr > 2.8d0) then
  write(lu,*) 'All bonds longer than 2.8 angstrom, this is probably wrong!'
  write(lu,*) 'The program will stop execution. To proceed, correct the '
  write(lu,*) 'input or use the "Expert" keyword to force execution'
  call AbEnd()
end if

99 continue

return

end subroutine DstChk
