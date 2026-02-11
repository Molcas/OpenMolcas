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
! Copyright (C) 1993, Roland Lindh                                     *
!***********************************************************************

subroutine Angles(Lbls,xyz,mCentr,rtrnc,Max_Center)
!***********************************************************************
!                                                                      *
! Object: to compute angles from a list of coordinates.                *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!***********************************************************************

use Molcas, only: LenIn
use Constants, only: Zero, One, deg2rad
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: mCentr, Max_Center
character(len=LenIn) Lbls(mCentr)
real(kind=wp) :: xyz(3,mCentr), rtrnc
integer(kind=iwp) :: ic, jc, kc
logical(kind=iwp) :: tp
real(kind=wp) :: Arg, Phi, r1, r2, x1, x2, x3, y1, y2, y3, z1, z2, z3

if (mCentr > Max_Center) return

tp = .false.
! The center atom
do ic=1,mCentr
  x1 = xyz(1,ic)
  y1 = xyz(2,ic)
  z1 = xyz(3,ic)
  do jc=1,mCentr
    if (jc == ic) cycle
    x2 = xyz(1,jc)
    y2 = xyz(2,jc)
    z2 = xyz(3,jc)
    r1 = sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
    if (r1 > rtrnc .or. r1 == Zero) cycle
    do kc=jc+1,mCentr
      if (kc == ic) cycle
      x3 = xyz(1,kc)
      y3 = xyz(2,kc)
      z3 = xyz(3,kc)
      r2 = sqrt((x3-x1)**2+(y3-y1)**2+(z3-z1)**2)
      if (r2 > rtrnc .or. r2 == Zero) cycle
      arg = ((x2-x1)*(x3-x1)+(y2-y1)*(y3-y1)+(z2-z1)*(z3-z1))/(r1*r2)
      if (abs(arg) > One) arg = sign(One,arg)
      Phi = acos(arg)/deg2rad
      if (.not. tp) then
        tp = .true.
        write(u6,*)
        write(u6,'(19X,A)') ' ************************************** '
        write(u6,'(19X,A)') ' *    Valence Bond Angles / degree    * '
        write(u6,'(19X,A)') ' ************************************** '
        write(u6,'(19X,A)') '       Atom centers                 Phi'
      end if
      write(u6,'(21X,3(I2,1X,A,2X),1X,F6.2)') jc,Lbls(jc),ic,Lbls(ic),kc,Lbls(kc),Phi
    end do
  end do
end do

return

end subroutine Angles
