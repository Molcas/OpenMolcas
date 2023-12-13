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

subroutine Dihedr(Lbls,xyz,mCentr,rtrnc,Max_Center)
!***********************************************************************
!                                                                      *
! Object: to compute dihedral angles from a list of coordinates.       *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!***********************************************************************

use Constants, only: Zero, One, deg2rad
use Definitions, only: wp, iwp, u6

implicit none
#include "Molcas.fh"
integer(kind=iwp), intent(in) :: mCentr, Max_Center
character(len=LenIn), intent(in) :: Lbls(mCentr)
real(kind=wp), intent(in) :: xyz(3,mCentr), rtrnc
integer(kind=iwp) :: ic, jc, kc, lc, Lu
real(kind=wp) :: arg, Bt(3,4), Coor(3,4), Dummy(1), Phi1, Phi12, Phi2, r1, r12, r2, r23, r3, Tau, x1, x12, x2, x23, x3, x4, y1, &
                 y12, y2, y23, y3, y4, z1, z12, z2, z23, z3, z4
character(len=8) :: Label
logical(kind=iwp) :: Typ
real(kind=wp), parameter :: Thr = 1.0e-12_wp

Lu = u6
Label = ' '
if (mCentr > Max_Center) return

Typ = .false.
do ic=1,mCentr
  x2 = xyz(1,ic)
  y2 = xyz(2,ic)
  z2 = xyz(3,ic)
  Coor(:,2) = xyz(:,ic)
  do jc=1,mCentr
    if (jc == ic) cycle
    x3 = xyz(1,jc)
    y3 = xyz(2,jc)
    z3 = xyz(3,jc)
    r2 = sqrt((x3-x2)**2+(y3-y2)**2+(z3-z2)**2)
    if ((r2 > rtrnc) .or. (r2 == Zero)) cycle
    !write(Lu,*)
    Coor(:,3) = xyz(:,jc)
    do kc=1,mCentr
      if (kc == ic) cycle
      if (kc == jc) cycle
      x1 = xyz(1,kc)
      y1 = xyz(2,kc)
      z1 = xyz(3,kc)
      r1 = sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
      if ((r1 > rtrnc) .or. (r1 == Zero)) cycle
      arg = ((x1-x2)*(x3-x2)+(y1-y2)*(y3-y2)+(z1-z2)*(z3-z2))/(r1*r2)
      if (abs(arg) > One) arg = sign(One,arg)
      if (One-abs(arg) < Thr) cycle
      Phi1 = acos(arg)/deg2rad
      x12 = (y2-y1)*(z3-z2)-(y3-y2)*(z2-z1)
      y12 = (z2-z1)*(x3-x2)-(z3-z2)*(x2-x1)
      z12 = (x2-x1)*(y3-y2)-(x3-x2)*(y2-y1)
      r12 = sqrt(x12**2+y12**2+z12**2)
      if (r12 == Zero) cycle
      Coor(:,1) = xyz(:,kc)
      do lc=kc+1,mCentr
        if (lc == ic) cycle
        if (lc == jc) cycle
        if (lc == kc) cycle
        x4 = xyz(1,lc)
        y4 = xyz(2,lc)
        z4 = xyz(3,lc)
        r3 = sqrt((x4-x3)**2+(y4-y3)**2+(z4-z3)**2)
        if ((r3 > rtrnc) .or. (r3 == Zero)) cycle
        arg = ((x2-x3)*(x4-x3)+(y2-y3)*(y4-y3)+(z2-z3)*(z4-z3))/(r2*r3)
        if (abs(arg) > One) arg = sign(One,arg)
        if (One-abs(arg) < Thr) cycle
        Phi2 = acos(arg)/deg2rad
        x23 = (y3-y2)*(z4-z3)-(y4-y3)*(z3-z2)
        y23 = (z3-z2)*(x4-x3)-(z4-z3)*(x3-x2)
        z23 = (x3-x2)*(y4-y3)-(x4-x3)*(y3-y2)
        r23 = sqrt(x23**2+y23**2+z23**2)
        if (r23 == Zero) cycle
        Coor(:,4) = xyz(:,lc)
        !arg = (x12*x23+y12*y23+z12*z23)/(r12*r23)
        !if (abs(arg) > One) arg = sign(One,arg)
        !Phi12 = acos(arg)/deg2rad
        call Trsn(Coor,4,Tau,Bt,.false.,.false.,Label,Dummy,.false.)
        Phi12 = Tau/deg2rad
        if (.not. Typ) then
          Typ = .true.
          write(Lu,*)
          write(Lu,'(10X,A)') ' ***************************************************************'
          write(Lu,'(10X,A)') ' *              Valence Dihedral Angles / Degree               *'
          write(Lu,'(10X,A)') ' ***************************************************************'
          write(Lu,'(7X,A)') '             Atom centers                       Phi1     Phi2     Theta '
        end if
        write(Lu,'(10X,4(I2,1X,A,2X),1X,3(F7.2,2X))') kc,Lbls(kc),ic,Lbls(ic),jc,Lbls(jc),lc,Lbls(lc),Phi1,Phi2,Phi12
      end do
    end do
  end do
end do

return

end subroutine Dihedr
